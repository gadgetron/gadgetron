/*****************************************
*  Standalone ISMRMRD Gadgetron Client  
*
* Author: Michael S. Hansen
* 
* Dependencies: ISMRMRD and Boost
*
*****************************************/

//TODO:
// -Blobs (for DICOM image support)
//  - First implementation is in, but testing needed
// -NIFTI and Analyze output
// -Check on potential threading problem with asio socket 
//    - having and reading and writing thread is supposedly not safe, but seems to work here
// -Static linking for standalone executable. 

#include <boost/program_options.hpp>
#include <boost/asio.hpp>

#include <ismrmrd/ismrmrd.h>
#include <ismrmrd/dataset.h>
#include <ismrmrd/meta.h>
#include <ismrmrd/xml.h>
#include <ismrmrd/waveform.h>

#include <streambuf>
#include <time.h>
#include <iomanip>
#include <iostream>
#include <exception>
#include <map>
#include <thread>
#include <chrono>
#include <condition_variable>

#include "NHLBICompression.h"
#include "GadgetronTimer.h"

#if defined GADGETRON_COMPRESSION_ZFP
#include <zfp.h>
#endif

std::string get_date_time_string()
{
    time_t rawtime;
    struct tm * timeinfo;
    time ( &rawtime );
    timeinfo = localtime ( &rawtime );

    std::stringstream str;
    str << timeinfo->tm_year+1900 << "-"
        << std::setw(2) << std::setfill('0') << timeinfo->tm_mon+1 << "-"
        << std::setw(2) << std::setfill('0') << timeinfo->tm_mday << " "
        << std::setw(2) << std::setfill('0') << timeinfo->tm_hour << ":"
        << std::setw(2) << std::setfill('0') << timeinfo->tm_min << ":"
        << std::setw(2) << std::setfill('0') << timeinfo->tm_sec;

    std::string ret = str.str();

    return ret;
}

struct NoiseStatistics
{
    bool status;
    uint16_t channels;
    float sigma_min;
    float sigma_max;
    float sigma_mean;
    float noise_dwell_time_us;
};

#if defined GADGETRON_COMPRESSION_ZFP
size_t compress_zfp_tolerance(float* in, size_t samples, size_t coils, double tolerance, char* buffer, size_t buf_size)
{
    zfp_type type = zfp_type_float;
    zfp_field* field = NULL;
    zfp_stream* zfp = NULL;
    bitstream* stream = NULL;
    size_t zfpsize = 0;

    zfp = zfp_stream_open(NULL);
    field = zfp_field_alloc();

    zfp_field_set_pointer(field, in);

    zfp_field_set_type(field, type);
    zfp_field_set_size_2d(field, samples, coils);
    zfp_stream_set_accuracy(zfp, tolerance);

    if (zfp_stream_maximum_size(zfp, field) > buf_size) {
        zfp_field_free(field);
        zfp_stream_close(zfp);
        stream_close(stream);
        throw std::runtime_error("Insufficient buffer space for compression");
    }

    stream = stream_open(buffer, buf_size);
    if (!stream) {
        zfp_field_free(field);
        zfp_stream_close(zfp);
        stream_close(stream);
        throw std::runtime_error("Cannot open compressed stream");
    }
    zfp_stream_set_bit_stream(zfp, stream);

    if (!zfp_write_header(zfp, field, ZFP_HEADER_FULL)) {
        zfp_field_free(field);
        zfp_stream_close(zfp);
        stream_close(stream);
        throw std::runtime_error("Unable to write compression header to stream");
    }

    zfpsize = zfp_compress(zfp, field);
    if (zfpsize == 0) {
        zfp_field_free(field);
        zfp_stream_close(zfp);
        stream_close(stream);
        throw std::runtime_error("Compression failed");
    }

    zfp_field_free(field);
    zfp_stream_close(zfp);
    stream_close(stream);
    return zfpsize;
}


size_t compress_zfp_precision(float* in, size_t samples, size_t coils, unsigned int precision, char* buffer, size_t buf_size)
{
  zfp_type type = zfp_type_float;
  zfp_field* field = NULL;
  zfp_stream* zfp = NULL;
  bitstream* stream = NULL;
  size_t zfpsize = 0;

  zfp = zfp_stream_open(NULL);
  field = zfp_field_alloc();

  zfp_field_set_pointer(field, in);

  zfp_field_set_type(field, type);
  zfp_field_set_size_2d(field, samples, coils);
  zfp_stream_set_precision(zfp, precision);

  if (zfp_stream_maximum_size(zfp, field) > buf_size) {
      zfp_field_free(field);
      zfp_stream_close(zfp);
      stream_close(stream);
      throw std::runtime_error("Insufficient buffer space for compression");
  }

  stream = stream_open(buffer, buf_size);
  if (!stream) {
      zfp_field_free(field);
      zfp_stream_close(zfp);
      stream_close(stream);
      throw std::runtime_error("Cannot open compressed stream");
  }
  zfp_stream_set_bit_stream(zfp, stream);

  if (!zfp_write_header(zfp, field, ZFP_HEADER_FULL)) {
      zfp_field_free(field);
      zfp_stream_close(zfp);
      stream_close(stream);
      throw std::runtime_error("Unable to write compression header to stream");
  }

  zfpsize = zfp_compress(zfp, field);
  if (zfpsize == 0) {
      zfp_field_free(field);
      zfp_stream_close(zfp);
      stream_close(stream);
      throw std::runtime_error("Compression failed");
  }
  
  zfp_field_free(field);
  zfp_stream_close(zfp);
  stream_close(stream);
  return zfpsize;
}
#endif //GADGETRON_COMPRESSION_ZFP

namespace po = boost::program_options;
using boost::asio::ip::tcp;


enum GadgetronMessageID {
    GADGET_MESSAGE_INT_ID_MIN                             =   0,
    GADGET_MESSAGE_CONFIG_FILE                            =   1,
    GADGET_MESSAGE_CONFIG_SCRIPT                          =   2,
    GADGET_MESSAGE_PARAMETER_SCRIPT                       =   3,
    GADGET_MESSAGE_CLOSE                                  =   4,
    GADGET_MESSAGE_TEXT                                   =   5,
    GADGET_MESSAGE_QUERY                                  =   6,
    GADGET_MESSAGE_RESPONSE                               =   7,
    GADGET_MESSAGE_INT_ID_MAX                             = 999,
    GADGET_MESSAGE_EXT_ID_MIN                             = 1000,
    GADGET_MESSAGE_ACQUISITION                            = 1001, /**< DEPRECATED */
    GADGET_MESSAGE_NEW_MEASUREMENT                        = 1002, /**< DEPRECATED */
    GADGET_MESSAGE_END_OF_SCAN                            = 1003, /**< DEPRECATED */
    GADGET_MESSAGE_IMAGE_CPLX_FLOAT                       = 1004, /**< DEPRECATED */
    GADGET_MESSAGE_IMAGE_REAL_FLOAT                       = 1005, /**< DEPRECATED */
    GADGET_MESSAGE_IMAGE_REAL_USHORT                      = 1006, /**< DEPRECATED */
    GADGET_MESSAGE_EMPTY                                  = 1007, /**< DEPRECATED */
    GADGET_MESSAGE_ISMRMRD_ACQUISITION                    = 1008,
    GADGET_MESSAGE_ISMRMRD_IMAGE_CPLX_FLOAT               = 1009, /**< DEPRECATED */
    GADGET_MESSAGE_ISMRMRD_IMAGE_REAL_FLOAT               = 1010, /**< DEPRECATED */
    GADGET_MESSAGE_ISMRMRD_IMAGE_REAL_USHORT              = 1011, /**< DEPRECATED */
    GADGET_MESSAGE_DICOM                                  = 1012, /**< DEPRECATED */
    GADGET_MESSAGE_CLOUD_JOB                              = 1013,
    GADGET_MESSAGE_GADGETCLOUD_JOB                        = 1014,
    GADGET_MESSAGE_ISMRMRD_IMAGEWITHATTRIB_CPLX_FLOAT     = 1015, /**< DEPRECATED */
    GADGET_MESSAGE_ISMRMRD_IMAGEWITHATTRIB_REAL_FLOAT     = 1016, /**< DEPRECATED */
    GADGET_MESSAGE_ISMRMRD_IMAGEWITHATTRIB_REAL_USHORT    = 1017, /**< DEPRECATED */
    GADGET_MESSAGE_DICOM_WITHNAME                         = 1018,
    GADGET_MESSAGE_DEPENDENCY_QUERY                       = 1019,
    GADGET_MESSAGE_ISMRMRD_IMAGE_REAL_SHORT               = 1020, /**< DEPRECATED */
    GADGET_MESSAGE_ISMRMRD_IMAGEWITHATTRIB_REAL_SHORT     = 1021, /**< DEPRECATED */
    GADGET_MESSAGE_ISMRMRD_IMAGE                          = 1022,
    GADGET_MESSAGE_RECONDATA                              = 1023,
    GADGET_MESSAGE_ISMRMRD_WAVEFORM                       = 1026,
    GADGET_MESSAGE_EXT_ID_MAX                             = 4096
};

std::mutex mtx;

struct GadgetMessageIdentifier
{
    uint16_t id;
};

struct GadgetMessageConfigurationFile
{
    char configuration_file[1024];
};

struct GadgetMessageScript
{
    uint32_t script_length;
};

class GadgetronClientException : public std::exception
{

public:
    GadgetronClientException(std::string msg)
        : msg_(msg)
    {

    }

    virtual ~GadgetronClientException() throw() {}

    virtual const char* what() const throw()
    {
        return msg_.c_str();
    }

protected:
    std::string msg_;
};

class GadgetronClientMessageReader
{
public:
    virtual ~GadgetronClientMessageReader() {}

    /**
    Function must be implemented to read a specific message.
    */
    virtual void read(tcp::socket* s) = 0;

};

class GadgetronClientResponseReader : public GadgetronClientMessageReader
{
    void read(tcp::socket *stream) override {

        uint64_t correlation_id = 0;
        uint64_t response_length = 0;

        boost::asio::read(*stream, boost::asio::buffer(&correlation_id, sizeof(correlation_id)));
        boost::asio::read(*stream, boost::asio::buffer(&response_length, sizeof(response_length)));

        std::vector<char> response(response_length + 1,0);

        boost::asio::read(*stream, boost::asio::buffer(response.data(),response_length));

        std::cout << response.data() << std::endl;
    }
};

class GadgetronClientTextReader : public GadgetronClientMessageReader
{
  
public:
  GadgetronClientTextReader()
  {
    
  }

  virtual ~GadgetronClientTextReader()
  {
    
  }

  virtual void read(tcp::socket* stream)
  {
    size_t recv_count = 0;
    
    typedef unsigned long long size_t_type;
    uint32_t len(0);
    boost::asio::read(*stream, boost::asio::buffer(&len, sizeof(uint32_t)));

    
    char* buf = NULL;
    try {
      buf = new char[len+1];
      memset(buf, '\0', len+1);
    } catch (std::runtime_error &err) {
      std::cerr << "TextReader, failed to allocate buffer" << std::endl;
      throw;
    }
       
    if (boost::asio::read(*stream, boost::asio::buffer(buf, len)) != len)
    {
      delete [] buf;
      throw GadgetronClientException("Incorrect number of bytes read for dependency query");  
    }

    std::string s(buf);
    std::cout << s;
    delete[] buf;
  }  
};


class GadgetronClientDependencyQueryReader : public GadgetronClientMessageReader
{
  
public:
  GadgetronClientDependencyQueryReader(std::string filename) : number_of_calls_(0) , filename_(filename)
  {
    
  }

  virtual ~GadgetronClientDependencyQueryReader()
  {
    
  }

  virtual void read(tcp::socket* stream)
  {
    size_t recv_count = 0;
    
    typedef unsigned long long size_t_type;
    size_t_type len(0);
    boost::asio::read(*stream, boost::asio::buffer(&len, sizeof(size_t_type)));

    char* buf = NULL;
    try {
      buf = new char[len];
      memset(buf, '\0', len);
      memcpy(buf, &len, sizeof(size_t_type));
    } catch (std::runtime_error &err) {
      std::cerr << "DependencyQueryReader, failed to allocate buffer" << std::endl;
      throw;
    }
    
    
    if (boost::asio::read(*stream, boost::asio::buffer(buf, len)) != len)
    {
      delete [] buf;
      throw GadgetronClientException("Incorrect number of bytes read for dependency query");  
    }
    
    std::ofstream outfile;
    outfile.open (filename_.c_str(), std::ios::out|std::ios::binary);

    if (outfile.good()) {
      outfile.write(buf, len);
      outfile.close();
      number_of_calls_++;
    } else {
      delete[] buf;
      throw GadgetronClientException("Unable to write dependency query to file");  
    }
    
    delete[] buf;
  }
  
  protected:
    size_t number_of_calls_;
    std::string filename_;
};


class GadgetronClientImageMessageReader : public GadgetronClientMessageReader
{

public:
    GadgetronClientImageMessageReader(std::string filename, std::string groupname)
        : file_name_(filename)
        , group_name_(groupname)
    {

    }

    ~GadgetronClientImageMessageReader() {
    } 

    template <typename T> 
    void read_data_attrib(tcp::socket* stream, const ISMRMRD::ImageHeader& h, ISMRMRD::Image<T>& im)
    {
        im.setHead(h);

        typedef unsigned long long size_t_type;

        //Read meta attributes
        size_t_type meta_attrib_length;
        boost::asio::read(*stream, boost::asio::buffer(&meta_attrib_length, sizeof(size_t_type)));

        if (meta_attrib_length>0)
        {
            std::string meta_attrib(meta_attrib_length, 0);
            boost::asio::read(*stream, boost::asio::buffer(const_cast<char*>(meta_attrib.c_str()), meta_attrib_length));
            im.setAttributeString(meta_attrib);
        }

        //Read image data
        boost::asio::read(*stream, boost::asio::buffer(im.getDataPtr(), im.getDataSize()));
        {
            if (!dataset_) {

                {
                    mtx.lock();
                    dataset_ = std::shared_ptr<ISMRMRD::Dataset>(new ISMRMRD::Dataset(file_name_.c_str(), group_name_.c_str(), true)); // create if necessary
                    mtx.unlock();
                }
            }

            std::stringstream st1;
            st1 << "image_" << h.image_series_index;
            std::string image_varname = st1.str();

            {
                mtx.lock();
                //TODO should this be wrapped in a try/catch?
                dataset_->appendImage(image_varname, im);
                mtx.unlock();
            }
        }
    }

    virtual void read(tcp::socket* stream) 
    {
        //Read the image headerfrom the socket
        ISMRMRD::ImageHeader h;
        boost::asio::read(*stream, boost::asio::buffer(&h,sizeof(ISMRMRD::ImageHeader)));

        if (h.data_type == ISMRMRD::ISMRMRD_USHORT)
        {
            ISMRMRD::Image<unsigned short> im;
            this->read_data_attrib(stream, h, im);
        }
        else if (h.data_type == ISMRMRD::ISMRMRD_SHORT)
        {
            ISMRMRD::Image<short> im;
            this->read_data_attrib(stream, h, im);
        }
        else if (h.data_type == ISMRMRD::ISMRMRD_UINT)
        {
            ISMRMRD::Image<unsigned int> im;
            this->read_data_attrib(stream, h, im);
        }
        else if (h.data_type == ISMRMRD::ISMRMRD_INT)
        {
            ISMRMRD::Image<int> im;
            this->read_data_attrib(stream, h, im);
        }
        else if (h.data_type == ISMRMRD::ISMRMRD_FLOAT)
        {
            ISMRMRD::Image<float> im;
            this->read_data_attrib(stream, h, im);
        }
        else if (h.data_type == ISMRMRD::ISMRMRD_DOUBLE)
        {
            ISMRMRD::Image<double> im;
            this->read_data_attrib(stream, h, im);
        }
        else if (h.data_type == ISMRMRD::ISMRMRD_CXFLOAT)
        {
            ISMRMRD::Image< std::complex<float> > im;
            this->read_data_attrib(stream, h, im);
        }
        else if (h.data_type == ISMRMRD::ISMRMRD_CXDOUBLE)
        {
            ISMRMRD::Image< std::complex<double> > im;
            this->read_data_attrib(stream, h, im);
        }
        else
        {
            throw GadgetronClientException("Invalide image data type ... ");
        }
    }

protected:
    std::string group_name_;
    std::string file_name_;
    std::shared_ptr<ISMRMRD::Dataset> dataset_;
};

// ----------------------------------------------------------------
// for the analyze image format
#ifdef DT_UNKNOWN
    #undef DT_UNKNOWN
#endif // DT_UNKNOWN

enum AnalyzeDataType
{
    DT_ANA_UNKNOWN=0,

    DT_NONE                    =0,
    DT_UNKNOWN                 =0,     /* what it says, dude           */
    DT_BINARY                  =1,     /* binary (1 bit/voxel)         */
    DT_UNSIGNED_CHAR           =2,     /* unsigned char (8 bits/voxel) */
    DT_SIGNED_SHORT            =4,     /* signed short (16 bits/voxel) */
    DT_UNSIGNED_SHORT          =5,
    DT_SIGNED_INT              =8,     /* signed int (32 bits/voxel)   */
    DT_UNSIGNED_INT            =9,
    DT_FLOAT                  =16,     /* float (32 bits/voxel)        */
    DT_COMPLEX                =32,     /* complex (64 bits/voxel)      */
    DT_DOUBLE                 =64,     /* double (64 bits/voxel)       */
    DT_RGB                   =128,     /* RGB triple (24 bits/voxel)   */
    DT_ALL                   =255,     /* not very useful (?)          */

                                /*----- another set of names for the same ---*/
    DT_UINT8                   =2,
    DT_INT16                   =4,
    DT_INT32                   =8,
    DT_FLOAT32                =16,
    DT_COMPLEX64              =32,
    DT_FLOAT64                =64,
    DT_RGB24                 =128,

                                /*------------------- new codes for NIFTI ---*/
    DT_INT8                  =256,     /* signed char (8 bits)         */
    DT_UINT16                =512,     /* unsigned short (16 bits)     */
    DT_UINT32                =768,     /* unsigned int (32 bits)       */
    DT_INT64                =1024,     /* long long (64 bits)          */
    DT_UINT64               =1280,     /* unsigned long long (64 bits) */
    DT_FLOAT128             =1536,     /* long double (128 bits)       */
    DT_COMPLEX128           =1792,     /* double pair (128 bits)       */
    DT_COMPLEX256           =2048,     /* long double pair (256 bits)  */
    DT_RGBA32               =2304,     /* 4 byte RGBA (32 bits/voxel)  */
};

AnalyzeDataType getDataTypeFromRTTI(const std::string& name)
{
    AnalyzeDataType analyzeDT = DT_ANA_UNKNOWN;

    if ( name == typeid(unsigned char).name() )
    {
        analyzeDT = DT_UNSIGNED_CHAR;
    }

    if ( name == typeid(short).name() )
    {
        analyzeDT = DT_SIGNED_SHORT;
    }

    if ( name == typeid(unsigned short).name() )
    {
        analyzeDT = DT_UINT16;
    }

    if ( name == typeid(int).name() )
    {
        analyzeDT = DT_SIGNED_INT;
    }

    if ( name == typeid(unsigned int).name() )
    {
        analyzeDT = DT_UINT32;
    }

    if ( name == typeid(float).name() )
    {
        analyzeDT = DT_FLOAT;
    }

    if ( name == typeid(double).name() )
    {
        analyzeDT = DT_DOUBLE;
    }

    if ( name == typeid(long double).name() )
    {
        analyzeDT = DT_FLOAT128;
    }

    if ( name == typeid(std::complex<float>).name() )
    {
        analyzeDT = DT_COMPLEX;
    }

    if ( name == typeid(std::complex<double>).name() )
    {
        analyzeDT = DT_COMPLEX128;
    }

    if ( name == typeid(std::complex<long double>).name() )
    {
        analyzeDT = DT_COMPLEX256;
    }

    return analyzeDT;
}

struct header_key
{
    int sizeof_hdr;
    char data_type[10];
    char db_name[18];
    int extents;
    short int session_error;
    char regular;
    char hkey_un0;
};

struct image_dimension
{
    short int dim[8];
    short int unused8;
    short int unused9;
    short int unused10;
    short int unused11;
    short int unused12;
    short int unused13;
    short int unused14;
    short int datatype;
    short int bitpix;
    short int dim_un0;
    float pixdim[8];
    float vox_offset;
    float funused1;
    float funused2;
    float funused3;
    float cal_max;
    float cal_min;
    float compressed;
    float verified;
    int glmax,glmin;
};

struct data_history
{
    char descrip[80];
    char aux_file[24];
    char orient;
    char originator[10];
    char generated[10];
    char scannum[10];
    char patient_id[10];
    char exp_date[10];
    char exp_time[10];
    char hist_un0[3];
    int views;
    int vols_added;
    int start_field;
    int field_skip;
    int omax, omin;
    int smax, smin;
};

// Analyze75 header has 348 bytes
struct dsr
{
    struct header_key hk;
    struct image_dimension dime;
    struct data_history hist;
};

class IOAnalyze
{
public:

    typedef dsr HeaderType;

    IOAnalyze() {}
    virtual ~IOAnalyze() {}

    template <typename T> void array2Header(const std::vector<size_t>& dim, const std::vector<float>& pixelSize, HeaderType& header)
    {
        try
        {
            // set everything to zero
            memset(&header, 0, sizeof(dsr));

            // header_key
            header.hk.sizeof_hdr = 348;
            size_t i;
            for (i=0; i<10; i++ ) header.hk.data_type[i] = 0;
            for (i=0; i<18; i++ ) header.hk.db_name[i] = 0;
            header.hk.extents = 16384;
            header.hk.session_error = 0;
            header.hk.regular = 'r';
            header.hk.hkey_un0 = 0;

            // image_dimension
            size_t NDim = dim.size();

            header.dime.dim[0] = (short)(NDim);
            header.dime.dim[1] = (short)(dim[0]);

            if ( NDim > 1 )
                header.dime.dim[2] = (short)(dim[1]);
            else
                header.dime.dim[2] = 1;

            if ( NDim > 2 )
                header.dime.dim[3] = (short)(dim[2]);
            else
                header.dime.dim[3] = 1;

            if ( NDim > 3 )
                header.dime.dim[4] = (short)(dim[3]);
            else
                header.dime.dim[4] = 1;

            if ( NDim > 4 )
                header.dime.dim[5] = (short)(dim[4]);
            else
                header.dime.dim[5] = 1;

            if ( NDim > 5 )
                header.dime.dim[6] = (short)(dim[5]);
            else
                header.dime.dim[6] = 1;

            if ( NDim > 6 )
                header.dime.dim[7] = (short)(dim[6]);
            else
                header.dime.dim[7] = 1;

            if ( NDim > 7 )
                header.dime.unused8 = (short)(dim[7]);
            else
                header.dime.unused8 = 1;

            if ( NDim > 8 )
                header.dime.unused9 = (short)(dim[8]);
            else
                header.dime.unused9 = 1;

            if ( NDim > 9 )
                header.dime.unused10 = (short)(dim[9]);
            else
                header.dime.unused10 = 1;

            header.dime.unused11 = 0;
            header.dime.unused12 = 0;
            header.dime.unused13 = 0;
            header.dime.unused14 = 0;

            std::string rttiID = std::string(typeid(T).name());
            header.dime.datatype = (short)getDataTypeFromRTTI(rttiID);
            header.dime.bitpix = (short)(8*sizeof(T));
            header.dime.dim_un0 = 0;

            // since the NDArray does not carry the pixel spacing
            header.dime.pixdim[0] = 0;
            if ( pixelSize.size() > 0 )
                header.dime.pixdim[1] = pixelSize[0];
            if ( pixelSize.size() > 1 )
                header.dime.pixdim[2] = pixelSize[1];
            if ( pixelSize.size() > 2 )
                header.dime.pixdim[3] = pixelSize[2];
            if ( pixelSize.size() > 3 )
                header.dime.pixdim[4] = pixelSize[3];
            if ( pixelSize.size() > 4 )
                header.dime.pixdim[5] = pixelSize[4];
            if ( pixelSize.size() > 5 )
                header.dime.pixdim[6] = pixelSize[5];
            if ( pixelSize.size() > 6 )
                header.dime.pixdim[7] = pixelSize[6];

            header.dime.vox_offset = 0;
            header.dime.funused1 = 0;
            header.dime.funused2 = 0;
            header.dime.funused3 = 0;
            header.dime.cal_max = 0;
            header.dime.cal_min = 0;
            header.dime.compressed = 0;
            header.dime.verified = 0;
            header.dime.glmax = 0;
            header.dime.glmin = 0;

            // data history
            for (i=0; i<80; i++ ) header.hist.descrip[i] = 0;
            for (i=0; i<24; i++ ) header.hist.aux_file[i] = 0;
            header.hist.orient = 0;
            for (i=0; i<10; i++ ) header.hist.originator[i] = 0;
            for (i=0; i<10; i++ ) header.hist.generated[i] = 0;
            for (i=0; i<10; i++ ) header.hist.scannum[i] = 0;
            for (i=0; i<10; i++ ) header.hist.patient_id[i] = 0;
            for (i=0; i<10; i++ ) header.hist.exp_date[i] = 0;
            for (i=0; i<10; i++ ) header.hist.exp_time[i] = 0;
            for (i=0; i<3; i++ ) header.hist.hist_un0[i] = 0;
            header.hist.views = 0;
            header.hist.vols_added = 0;
            header.hist.start_field = 0;
            header.hist.field_skip = 0;
            header.hist.omax = 0;
            header.hist.omin = 0;
            header.hist.smax = 0;
            header.hist.smin = 0;
        }
        catch(...)
        {
            throw GadgetronClientException("Errors in IOAnalyze::array2Analyze(dim, header) ... ");
        }
    }
};

class GadgetronClientAnalyzeImageMessageReader : public GadgetronClientMessageReader
{

public:

    GadgetronClientAnalyzeImageMessageReader(const std::string& prefix = std::string("Image")) : prefix_(prefix)
    {

    }

    ~GadgetronClientAnalyzeImageMessageReader() {
    } 

    template <typename T>
    void read_data_attrib(tcp::socket* stream, const ISMRMRD::ImageHeader& h, ISMRMRD::Image<T>& im)
    {
        im.setHead(h);

        std::cout << "Receiving image : " << h.image_series_index << " - " << h.image_index << std::endl;

        typedef unsigned long long size_t_type;

        std::ostringstream ostr;

        if (!prefix_.empty())
        {
            ostr << prefix_ << "_";
        }

        ostr << "SLC" << h.slice << "_"
            << "CON" << h.contrast << "_"
            << "PHS" << h.phase << "_"
            << "REP" << h.repetition << "_"
            << "SET" << h.set << "_"
            << "AVE" << h.average << "_"
            << h.image_index
            << "_" << h.image_series_index;

        std::string filename = ostr.str();

        //Read meta attributes
        size_t_type meta_attrib_length;
        boost::asio::read(*stream, boost::asio::buffer(&meta_attrib_length, sizeof(size_t_type)));

        if (meta_attrib_length > 0)
        {
            std::string meta_attrib(meta_attrib_length, 0);
            boost::asio::read(*stream, boost::asio::buffer(const_cast<char*>(meta_attrib.c_str()), meta_attrib_length));

            // deserialize the meta attribute
            ISMRMRD::MetaContainer imgAttrib;
            ISMRMRD::deserialize(meta_attrib.c_str(), imgAttrib);

            std::stringstream st3;
            st3 << filename << ".attrib";
            std::string meta_varname = st3.str();

            std::ofstream outfile;
            outfile.open(meta_varname.c_str(), std::ios::out | std::ios::binary);
            outfile.write(meta_attrib.c_str(), meta_attrib_length);
            outfile.close();
        }

        //Read data
        boost::asio::read(*stream, boost::asio::buffer(im.getDataPtr(), im.getDataSize()));

        // analyze header
        std::stringstream st1;
        st1 << filename << ".hdr";
        std::string head_varname = st1.str();

        std::vector<size_t> dim(4, 1);
        dim[0] = h.matrix_size[0];
        dim[1] = h.matrix_size[1];
        dim[2] = h.matrix_size[2];
        dim[3] = h.channels;

        std::vector<float> pixelSize(4, 1);
        pixelSize[0] = h.field_of_view[0] / h.matrix_size[0];
        pixelSize[1] = h.field_of_view[1] / h.matrix_size[1];
        pixelSize[2] = h.field_of_view[2] / h.matrix_size[2];

        IOAnalyze hdr;
        dsr header;
        hdr.array2Header<T>(dim, pixelSize, header);

        std::ofstream outfileHeader;
        outfileHeader.open(head_varname.c_str(), std::ios::out | std::ios::binary);
        outfileHeader.write(reinterpret_cast<const char*>(&header), sizeof(dsr));
        outfileHeader.close();

        // data
        std::stringstream st2;
        st2 << filename << ".img";
        std::string img_varname = st2.str();

        std::ofstream outfileData;
        outfileData.open(img_varname.c_str(), std::ios::out | std::ios::binary);
        outfileData.write(reinterpret_cast<const char*>(im.getDataPtr()), sizeof(T)*dim[0] * dim[1] * dim[2] * dim[3]);
        outfileData.close();
    }

    virtual void read(tcp::socket* stream) 
    {
        //Read the image headerfrom the socket
        ISMRMRD::ImageHeader h;
        boost::asio::read(*stream, boost::asio::buffer(&h,sizeof(ISMRMRD::ImageHeader)));

        if (h.data_type == ISMRMRD::ISMRMRD_USHORT)
        {
            ISMRMRD::Image<unsigned short> im;
            this->read_data_attrib(stream, h, im);
        }
        else if (h.data_type == ISMRMRD::ISMRMRD_SHORT)
        {
            ISMRMRD::Image<short> im;
            this->read_data_attrib(stream, h, im);
        }
        else if (h.data_type == ISMRMRD::ISMRMRD_UINT)
        {
            ISMRMRD::Image<unsigned int> im;
            this->read_data_attrib(stream, h, im);
        }
        else if (h.data_type == ISMRMRD::ISMRMRD_INT)
        {
            ISMRMRD::Image<int> im;
            this->read_data_attrib(stream, h, im);
        }
        else if (h.data_type == ISMRMRD::ISMRMRD_FLOAT)
        {
            ISMRMRD::Image<float> im;
            this->read_data_attrib(stream, h, im);
        }
        else if (h.data_type == ISMRMRD::ISMRMRD_DOUBLE)
        {
            ISMRMRD::Image<double> im;
            this->read_data_attrib(stream, h, im);
        }
        else if (h.data_type == ISMRMRD::ISMRMRD_CXFLOAT)
        {
            ISMRMRD::Image< std::complex<float> > im;
            this->read_data_attrib(stream, h, im);
        }
        else if (h.data_type == ISMRMRD::ISMRMRD_CXDOUBLE)
        {
            ISMRMRD::Image< std::complex<double> > im;
            this->read_data_attrib(stream, h, im);
        }
        else
        {
            throw GadgetronClientException("Invalide image data type ... ");
        }
    }

protected:

    std::string prefix_;
};

// ----------------------------------------------------------------

#define MAX_BLOBS_LOG_10    6

class GadgetronClientBlobMessageReader 
    : public GadgetronClientMessageReader
{

public:
    GadgetronClientBlobMessageReader(std::string fileprefix, std::string filesuffix)
        : number_of_calls_(0)
        , file_prefix(fileprefix)
        , file_suffix(filesuffix)

    {

    }

    virtual ~GadgetronClientBlobMessageReader() {}

    virtual void read(tcp::socket* socket) 
    {

        // MUST READ 32-bits
        uint32_t nbytes;
        boost::asio::read(*socket, boost::asio::buffer(&nbytes,sizeof(uint32_t)));

        std::vector<char> data(nbytes,0);
        boost::asio::read(*socket, boost::asio::buffer(&data[0],nbytes));

        unsigned long long fileNameLen;
        boost::asio::read(*socket, boost::asio::buffer(&fileNameLen,sizeof(unsigned long long)));

        std::string filenameBuf(fileNameLen,0);
        boost::asio::read(*socket, boost::asio::buffer(const_cast<char*>(filenameBuf.c_str()),fileNameLen));

        typedef unsigned long long size_t_type;

        size_t_type meta_attrib_length;
        boost::asio::read(*socket, boost::asio::buffer(&meta_attrib_length, sizeof(size_t_type)));

        std::string meta_attrib;
        if (meta_attrib_length > 0)
        {
            std::string meta_attrib_socket(meta_attrib_length, 0);
            boost::asio::read(*socket, boost::asio::buffer(const_cast<char*>(meta_attrib_socket.c_str()), meta_attrib_length));
            meta_attrib = meta_attrib_socket;
        }

        std::stringstream filename;
        std::string filename_attrib;

        // Create the filename: (prefix_%06.suffix)
        filename << file_prefix << "_";
        filename << std::setfill('0') << std::setw(MAX_BLOBS_LOG_10) << number_of_calls_;
        filename_attrib = filename.str();
        filename << "." << file_suffix;
        filename_attrib.append("_attrib.xml");

        std::cout << "Writing image " << filename.str() << std::endl;

        std::ofstream outfile;
        outfile.open(filename.str().c_str(), std::ios::out | std::ios::binary);

        std::ofstream outfile_attrib;
        if (meta_attrib_length > 0)
        {
            outfile_attrib.open(filename_attrib.c_str(), std::ios::out | std::ios::binary);
        }

        if (outfile.good())
        {
            /* write 'size' bytes starting at 'data's pointer */
            outfile.write(&data[0], nbytes);
            outfile.close();

            if (meta_attrib_length > 0)
            {
                outfile_attrib.write(meta_attrib.c_str(), meta_attrib.length());
                outfile_attrib.close();
            }

            number_of_calls_++;
        }
        else
        {
            throw GadgetronClientException("Unable to write blob to output file\n");
        }
    }

protected:
    size_t number_of_calls_;
    std::string file_prefix;
    std::string file_suffix;

};

class GadgetronClientConnector
{

public:
    GadgetronClientConnector() 
        : socket_(0)
        , timeout_ms_(10000)
        , uncompressed_bytes_sent_(0)
        , compressed_bytes_sent_(0)
        , header_bytes_sent_(0)
    {

    }

    virtual ~GadgetronClientConnector() 
    {
        if (socket_) {
            socket_->close();
            delete socket_;
        }
    }

    double compression_ratio()
    {
        if (compressed_bytes_sent_ <= 0) {
            return 1.0;
        }

        return uncompressed_bytes_sent_/compressed_bytes_sent_;
    }

    double get_bytes_transmitted()
    {
        if (compressed_bytes_sent_ <= 0) {
            return header_bytes_sent_ + uncompressed_bytes_sent_;
        } else {
            return header_bytes_sent_ + compressed_bytes_sent_;
        }
    }
    
    void set_timeout(unsigned int t)
    {
        timeout_ms_ = t;
    }

    void read_task()
    {
        if (!socket_) {
            throw GadgetronClientException("Unable to create socket.");
        }

        while (socket_->is_open()) {

            GadgetMessageIdentifier id;
            boost::asio::read(*socket_, boost::asio::buffer(&id, sizeof(GadgetMessageIdentifier)));

            if (id.id == GADGET_MESSAGE_CLOSE) {
                break;
            }

            GadgetronClientMessageReader* r = find_reader(id.id);

            if (!r) {
                std::cout << "Message received with ID: " << id.id << std::endl;
                throw GadgetronClientException("Unknown Message ID");
            } else {
                r->read(socket_);
            }
        }
    }

    void wait() {
        reader_thread_.join();
    }

    void connect(std::string hostname, std::string port)
    {
        tcp::resolver resolver(io_service);
        // numeric_service flag is required to send data if the Linux machine has no internet connection (in this case the loopback device is the only network device with an address).
        // https://stackoverflow.com/questions/5971242/how-does-boost-asios-hostname-resolution-work-on-linux-is-it-possible-to-use-n
        // https://www.boost.org/doc/libs/1_65_0/doc/html/boost_asio/reference/ip__basic_resolver_query.html
        tcp::resolver::query query(tcp::v4(), hostname.c_str(), port.c_str(), boost::asio::ip::resolver_query_base::numeric_service);
        tcp::resolver::iterator endpoint_iterator = resolver.resolve(query);
        tcp::resolver::iterator end;

        socket_ = new tcp::socket(io_service);
        if (!socket_) {
            throw GadgetronClientException("Unable to create socket.");
        }

        std::condition_variable cv;
        std::mutex cv_m;
        
        boost::system::error_code error = boost::asio::error::host_not_found;
        std::thread t([&](){
                //TODO:
                //For newer versions of Boost, we should use
                //   boost::asio::connect(*socket_, iterator);
                while (error && endpoint_iterator != end) {
                    socket_->close();
                    socket_->connect(*endpoint_iterator++, error);
                }
                cv.notify_all();
            });

        {
            std::unique_lock<std::mutex> lk(cv_m);
            if (std::cv_status::timeout == cv.wait_until(lk, std::chrono::system_clock::now() +std::chrono::milliseconds(timeout_ms_)) ) {
                socket_->close();
             }
        }

        t.join();
        if (error)
            throw GadgetronClientException("Error connecting using socket.");

        reader_thread_ = std::thread([&](){this->read_task();});
    }

    void send_gadgetron_close() { 
        if (!socket_) {
            throw GadgetronClientException("Invalid socket.");
        }
        GadgetMessageIdentifier id;
        id.id = GADGET_MESSAGE_CLOSE;    
        header_bytes_sent_ += boost::asio::write(*socket_, boost::asio::buffer(&id, sizeof(GadgetMessageIdentifier)));
    }

    void send_gadgetron_info_query(const std::string &query, uint64_t correlation_id = 0) {
        GadgetMessageIdentifier id{ 6 }; // 6 = QUERY; Deal with it.

        header_bytes_sent_ += boost::asio::write(*socket_, boost::asio::buffer(&id, sizeof(id)));

        uint64_t reserved = 0;

        header_bytes_sent_ += boost::asio::write(*socket_, boost::asio::buffer(&reserved, sizeof(reserved)));
        header_bytes_sent_ += boost::asio::write(*socket_, boost::asio::buffer(&correlation_id, sizeof(correlation_id)));

        uint64_t query_length = query.size();

        header_bytes_sent_ += boost::asio::write(*socket_, boost::asio::buffer(&query_length, sizeof(query_length)));
        header_bytes_sent_ += boost::asio::write(*socket_, boost::asio::buffer(query));
    }

    void send_gadgetron_configuration_file(std::string config_xml_name) {

        if (!socket_) {
            throw GadgetronClientException("Invalid socket.");
        }

        GadgetMessageIdentifier id;
        id.id = GADGET_MESSAGE_CONFIG_FILE;

        GadgetMessageConfigurationFile ini;
        memset(&ini,0,sizeof(GadgetMessageConfigurationFile));
        strncpy(ini.configuration_file, config_xml_name.c_str(),config_xml_name.size());

        header_bytes_sent_ += boost::asio::write(*socket_, boost::asio::buffer(&id, sizeof(GadgetMessageIdentifier)));
        header_bytes_sent_ += boost::asio::write(*socket_, boost::asio::buffer(&ini, sizeof(GadgetMessageConfigurationFile)));

    }

    void send_gadgetron_configuration_script(std::string xml_string)
    {
        if (!socket_) {
            throw GadgetronClientException("Invalid socket.");
        }

        GadgetMessageIdentifier id;
        id.id = GADGET_MESSAGE_CONFIG_SCRIPT;

        GadgetMessageScript conf;
        conf.script_length = (uint32_t)xml_string.size();

        header_bytes_sent_ += boost::asio::write(*socket_, boost::asio::buffer(&id, sizeof(GadgetMessageIdentifier)));
        header_bytes_sent_ += boost::asio::write(*socket_, boost::asio::buffer(&conf, sizeof(GadgetMessageScript)));
        header_bytes_sent_ += boost::asio::write(*socket_, boost::asio::buffer(xml_string.c_str(), conf.script_length));
    }

    void  send_gadgetron_parameters(std::string xml_string)
    {
        if (!socket_) {
            throw GadgetronClientException("Invalid socket.");
        }

        GadgetMessageIdentifier id;
        id.id = GADGET_MESSAGE_PARAMETER_SCRIPT;

        GadgetMessageScript conf;
        conf.script_length = (uint32_t)xml_string.size();

        header_bytes_sent_ += boost::asio::write(*socket_, boost::asio::buffer(&id, sizeof(GadgetMessageIdentifier)));
        header_bytes_sent_ += boost::asio::write(*socket_, boost::asio::buffer(&conf, sizeof(GadgetMessageScript)));
        header_bytes_sent_ += boost::asio::write(*socket_, boost::asio::buffer(xml_string.c_str(), conf.script_length));
    }

    void send_ismrmrd_acquisition(ISMRMRD::Acquisition& acq) 
    {
        if (!socket_) {
            throw GadgetronClientException("Invalid socket.");
        }

        GadgetMessageIdentifier id;
        id.id = GADGET_MESSAGE_ISMRMRD_ACQUISITION;;

        header_bytes_sent_ += boost::asio::write(*socket_, boost::asio::buffer(&id, sizeof(GadgetMessageIdentifier)));
        header_bytes_sent_ += boost::asio::write(*socket_, boost::asio::buffer(&acq.getHead(), sizeof(ISMRMRD::AcquisitionHeader)));

        unsigned long trajectory_elements = acq.getHead().trajectory_dimensions*acq.getHead().number_of_samples;
        unsigned long data_elements = acq.getHead().active_channels*acq.getHead().number_of_samples;

        if (trajectory_elements) {
            header_bytes_sent_ += boost::asio::write(*socket_, boost::asio::buffer(&acq.getTrajPtr()[0], sizeof(float)*trajectory_elements));
        }


        if (data_elements) {
            uncompressed_bytes_sent_ +=boost::asio::write(*socket_, boost::asio::buffer(&acq.getDataPtr()[0], 2*sizeof(float)*data_elements));
        }
    }


    void send_ismrmrd_compressed_acquisition_precision(ISMRMRD::Acquisition& acq, unsigned int compression_precision) 
    {
        if (!socket_) {
            throw GadgetronClientException("Invalid socket.");
        }

        GadgetMessageIdentifier id;
        id.id = GADGET_MESSAGE_ISMRMRD_ACQUISITION;

        ISMRMRD::AcquisitionHeader h = acq.getHead(); //We will make a copy because we will be setting some flags
        h.setFlag(ISMRMRD::ISMRMRD_ACQ_COMPRESSION2);

        header_bytes_sent_ += boost::asio::write(*socket_, boost::asio::buffer(&id, sizeof(GadgetMessageIdentifier)));
        header_bytes_sent_ += boost::asio::write(*socket_, boost::asio::buffer(&h, sizeof(ISMRMRD::AcquisitionHeader)));

        unsigned long trajectory_elements = acq.getHead().trajectory_dimensions*acq.getHead().number_of_samples;
        unsigned long data_elements = acq.getHead().active_channels*acq.getHead().number_of_samples;

        if (trajectory_elements) {
            header_bytes_sent_ += boost::asio::write(*socket_, boost::asio::buffer(&acq.getTrajPtr()[0], sizeof(float)*trajectory_elements));
        }


        if (data_elements) {
            std::vector<float> input_data((float*)&acq.getDataPtr()[0], (float*)&acq.getDataPtr()[0] + acq.getHead().active_channels*acq.getHead().number_of_samples*2);

            CompressedBuffer<float> comp_buffer(input_data, -1.0, compression_precision);
            std::vector<uint8_t> serialized_buffer = comp_buffer.serialize();
 
            compressed_bytes_sent_ += serialized_buffer.size();
            uncompressed_bytes_sent_ += data_elements*2*sizeof(float);
                            
            uint32_t bs = (uint32_t)serialized_buffer.size();
            boost::asio::write(*socket_, boost::asio::buffer(&bs, sizeof(uint32_t)));
            boost::asio::write(*socket_, boost::asio::buffer(&serialized_buffer[0], serialized_buffer.size()));
        }
        
    }


    void send_ismrmrd_compressed_acquisition_tolerance(ISMRMRD::Acquisition& acq, float compression_tolerance, NoiseStatistics& stat) 
    {
        if (!socket_) {
            throw GadgetronClientException("Invalid socket.");
        }

        GadgetMessageIdentifier id;
        id.id = GADGET_MESSAGE_ISMRMRD_ACQUISITION;

        ISMRMRD::AcquisitionHeader h = acq.getHead(); //We will make a copy because we will be setting some flags
        h.setFlag(ISMRMRD::ISMRMRD_ACQ_COMPRESSION2);

        header_bytes_sent_ += boost::asio::write(*socket_, boost::asio::buffer(&id, sizeof(GadgetMessageIdentifier)));
        header_bytes_sent_ += boost::asio::write(*socket_, boost::asio::buffer(&h, sizeof(ISMRMRD::AcquisitionHeader)));

        unsigned long trajectory_elements = acq.getHead().trajectory_dimensions*acq.getHead().number_of_samples;
        unsigned long data_elements = acq.getHead().active_channels*acq.getHead().number_of_samples;

        if (trajectory_elements) {
            header_bytes_sent_ += boost::asio::write(*socket_, boost::asio::buffer(&acq.getTrajPtr()[0], sizeof(float)*trajectory_elements));
        }


        if (data_elements) {
            std::vector<float> input_data((float*)&acq.getDataPtr()[0], (float*)&acq.getDataPtr()[0] + acq.getHead().active_channels* acq.getHead().number_of_samples*2);

            float local_tolerance = compression_tolerance;
            float sigma = stat.sigma_min; //We use the minimum sigma of all channels to "cap" the error
            if (stat.status && sigma > 0 && stat.noise_dwell_time_us && acq.getHead().sample_time_us) {
                local_tolerance = local_tolerance*stat.sigma_min*acq.getHead().sample_time_us*std::sqrt(stat.noise_dwell_time_us/acq.getHead().sample_time_us);
            }

            CompressedBuffer<float> comp_buffer(input_data, local_tolerance);
            std::vector<uint8_t> serialized_buffer = comp_buffer.serialize();

            compressed_bytes_sent_ += serialized_buffer.size();
            uncompressed_bytes_sent_ += data_elements*2*sizeof(float);

            uint32_t bs = (uint32_t)serialized_buffer.size();
            boost::asio::write(*socket_, boost::asio::buffer(&bs, sizeof(uint32_t)));
            boost::asio::write(*socket_, boost::asio::buffer(&serialized_buffer[0], serialized_buffer.size()));
        }
    }

    void send_ismrmrd_zfp_compressed_acquisition_precision(ISMRMRD::Acquisition& acq, unsigned int compression_precision) 
    {

#if defined GADGETRON_COMPRESSION_ZFP
        if (!socket_) {
            throw GadgetronClientException("Invalid socket.");
        }

        GadgetMessageIdentifier id;
        //TODO: switch data type
        id.id = GADGET_MESSAGE_ISMRMRD_ACQUISITION;

        ISMRMRD::AcquisitionHeader h = acq.getHead(); //We will make a copy because we will be setting some flags
        h.setFlag(ISMRMRD::ISMRMRD_ACQ_COMPRESSION1);

        header_bytes_sent_ += boost::asio::write(*socket_, boost::asio::buffer(&id, sizeof(GadgetMessageIdentifier)));
        header_bytes_sent_ += boost::asio::write(*socket_, boost::asio::buffer(&h, sizeof(ISMRMRD::AcquisitionHeader)));

        unsigned long trajectory_elements = acq.getHead().trajectory_dimensions*acq.getHead().number_of_samples;
        unsigned long data_elements = acq.getHead().active_channels*acq.getHead().number_of_samples;

        if (trajectory_elements) {
            header_bytes_sent_ += boost::asio::write(*socket_, boost::asio::buffer(&acq.getTrajPtr()[0], sizeof(float)*trajectory_elements));
        }


        if (data_elements) {
            size_t comp_buffer_size = 4*sizeof(float)*data_elements;
            char* comp_buffer = new char[comp_buffer_size];
            size_t compressed_size = 0;
            try {
                compressed_size = compress_zfp_precision((float*)&acq.getDataPtr()[0],
                                                         acq.getHead().number_of_samples*2, acq.getHead().active_channels,
                                                         compression_precision, comp_buffer, comp_buffer_size);

                compressed_bytes_sent_ += compressed_size;
                uncompressed_bytes_sent_ += data_elements*2*sizeof(float);
                float compression_ratio = (1.0*data_elements*2*sizeof(float))/(float)compressed_size;
                //std::cout << "Compression ratio: " << compression_ratio << std::endl;
                
            } catch (...) {
                delete [] comp_buffer;
                std::cout << "Compression failure caught" << std::endl;
                throw;
            }


            //TODO: Write compressed buffer
            uint32_t bs = (uint32_t)compressed_size;
            boost::asio::write(*socket_, boost::asio::buffer(&bs, sizeof(uint32_t)));
            boost::asio::write(*socket_, boost::asio::buffer(comp_buffer, compressed_size));

            delete [] comp_buffer;
        }

#else //GADGETRON_COMPRESSION_ZFP
        throw GadgetronClientException("Attempting to do ZFP compression, but ZFP not available");
#endif //GADGETRON_COMPRESSION_ZFP
    }

    void send_ismrmrd_zfp_compressed_acquisition_tolerance(ISMRMRD::Acquisition& acq, float compression_tolerance, NoiseStatistics& stat) 
    {
#if defined GADGETRON_COMPRESSION_ZFP
        if (!socket_) {
            throw GadgetronClientException("Invalid socket.");
        }

        GadgetMessageIdentifier id;
        //TODO: switch data type
        id.id = GADGET_MESSAGE_ISMRMRD_ACQUISITION;

        ISMRMRD::AcquisitionHeader h = acq.getHead(); //We will make a copy because we will be setting some flags
        h.setFlag(ISMRMRD::ISMRMRD_ACQ_COMPRESSION1);

        header_bytes_sent_ += boost::asio::write(*socket_, boost::asio::buffer(&id, sizeof(GadgetMessageIdentifier)));
        header_bytes_sent_ += boost::asio::write(*socket_, boost::asio::buffer(&h, sizeof(ISMRMRD::AcquisitionHeader)));

        unsigned long trajectory_elements = acq.getHead().trajectory_dimensions*acq.getHead().number_of_samples;
        unsigned long data_elements = acq.getHead().active_channels*acq.getHead().number_of_samples;

        if (trajectory_elements) {
            header_bytes_sent_ += boost::asio::write(*socket_, boost::asio::buffer(&acq.getTrajPtr()[0], sizeof(float)*trajectory_elements));
        }

        float local_tolerance = compression_tolerance;
        float sigma = stat.sigma_min; //We use the minimum sigma of all channels to "cap" the error
        if (stat.status && sigma > 0 && stat.noise_dwell_time_us && acq.getHead().sample_time_us) {
            local_tolerance = local_tolerance*stat.sigma_min*acq.getHead().sample_time_us*std::sqrt(stat.noise_dwell_time_us/acq.getHead().sample_time_us);
        }

        if (data_elements) {
            size_t comp_buffer_size = 4*sizeof(float)*data_elements;
            char* comp_buffer = new char[comp_buffer_size];
            size_t compressed_size = 0;
            try {
                compressed_size = compress_zfp_tolerance((float*)&acq.getDataPtr()[0],
                                                         acq.getHead().number_of_samples*2, acq.getHead().active_channels,
                                                         local_tolerance, comp_buffer, comp_buffer_size);

                compressed_bytes_sent_ += compressed_size;
                uncompressed_bytes_sent_ += data_elements*2*sizeof(float);
                float compression_ratio = (1.0*data_elements*2*sizeof(float))/(float)compressed_size;
                //std::cout << "Compression ratio: " << compression_ratio << std::endl;
                
            } catch (...) {
                delete [] comp_buffer;
                std::cout << "Compression failure caught" << std::endl;
                throw;
            }


            //TODO: Write compressed buffer
            uint32_t bs = (uint32_t)compressed_size;
            boost::asio::write(*socket_, boost::asio::buffer(&bs, sizeof(uint32_t)));
            boost::asio::write(*socket_, boost::asio::buffer(comp_buffer, compressed_size));

            delete [] comp_buffer;
        }
#else //GADGETRON_COMPRESSION_ZFP
        throw GadgetronClientException("Attempting to do ZFP compression, but ZFP not available");
#endif //GADGETRON_COMPRESSION_ZFP

    }

    void send_ismrmrd_waveform(ISMRMRD::Waveform& wav)
    {
        if (!socket_)
        {
            throw GadgetronClientException("Invalid socket.");
        }

        GadgetMessageIdentifier id;
        id.id = GADGET_MESSAGE_ISMRMRD_WAVEFORM;;

        header_bytes_sent_ += boost::asio::write(*socket_, boost::asio::buffer(&id, sizeof(GadgetMessageIdentifier)));
        header_bytes_sent_ += boost::asio::write(*socket_, boost::asio::buffer(&wav.head, sizeof(ISMRMRD::ISMRMRD_WaveformHeader)));

        unsigned long data_elements = wav.head.channels*wav.head.number_of_samples;

        if (data_elements)
        {
            header_bytes_sent_ += boost::asio::write(*socket_, boost::asio::buffer(wav.begin_data(), sizeof(uint32_t)*data_elements));
        }
    }

    void register_reader(unsigned short slot, std::shared_ptr<GadgetronClientMessageReader> r) {
        readers_[slot] = r;
    }

protected:
    typedef std::map<unsigned short, std::shared_ptr<GadgetronClientMessageReader> > maptype;

    GadgetronClientMessageReader* find_reader(unsigned short r)
    {
        GadgetronClientMessageReader* ret = 0;

        maptype::iterator it = readers_.find(r);

        if (it != readers_.end()) {
            ret = it->second.get();
        }

        return ret;
    }

    boost::asio::io_service io_service;
    tcp::socket* socket_;
    std::thread reader_thread_;
    maptype readers_;
    unsigned int timeout_ms_;
    double header_bytes_sent_;
    double uncompressed_bytes_sent_;
    double compressed_bytes_sent_;
};


class GadgetronClientQueryToStringReader : public GadgetronClientMessageReader
{
  
public:
  GadgetronClientQueryToStringReader(std::string& result) : result_(result)
  {
    
  }

  virtual ~GadgetronClientQueryToStringReader()
  {
    
  }

  virtual void read(tcp::socket* stream)
  {
    size_t recv_count = 0;
    
    typedef unsigned long long size_t_type;
    size_t_type len(0);
    boost::asio::read(*stream, boost::asio::buffer(&len, sizeof(size_t_type)));

    std::vector<char> temp(len,0);
    if (boost::asio::read(*stream, boost::asio::buffer(temp.data(), len)) != len)
    {
      throw GadgetronClientException("Incorrect number of bytes read for dependency query");
    }
    result_ = std::string(temp.data(),len);
    
  }
  
  protected:
    std::string& result_;
};


NoiseStatistics get_noise_statistics(std::string dependency_name, std::string host_name, std::string port, unsigned int timeout_ms)
{
    GadgetronClientConnector con;
    con.set_timeout(timeout_ms);
    std::string result;
    NoiseStatistics stat;

    con.register_reader(GADGET_MESSAGE_DEPENDENCY_QUERY, std::make_shared<GadgetronClientQueryToStringReader>(result));
    
    std::string xml_config;
    
    xml_config += "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    xml_config += "      <gadgetronStreamConfiguration xsi:schemaLocation=\"http://gadgetron.sf.net/gadgetron gadgetron.xsd\"\n";
    xml_config += "        xmlns=\"http://gadgetron.sf.net/gadgetron\"\n";
    xml_config += "      xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\">\n";
    xml_config += "\n";
    xml_config += "    <writer>\n";
    xml_config += "      <slot>1019</slot>\n";
    xml_config += "      <dll>gadgetron_mricore</dll>\n";
    xml_config += "      <classname>DependencyQueryWriter</classname>\n";
    xml_config += "    </writer>\n";
    xml_config += "\n";
    xml_config += "    <gadget>\n";
    xml_config += "      <name>NoiseSummary</name>\n";
    xml_config += "      <dll>gadgetron_mricore</dll>\n";
    xml_config += "      <classname>NoiseSummaryGadget</classname>\n";
    xml_config += "\n";
    xml_config += "      <property>\n";
    xml_config += "         <name>noise_file</name>\n";
    xml_config += "         <value>" + dependency_name + "</value>\n";
    xml_config += "      </property>\n";
    xml_config += "    </gadget>\n";
    xml_config += "\n";
    xml_config += "</gadgetronStreamConfiguration>\n";

    try {
        con.connect(host_name,port);
        con.send_gadgetron_configuration_script(xml_config);
        con.send_gadgetron_close();
        con.wait();
    } catch (...) {
        std::cerr << "Unable to retrive noise statistics from server" << std::endl;
        stat.status = false;
    }

    try {
        ISMRMRD::MetaContainer meta;
        ISMRMRD::deserialize(result.c_str(), meta);
        stat.status = meta.as_str("status") == std::string("success");
        stat.channels = meta.as_long("channels");
        stat.sigma_min = meta.as_double("min_sigma");
        stat.sigma_max = meta.as_double("max_sigma");
        stat.sigma_mean = meta.as_double("mean_sigma");
        stat.noise_dwell_time_us = meta.as_double("noise_dwell_time_us");
    } catch (...) {
        stat.status = false;
    }

    return stat;
}

void send_ismrmrd_acq(GadgetronClientConnector& con, ISMRMRD::Acquisition& acq_tmp, 
    unsigned int compression_precision, bool use_zfp_compression, float compression_tolerance, NoiseStatistics& noise_stats)
{
    try
    {
        if (compression_precision > 0)
        {
            if (use_zfp_compression) {
                con.send_ismrmrd_zfp_compressed_acquisition_precision(acq_tmp, compression_precision);
            }
            else {
                con.send_ismrmrd_compressed_acquisition_precision(acq_tmp, compression_precision);
            }
        }
        else if (compression_tolerance > 0.0)
        {
            if (use_zfp_compression) {
                con.send_ismrmrd_zfp_compressed_acquisition_tolerance(acq_tmp, compression_tolerance, noise_stats);
            }
            else {
                con.send_ismrmrd_compressed_acquisition_tolerance(acq_tmp, compression_tolerance, noise_stats);
            }
        }
        else
        {
            con.send_ismrmrd_acquisition(acq_tmp);
        }
    }
    catch(...)
    {
        throw GadgetronClientException("send_ismrmrd_acq failed ... ");
    }
}

int main(int argc, char **argv)
{

    std::string host_name;
    std::string port;
    std::string in_filename;
    std::string out_filename;
    std::string hdf5_in_group;
    std::string hdf5_out_group;
    std::string config_file;
    std::string config_file_local;
    std::string config_xml_local;
    unsigned int loops;
    unsigned int timeout_ms;
    std::string out_fileformat;
    bool open_input_file = true;
    unsigned int compression_precision = 0;
    float compression_tolerance = 0.0;
    bool use_zfp_compression = false;
    bool verbose = false;
    Gadgetron::GadgetronTimer timer(false);

    po::options_description desc("Allowed options");

    desc.add_options()
        ("help,h", "Produce help message")
        ("query,q", "Dependency query mode")
        ("verbose,v", "Verbose mode")
        ("info,Q", po::value<std::string>(), "Query Gadgetron information")
        ("port,p", po::value<std::string>(&port)->default_value("9002"), "Port")
        ("address,a", po::value<std::string>(&host_name)->default_value("localhost"), "Address (hostname) of Gadgetron host")
        ("filename,f", po::value<std::string>(&in_filename), "Input file")
        ("outfile,o", po::value<std::string>(&out_filename)->default_value("out.h5"), "Output file")
        ("in-group,g", po::value<std::string>(&hdf5_in_group)->default_value("/dataset"), "Input data group")
        ("out-group,G", po::value<std::string>(&hdf5_out_group)->default_value(get_date_time_string()), "Output group name")  
        ("config,c", po::value<std::string>(&config_file)->default_value("default.xml"), "Configuration file (remote)")
        ("config-local,C", po::value<std::string>(&config_file_local), "Configuration file (local)")
        ("loops,l", po::value<unsigned int>(&loops)->default_value(1), "Loops")
        ("timeout,t", po::value<unsigned int>(&timeout_ms)->default_value(10000), "Timeout [ms]")
        ("outformat,F", po::value<std::string>(&out_fileformat)->default_value("h5"), "Out format, h5 for hdf5 and hdr for analyze image")
        ("precision,P", po::value<unsigned int>(&compression_precision)->default_value(0), "Compression precision (bits)")
        ("tolerance,T", po::value<float>(&compression_tolerance)->default_value(0.0), "Compression tolerance (fraction of sigma, if no noise stats, assume sigma 1)")
#if defined GADGETRON_COMPRESSION_ZFP
        ("ZFP,Z", po::value<bool>(&use_zfp_compression)->default_value(false), "Use ZFP library for compression");
#endif //GADGETRON_COMPRESSION_ZFP
        ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << desc << std::endl;
        return 0;
    }

    if (!vm.count("filename") && !vm.count("query") && !vm.count("info")) {
        std::cout << std::endl << std::endl << "\tYou must supply a filename" << std::endl << std::endl;
        std::cout << desc << std::endl;
        return -1;
    }

    if (vm.count("query") || vm.count("info")) {
        open_input_file = false;
    }

    if (vm.count("verbose")) {
        verbose = true;
    }

    if (vm.count("config-local")) {
        std::ifstream t(config_file_local.c_str());
        if (t) {
            //Read in the file.
            config_xml_local = std::string((std::istreambuf_iterator<char>(t)),
                std::istreambuf_iterator<char>());
        } else {
            std::cout << "Unable to read local xml configuration: " << config_file_local  << std::endl;
            return -1;
        }
    }

    if (compression_precision > 0 && compression_tolerance > 0.0) {
       std::cout << "You cannot supply both compression precision (P) and compression tolerance (T) at the same time" << std::endl;
       return -1;
    }

    //Let's check if the files exist:
    std::string hdf5_xml_varname = std::string(hdf5_in_group) + std::string("/xml");
    std::string hdf5_data_varname = std::string(hdf5_in_group) + std::string("/data");
    std::string hdf5_waveform_varname = std::string(hdf5_in_group) + std::string("/waveform");

    //TODO:
    // Add check to see if input file exists

    //Let's open the input file
    std::shared_ptr<ISMRMRD::Dataset> ismrmrd_dataset;
    std::string xml_config;
    if (open_input_file) {
      ismrmrd_dataset = std::shared_ptr<ISMRMRD::Dataset>(new ISMRMRD::Dataset(in_filename.c_str(), hdf5_in_group.c_str(), false));
      // Read the header
      ismrmrd_dataset->readHeader(xml_config);
    }

    if (!vm.count("query") && !vm.count("info")) {
      std::cout << "Gadgetron ISMRMRD client" << std::endl;
      std::cout << "  -- host            :      " << host_name << std::endl;
      std::cout << "  -- port            :      " << port << std::endl;
      std::cout << "  -- hdf5 file  in   :      " << in_filename << std::endl;
      std::cout << "  -- hdf5 group in   :      " << hdf5_in_group << std::endl;
      std::cout << "  -- conf            :      " << config_file << std::endl;
      std::cout << "  -- loop            :      " << loops << std::endl;
      std::cout << "  -- hdf5 file out   :      " << out_filename << std::endl;
      std::cout << "  -- hdf5 group out  :      " << hdf5_out_group << std::endl;
    }


    //Let's figure out if this measurement has dependencies
    NoiseStatistics noise_stats; noise_stats.status = false;
    if (!vm.count("query") && !vm.count("info")) {
        ISMRMRD::IsmrmrdHeader h;
        ISMRMRD::deserialize(xml_config.c_str(),h);

        std::string noise_id;
        if (h.measurementInformation.is_present() &&
            (h.measurementInformation().measurementDependency.size() > 0)) {
            std::cout << "This measurement has dependent measurements" << std::endl;
            for (auto d: h.measurementInformation().measurementDependency) {
                std::cout << "  " << d.dependencyType << " : " << d.measurementID << std::endl;
                if (d.dependencyType == "Noise") {
                    noise_id = d.measurementID;
                }
            }
        }
        
        if (!noise_id.empty()) {
            std::cout << "Querying the Gadgetron instance for the dependent measurement: " << noise_id << std::endl;
            noise_stats = get_noise_statistics(std::string("GadgetronNoiseCovarianceMatrix_") + noise_id, host_name, port, timeout_ms);
            if (!noise_stats.status) {
                std::cout << "WARNING: Dependent noise measurement not found on Gadgetron server. Was the noise data processed?" << std::endl;
                if (compression_tolerance > 0.0) {
                    std::cout << "  !!!!!! COMPRESSION TOLERANCE LEVEL SPECIFIED, BUT IT IS NOT POSSIBLE TO DETERMINE SIGMA. ASSIMUMING SIGMA == 1 !!!!!!" << std::endl;
                }
            } else {
                std::cout << "Noise level: Min sigma = " << noise_stats.sigma_min << ", Mean sigma = " << noise_stats.sigma_mean << ", Max sigma = " << noise_stats.sigma_max << std::endl; 
            }
        }
    }

    GadgetronClientConnector con;
    con.set_timeout(timeout_ms);

    if ( out_fileformat == "hdr" )
    {
        con.register_reader(GADGET_MESSAGE_ISMRMRD_IMAGE, std::shared_ptr<GadgetronClientMessageReader>(new GadgetronClientAnalyzeImageMessageReader(hdf5_out_group)));
    }
    else
    {
        con.register_reader(GADGET_MESSAGE_ISMRMRD_IMAGE, std::shared_ptr<GadgetronClientMessageReader>(new GadgetronClientImageMessageReader(out_filename, hdf5_out_group)));
    }

    con.register_reader(GADGET_MESSAGE_DICOM_WITHNAME, std::shared_ptr<GadgetronClientMessageReader>(new GadgetronClientBlobMessageReader(std::string(hdf5_out_group), std::string("dcm"))));

    con.register_reader(GADGET_MESSAGE_DEPENDENCY_QUERY, std::shared_ptr<GadgetronClientMessageReader>(new GadgetronClientDependencyQueryReader(std::string(out_filename))));
    con.register_reader(GADGET_MESSAGE_TEXT, std::shared_ptr<GadgetronClientMessageReader>(new GadgetronClientTextReader()));
    con.register_reader(7, std::shared_ptr<GadgetronClientResponseReader>(new GadgetronClientResponseReader()));

    try
    {
        timer.start();
        con.connect(host_name,port);

        if (vm.count("info")) {
            con.send_gadgetron_info_query(vm["info"].as<std::string>());
        }
        else if (vm.count("config-local"))
        {
            con.send_gadgetron_configuration_script(config_xml_local);
        }
        else
        {
            con.send_gadgetron_configuration_file(config_file);
        }

        if (open_input_file)
        {
            con.send_gadgetron_parameters(xml_config);

            uint32_t acquisitions = 0;
            {
                mtx.lock();
                acquisitions = ismrmrd_dataset->getNumberOfAcquisitions();
                mtx.unlock();
            }

            uint32_t waveforms = 0;
            {
                mtx.lock();
                waveforms = ismrmrd_dataset->getNumberOfWaveforms();
                mtx.unlock();
            }

            if(verbose)
            {
                std::cout << "Find " << acquisitions << " ismrmrd acquisitions" << std::endl;
                std::cout << "Find " << waveforms << " ismrmrd waveforms" << std::endl;
            }

            ISMRMRD::Acquisition acq_tmp;
            ISMRMRD::Waveform wav_tmp;

            uint32_t i(0), j(0); // i : index over the acquisition; j : index over the waveform

            if(waveforms>0)
            {
                {
                    std::lock_guard<std::mutex> scoped_lock(mtx);
                    ismrmrd_dataset->readAcquisition(i, acq_tmp);
                }

                {
                    std::lock_guard<std::mutex> scoped_lock(mtx);
                    ismrmrd_dataset->readWaveform(j, wav_tmp);
                }

                while(i<acquisitions && j<waveforms)
                {
                    while(wav_tmp.head.time_stamp < acq_tmp.getHead().acquisition_time_stamp)
                    {
                        con.send_ismrmrd_waveform(wav_tmp);

                        if (verbose)
                        {
                            std::cout << "--> Send out ismrmrd waveform : " << j << " - " << wav_tmp.head.scan_counter 
                                << " - " << wav_tmp.head.time_stamp 
                                << " - " << wav_tmp.head.channels
                                << " - " << wav_tmp.head.number_of_samples
                                << " - " << wav_tmp.head.waveform_id
                                << std::endl;
                        }

                        j++;

                        if(j<waveforms)
                        {
                            std::lock_guard<std::mutex> scoped_lock(mtx);
                            ismrmrd_dataset->readWaveform(j, wav_tmp);
                        }
                        else
                        {
                            break;
                        }
                    }

                    while (acq_tmp.getHead().acquisition_time_stamp <= wav_tmp.head.time_stamp)
                    {
                        send_ismrmrd_acq(con, acq_tmp, compression_precision, use_zfp_compression, compression_tolerance, noise_stats);

                        if (verbose)
                        {
                            std::cout << "==> Send out ismrmrd acq : " << i << " - " << acq_tmp.getHead().scan_counter << " - " << acq_tmp.getHead().acquisition_time_stamp << std::endl;
                        }

                        i++;

                        if(i<acquisitions)
                        {
                            std::lock_guard<std::mutex> scoped_lock(mtx);
                            ismrmrd_dataset->readAcquisition(i, acq_tmp);
                        }
                        else
                        {
                            break;
                        }
                    }

                    if(j==waveforms && i<acquisitions)
                    {
                        send_ismrmrd_acq(con, acq_tmp, compression_precision, use_zfp_compression, compression_tolerance, noise_stats);

                        for (uint32_t ia=i+1; ia<acquisitions; ia++)
                        {
                            {
                                std::lock_guard<std::mutex> scoped_lock(mtx);
                                ismrmrd_dataset->readAcquisition(ia, acq_tmp);
                            }

                            send_ismrmrd_acq(con, acq_tmp, compression_precision, use_zfp_compression, compression_tolerance, noise_stats);
                        }
                    }

                    if (i==acquisitions && j<waveforms)
                    {
                        con.send_ismrmrd_waveform(wav_tmp);

                        for (uint32_t iw = j + 1; iw<waveforms; iw++)
                        {
                            {
                                std::lock_guard<std::mutex> scoped_lock(mtx);
                                ismrmrd_dataset->readWaveform(iw, wav_tmp);
                            }

                            con.send_ismrmrd_waveform(wav_tmp);
                        }
                    }
                }
            }
            else
            {
                for (i=0; i<acquisitions; i++)
                {
                    {
                        std::lock_guard<std::mutex> scoped_lock(mtx);
                        ismrmrd_dataset->readAcquisition(i, acq_tmp);
                    }

                    send_ismrmrd_acq(con, acq_tmp, compression_precision, use_zfp_compression, compression_tolerance, noise_stats);
                }
            }
        }

        if (compression_precision > 0 || compression_tolerance > 0.0) {
            std::cout << "Compression ratio: " << con.compression_ratio() << std::endl;
        }

        if (verbose) {
            double transmission_time_s = timer.stop()/1e6;
            double transmitted_mb = con.get_bytes_transmitted()/(1024*1024);
            std::cout << "Time sending: " << transmission_time_s << "s" << std::endl;
            std::cout << "Data sent: " << transmitted_mb << "MB" << std::endl;
            std::cout << "Transmission rate: " << transmitted_mb/transmission_time_s << "MB/s" << std::endl;
        }

        con.send_gadgetron_close();
        con.wait();
    }
    catch (std::exception& ex)
    {
        std::cerr << "Error caught: " << ex.what() << std::endl;
        return -1;
    }

    return 0;
}
