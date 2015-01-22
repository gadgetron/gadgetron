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
// -Add command line switch for controlling verbosity of output
// -Static linking for standalone executable. 

#include <boost/program_options.hpp>
#include <boost/asio.hpp>
#include <boost/thread/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/shared_ptr.hpp>

#include <ismrmrd/ismrmrd.h>
#include <ismrmrd/dataset.h>
#include <ismrmrd/meta.h>

#include <fstream>
#include <streambuf>
#include <time.h>
#include <iomanip>
#include <sstream>
#include <iostream>
#include <exception>
#include <map>


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


namespace po = boost::program_options;
using boost::asio::ip::tcp;


enum GadgetronMessageID {
    GADGET_MESSAGE_INT_ID_MIN                             =   0,
    GADGET_MESSAGE_CONFIG_FILE                            =   1,
    GADGET_MESSAGE_CONFIG_SCRIPT                          =   2,
    GADGET_MESSAGE_PARAMETER_SCRIPT                       =   3,
    GADGET_MESSAGE_CLOSE                                  =   4,
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
    GADGET_MESSAGE_ISMRMRD_IMAGE_CPLX_FLOAT               = 1009,
    GADGET_MESSAGE_ISMRMRD_IMAGE_REAL_FLOAT               = 1010,
    GADGET_MESSAGE_ISMRMRD_IMAGE_REAL_USHORT              = 1011,
    GADGET_MESSAGE_DICOM                                  = 1012,
    GADGET_MESSAGE_CLOUD_JOB                              = 1013,
    GADGET_MESSAGE_GADGETCLOUD_JOB                        = 1014,
    GADGET_MESSAGE_ISMRMRD_IMAGEWITHATTRIB_CPLX_FLOAT     = 1015,
    GADGET_MESSAGE_ISMRMRD_IMAGEWITHATTRIB_REAL_FLOAT     = 1016,
    GADGET_MESSAGE_ISMRMRD_IMAGEWITHATTRIB_REAL_USHORT    = 1017,
    GADGET_MESSAGE_DICOM_WITHNAME                         = 1018,
    GADGET_MESSAGE_DEPENDENCY_QUERY                       = 1019,
    GADGET_MESSAGE_EXT_ID_MAX                             = 4096
};

boost::mutex mtx;

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


template <typename T> class GadgetronClientImageMessageReader 
    : public GadgetronClientMessageReader
{

public:
    GadgetronClientImageMessageReader(std::string filename, std::string groupname)
        : file_name_(filename)
        , group_name_(groupname)
    {

    }

    ~GadgetronClientImageMessageReader() {
    } 

    virtual void read(tcp::socket* stream) 
    {
        //std::cout << "Receiving image." << std::endl;
        //Read the image from the socket
        ISMRMRD::ImageHeader h;
        boost::asio::read(*stream, boost::asio::buffer(&h,sizeof(ISMRMRD::ImageHeader)));

        // TODO check the datatype!
        ISMRMRD::Image<T> im; 
        im.setHead(h);
        boost::asio::read(*stream, boost::asio::buffer(im.getDataPtr(), im.getDataSize()));
        {
            if (!dataset_) {

                {
                    mtx.lock();
                    dataset_ = boost::shared_ptr<ISMRMRD::Dataset>(new ISMRMRD::Dataset(file_name_.c_str(), group_name_.c_str(), true)); // create if necessary 
                    mtx.unlock();
                }
            }

            std::stringstream st1;
            st1 << "image_" << h.image_series_index;
            std::string image_varname = st1.str();

            {
                mtx.lock();
                // TODO should this be wrapped in a try/catch?
                dataset_->appendImage(image_varname, im);
                mtx.unlock();
            }

        }
    }

protected:
    std::string group_name_;
    std::string file_name_;
    boost::shared_ptr<ISMRMRD::Dataset> dataset_;
};

template <typename T> class GadgetronClientAttribImageMessageReader 
    : public GadgetronClientMessageReader
{

public:
    GadgetronClientAttribImageMessageReader(std::string filename, std::string groupname)
        : file_name_(filename)
        , group_name_(groupname)
    {

    }

    ~GadgetronClientAttribImageMessageReader() {
    } 

    virtual void read(tcp::socket* stream) 
    {
        //std::cout << "Receiving image with attributes." << std::endl;
        //Read the image headerfrom the socket
        ISMRMRD::ImageHeader h;
        boost::asio::read(*stream, boost::asio::buffer(&h,sizeof(ISMRMRD::ImageHeader)));
        ISMRMRD::Image<T> im;
        im.setHead(h);

        typedef unsigned long long size_t_type;

        //Read meta attributes
        size_t_type meta_attrib_length;
        boost::asio::read(*stream, boost::asio::buffer(&meta_attrib_length, sizeof(size_t_type)));

        std::string meta_attrib(meta_attrib_length,0);
        boost::asio::read(*stream, boost::asio::buffer(const_cast<char*>(meta_attrib.c_str()), meta_attrib_length));
        im.setAttributeString(meta_attrib);

        //Read image data
        boost::asio::read(*stream, boost::asio::buffer(im.getDataPtr(), im.getDataSize()));
        {
            if (!dataset_) {

                {
                    mtx.lock();
                    dataset_ = boost::shared_ptr<ISMRMRD::Dataset>(new ISMRMRD::Dataset(file_name_.c_str(), group_name_.c_str(), true)); // create if necessary 
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

protected:
    std::string group_name_;
    std::string file_name_;
    boost::shared_ptr<ISMRMRD::Dataset> dataset_;
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
            if ( pixelSize.size() > 1 )
                header.dime.pixdim[1] = pixelSize[0];
            if ( pixelSize.size() > 2 )
                header.dime.pixdim[2] = pixelSize[1];
            if ( pixelSize.size() > 3 )
                header.dime.pixdim[3] = pixelSize[2];
            if ( pixelSize.size() > 4 )
                header.dime.pixdim[4] = pixelSize[3];
            if ( pixelSize.size() > 5 )
                header.dime.pixdim[5] = pixelSize[4];
            if ( pixelSize.size() > 6 )
                header.dime.pixdim[6] = pixelSize[5];
            if ( pixelSize.size() > 7 )
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

template <typename T> class GadgetronClientAnalyzeImageMessageReader 
    : public GadgetronClientMessageReader
{

public:
    GadgetronClientAnalyzeImageMessageReader(const std::string& prefix=std::string("Image")) : prefix_(prefix)
    {

    }

    ~GadgetronClientAnalyzeImageMessageReader() {
    } 

    virtual void read(tcp::socket* stream) 
    {
        using namespace ISMRMRD;

        //Read the image from the socket
        ISMRMRD::ImageHeader h;
        boost::asio::read(*stream, boost::asio::buffer(&h,sizeof(ISMRMRD::ImageHeader)));
        ISMRMRD::Image<T> im; 
        im.setHead(h);
        boost::asio::read(*stream, boost::asio::buffer(im.getDataPtr(), im.getDataSize()));

        std::cout << "Receiving image : " << h.image_series_index << " - " << h.image_index << std::endl;
        {
            // analyze header
            std::stringstream st1;
            st1 << prefix_ << "_" << h.image_series_index << "_" << h.image_index << ".hdr";
            std::string head_varname = st1.str();

            std::vector<size_t> dim(3);
            dim[0] = h.matrix_size[0];
            dim[1] = h.matrix_size[1];
            dim[2] = h.matrix_size[2];

            std::vector<float> pixelSize(3);
            pixelSize[0] = h.field_of_view[0]/h.matrix_size[0];
            pixelSize[1] = h.field_of_view[1]/h.matrix_size[1];
            pixelSize[2] = h.field_of_view[2]/h.matrix_size[2];

            IOAnalyze hdr;
            dsr header;
            hdr.array2Header<T>(dim, pixelSize, header);

            std::ofstream outfileHeader;
            outfileHeader.open (head_varname.c_str(), std::ios::out|std::ios::binary);
            outfileHeader.write(reinterpret_cast<const char*>(&header), sizeof(dsr));
            outfileHeader.close();

            // data
            std::stringstream st2;
            st2 << prefix_ << "_" << h.image_series_index << "_" << h.image_index << ".img";
            std::string img_varname = st2.str();

            std::ofstream outfileData;
            outfileData.open (img_varname.c_str(), std::ios::out|std::ios::binary);
            outfileData.write(reinterpret_cast<const char*>(im.getDataPtr()), sizeof(T)*dim[0]*dim[1]*dim[2]);
            outfileData.close();
        }
    }

protected:

    std::string prefix_;
};

template <typename T> class GadgetronClientAttribAnalyzeImageMessageReader 
    : public GadgetronClientMessageReader
{

public:
    GadgetronClientAttribAnalyzeImageMessageReader(const std::string& prefix=std::string("Image")) : prefix_(prefix)
    {

    }

    ~GadgetronClientAttribAnalyzeImageMessageReader() {
    } 

    virtual void read(tcp::socket* stream) 
    {
        //Read the image headerfrom the socket
        ISMRMRD::ImageHeader h;
        boost::asio::read(*stream, boost::asio::buffer(&h,sizeof(ISMRMRD::ImageHeader)));
        ISMRMRD::Image<T> im; 
        im.setHead(h);

        std::cout << "Receiving image with attributes : " << h.image_series_index << " - " << h.image_index << std::endl;

        typedef unsigned long long size_t_type;

        //Read meta attributes
        size_t_type meta_attrib_length;
        boost::asio::read(*stream, boost::asio::buffer(&meta_attrib_length, sizeof(size_t_type)));

        std::string meta_attrib(meta_attrib_length,0);
        boost::asio::read(*stream, boost::asio::buffer(const_cast<char*>(meta_attrib.c_str()), meta_attrib_length));

        //Read image data
        boost::asio::read(*stream, boost::asio::buffer(im.getDataPtr(), im.getDataSize()));
        {
            // deserialize the meta attribute
            ISMRMRD::MetaContainer imgAttrib;
            ISMRMRD::deserialize(meta_attrib.c_str(), imgAttrib);

            size_t n;
            size_t num = imgAttrib.length("GADGETRON_DataRole");

            std::vector<std::string> dataRole;
            if ( num == 0 )
            {
                dataRole.push_back("GADGETRON_Image");
            }
            else
            {
                dataRole.resize(num);
                for ( n=0; n<num; n++ )
                {
                    dataRole[n] = std::string( imgAttrib.as_str("GADGETRON_DataRole", n) );
                }
            }

            long imageNumber = imgAttrib.as_long("GADGETRON_ImageNumber", 0);

            long cha, slc, e2, con, phs, rep, set, ave;
            cha = imgAttrib.as_long("CHA",          0);
            slc = imgAttrib.as_long("SLC",          0);
            e2  = imgAttrib.as_long("E2",           0);
            con = imgAttrib.as_long("CON",          0);
            phs = imgAttrib.as_long("PHS",          0);
            rep = imgAttrib.as_long("REP",          0);
            set = imgAttrib.as_long("SET",          0);
            ave = imgAttrib.as_long("AVE",          0);

            std::ostringstream ostr;

            if ( !prefix_.empty() )
            {
                ostr << prefix_ << "_";
            }

            for ( n=0; n<dataRole.size(); n++ )
            {
                ostr << dataRole[n] << "_";
            }

            ostr << "SLC" << slc << "_"
                 << "E2"  << e2  << "_"
                 << "CON" << con << "_"
                 << "PHS" << phs << "_"
                 << "REP" << rep << "_"
                 << "SET" << set << "_"
                 << "AVE" << ave << "_"
                 << "CHA" << cha << "_" 
                 << "ImageSeries" << h.image_series_index;

            std::string filename = ostr.str();

            // analyze header
            std::stringstream st1;
            st1 << filename << ".hdr";
            std::string head_varname = st1.str();

            std::vector<size_t> dim(3);
            dim[0] = h.matrix_size[0];
            dim[1] = h.matrix_size[1];
            dim[2] = h.matrix_size[2];

            std::vector<float> pixelSize(3);
            pixelSize[0] = h.field_of_view[0]/h.matrix_size[0];
            pixelSize[1] = h.field_of_view[1]/h.matrix_size[1];
            pixelSize[2] = h.field_of_view[2]/h.matrix_size[2];

            IOAnalyze hdr;
            dsr header;
            hdr.array2Header<T>(dim, pixelSize, header);

            std::ofstream outfileHeader;
            outfileHeader.open (head_varname.c_str(), std::ios::out|std::ios::binary);
            outfileHeader.write(reinterpret_cast<const char*>(&header), sizeof(dsr));
            outfileHeader.close();

            // data
            std::stringstream st2;
            st2 << filename << ".img";
            std::string img_varname = st2.str();

            std::ofstream outfileData;
            outfileData.open (img_varname.c_str(), std::ios::out|std::ios::binary);
            outfileData.write(reinterpret_cast<const char*>(im.getDataPtr()), sizeof(T)*dim[0]*dim[1]*dim[2]);
            outfileData.close();

            // attribute
            std::stringstream st3;
            st3 << filename << ".attrib";
            std::string meta_varname = st3.str();

            std::ofstream outfile;
            outfile.open (meta_varname.c_str(), std::ios::out|std::ios::binary);
            outfile.write(meta_attrib.c_str(), meta_attrib_length);
            outfile.close();
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

        std::stringstream filename;

        // Create the filename: (prefix_%06.suffix)
        filename << file_prefix << "_";
        filename << std::setfill('0') << std::setw(MAX_BLOBS_LOG_10) << number_of_calls_;
        filename << "." << file_suffix;

        std::ofstream outfile;
        outfile.open (filename.str().c_str(), std::ios::out|std::ios::binary);

        std::cout << "Writing image " << filename.str() << std::endl;

        if (outfile.good()) {
            /* write 'size' bytes starting at 'data's pointer */
            outfile.write(&data[0], nbytes);
            outfile.close();
            number_of_calls_++;
        } else {
            throw GadgetronClientException("Unable to write blob to output file\n");
        }
    }

protected:
    size_t number_of_calls_;
    std::string file_prefix;
    std::string file_suffix;

};

class GadgetronClientBlobAttribMessageReader 
    : public GadgetronClientMessageReader
{

public:
    GadgetronClientBlobAttribMessageReader(std::string fileprefix, std::string filesuffix)
        : number_of_calls_(0)
        , file_prefix(fileprefix)
        , file_suffix(filesuffix)

    {

    }

    virtual ~GadgetronClientBlobAttribMessageReader() {}

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

        std::string meta_attrib(meta_attrib_length-sizeof(size_t_type),0);
        boost::asio::read(*socket, boost::asio::buffer(const_cast<char*>(meta_attrib.c_str()), meta_attrib_length-sizeof(size_t_type)));


        std::string filename_image, filename_attrib;

        // Create the filename: (prefix_%06.suffix)
        if ( file_prefix.empty() )
        {
            filename_image =  filenameBuf + "." + file_suffix;
            filename_attrib =  filenameBuf + "_attrib.xml";
        }
        else
        {
            filename_image = file_prefix + "_" + filenameBuf + "." + file_suffix;
            filename_attrib = file_prefix + "_" + filenameBuf + "_attrib.xml";
        }

        std::cout << "Writing image " << filename_image.c_str() << std::endl;

        std::ofstream outfile;
        outfile.open (filename_image.c_str(), std::ios::out|std::ios::binary);

        std::ofstream outfile_attrib;
        outfile_attrib.open (filename_attrib.c_str(), std::ios::out|std::ios::binary);

        if (outfile.good())
        {
            /* write 'size' bytes starting at 'data's pointer */
            outfile.write(&data[0], nbytes);
            outfile.close();

            outfile_attrib.write(meta_attrib.c_str(), meta_attrib.length());
            outfile_attrib.close();

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
    {

    }

    virtual ~GadgetronClientConnector() 
    {
        if (socket_) {
            socket_->close();
            delete socket_;
        }
    }

    void read_task()
    {
        if (!socket_) {
            throw GadgetronClientException("Unable to create socket.");
        }

        GadgetMessageIdentifier id;
        while (socket_->is_open()) {
            boost::asio::read(*socket_, boost::asio::buffer(&id,sizeof(GadgetMessageIdentifier)));

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
        tcp::resolver::query query(tcp::v4(), hostname.c_str(), port.c_str());
        tcp::resolver::iterator endpoint_iterator = resolver.resolve(query);
        tcp::resolver::iterator end;

        socket_ = new tcp::socket(io_service);

        if (!socket_) {
            throw GadgetronClientException("Unable to create socket.");
        }

        //TODO:
        //For newer versions of Boost, we should use
        //   boost::asio::connect(*socket_, iterator);

        boost::system::error_code error = boost::asio::error::host_not_found;
        while (error && endpoint_iterator != end) {
            socket_->close();
            socket_->connect(*endpoint_iterator++, error);
        }
        if (error)
            throw GadgetronClientException("Error connecting using socket.");

        reader_thread_ = boost::thread(boost::bind(&GadgetronClientConnector::read_task, this));

    }

    void send_gadgetron_close() { 
        if (!socket_) {
            throw GadgetronClientException("Invalid socket.");
        }
        GadgetMessageIdentifier id;
        id.id = GADGET_MESSAGE_CLOSE;    
        boost::asio::write(*socket_, boost::asio::buffer(&id, sizeof(GadgetMessageIdentifier)));
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

        boost::asio::write(*socket_, boost::asio::buffer(&id, sizeof(GadgetMessageIdentifier)));
        boost::asio::write(*socket_, boost::asio::buffer(&ini, sizeof(GadgetMessageConfigurationFile)));

    }

    void send_gadgetron_configuration_script(std::string xml_string)
    {
        if (!socket_) {
            throw GadgetronClientException("Invalid socket.");
        }

        GadgetMessageIdentifier id;
        id.id = GADGET_MESSAGE_CONFIG_SCRIPT;

        GadgetMessageScript conf;
        conf.script_length = (uint32_t)xml_string.size()+1;

        boost::asio::write(*socket_, boost::asio::buffer(&id, sizeof(GadgetMessageIdentifier)));
        boost::asio::write(*socket_, boost::asio::buffer(&conf, sizeof(GadgetMessageScript)));
        boost::asio::write(*socket_, boost::asio::buffer(xml_string.c_str(), conf.script_length));    

    }


    void  send_gadgetron_parameters(std::string xml_string)
    {
        if (!socket_) {
            throw GadgetronClientException("Invalid socket.");
        }

        GadgetMessageIdentifier id;
        id.id = GADGET_MESSAGE_PARAMETER_SCRIPT;

        GadgetMessageScript conf;
        conf.script_length = (uint32_t)xml_string.size()+1;

        boost::asio::write(*socket_, boost::asio::buffer(&id, sizeof(GadgetMessageIdentifier)));
        boost::asio::write(*socket_, boost::asio::buffer(&conf, sizeof(GadgetMessageScript)));
        boost::asio::write(*socket_, boost::asio::buffer(xml_string.c_str(), conf.script_length));    
    }

    void send_ismrmrd_acquisition(ISMRMRD::Acquisition& acq) 
    {
        if (!socket_) {
            throw GadgetronClientException("Invalid socket.");
        }

        GadgetMessageIdentifier id;
        id.id = GADGET_MESSAGE_ISMRMRD_ACQUISITION;;

        boost::asio::write(*socket_, boost::asio::buffer(&id, sizeof(GadgetMessageIdentifier)));
        boost::asio::write(*socket_, boost::asio::buffer(&acq.getHead(), sizeof(ISMRMRD::AcquisitionHeader)));

        unsigned long trajectory_elements = acq.getHead().trajectory_dimensions*acq.getHead().number_of_samples;
        unsigned long data_elements = acq.getHead().active_channels*acq.getHead().number_of_samples;

        if (trajectory_elements) {
            boost::asio::write(*socket_, boost::asio::buffer(&acq.getTrajPtr()[0], sizeof(float)*trajectory_elements));
        }


        if (data_elements) {
            boost::asio::write(*socket_, boost::asio::buffer(&acq.getDataPtr()[0], 2*sizeof(float)*data_elements));
        }
    }

    void register_reader(unsigned short slot, boost::shared_ptr<GadgetronClientMessageReader> r) {
        readers_[slot] = r;
    }

protected:
    typedef std::map<unsigned short, boost::shared_ptr<GadgetronClientMessageReader> > maptype;

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
    boost::thread reader_thread_;
    maptype readers_;


};


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
    std::string out_fileformat;

    po::options_description desc("Allowed options");

    desc.add_options()
        ("help,h", "produce help message")
        ("port,p", po::value<std::string>(&port)->default_value("9002"), "Port")
        ("address,a", po::value<std::string>(&host_name)->default_value("localhost"), "Address (hostname) of Gadgetron host")
        ("filename,f", po::value<std::string>(&in_filename), "Input file")
        ("outfile,o", po::value<std::string>(&out_filename)->default_value("out.h5"), "Output file")
        ("in-group,g", po::value<std::string>(&hdf5_in_group)->default_value("/dataset"), "Input data group")
        ("out-group,G", po::value<std::string>(&hdf5_out_group)->default_value(get_date_time_string()), "Output group name")  
        ("config,c", po::value<std::string>(&config_file)->default_value("default.xml"), "Configuration file (remote)")
        ("config-local,C", po::value<std::string>(&config_file_local), "Configuration file (local)")
        ("loops,l", po::value<unsigned int>(&loops)->default_value(1), "Loops")
        ("outformat,F", po::value<std::string>(&out_fileformat)->default_value("h5"), "Out format, h5 for hdf5 and hdr for analyze image")
        ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << desc << std::endl;
        return 0;
    }

    if (!vm.count("filename")) {
        std::cout << std::endl << std::endl << "\tYou must supply a filename" << std::endl << std::endl;
        std::cout << desc << std::endl;
        return -1;
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

    std::cout << "Gadgetron ISMRMRD client" << std::endl;

    //Let's check if the files exist:
    std::string hdf5_xml_varname = std::string(hdf5_in_group) + std::string("/xml");
    std::string hdf5_data_varname = std::string(hdf5_in_group) + std::string("/data");


    //TODO:
    // Add check to see if input file exists

    //Let's open the input file
    ISMRMRD::Dataset ismrmrd_dataset(in_filename.c_str(), hdf5_in_group.c_str(), false);
    // Read the header
    std::string xml_config;
    ismrmrd_dataset.readHeader(xml_config);


    std::cout << "  -- host            :      " << host_name << std::endl;
    std::cout << "  -- port            :      " << port << std::endl;
    std::cout << "  -- hdf5 file  in   :      " << in_filename << std::endl;
    std::cout << "  -- hdf5 group in   :      " << hdf5_in_group << std::endl;
    std::cout << "  -- conf            :      " << config_file << std::endl;
    std::cout << "  -- loop            :      " << loops << std::endl;
    std::cout << "  -- hdf5 file out   :      " << out_filename << std::endl;
    std::cout << "  -- hdf5 group out  :      " << hdf5_out_group << std::endl;


    GadgetronClientConnector con;

    if ( out_fileformat == "hdr" )
    {
        con.register_reader(GADGET_MESSAGE_ISMRMRD_IMAGE_REAL_USHORT, boost::shared_ptr<GadgetronClientMessageReader>(new GadgetronClientAnalyzeImageMessageReader<uint16_t>(hdf5_out_group)));
        con.register_reader(GADGET_MESSAGE_ISMRMRD_IMAGE_REAL_FLOAT, boost::shared_ptr<GadgetronClientMessageReader>(new GadgetronClientAnalyzeImageMessageReader<float>(hdf5_out_group)));
        con.register_reader(GADGET_MESSAGE_ISMRMRD_IMAGE_CPLX_FLOAT, boost::shared_ptr<GadgetronClientMessageReader>(new GadgetronClientAnalyzeImageMessageReader< std::complex<float> >(hdf5_out_group)));

        //Image with attributes 
        con.register_reader(GADGET_MESSAGE_ISMRMRD_IMAGEWITHATTRIB_REAL_USHORT, boost::shared_ptr<GadgetronClientMessageReader>(new GadgetronClientAttribAnalyzeImageMessageReader<uint16_t>(hdf5_out_group)));
        con.register_reader(GADGET_MESSAGE_ISMRMRD_IMAGEWITHATTRIB_REAL_FLOAT, boost::shared_ptr<GadgetronClientMessageReader>(new GadgetronClientAttribAnalyzeImageMessageReader<float>(hdf5_out_group)));
        con.register_reader(GADGET_MESSAGE_ISMRMRD_IMAGEWITHATTRIB_CPLX_FLOAT, boost::shared_ptr<GadgetronClientMessageReader>(new GadgetronClientAttribAnalyzeImageMessageReader< std::complex<float> >(hdf5_out_group)));
    }
    else
    {
        con.register_reader(GADGET_MESSAGE_ISMRMRD_IMAGE_REAL_USHORT, boost::shared_ptr<GadgetronClientMessageReader>(new GadgetronClientImageMessageReader<uint16_t>(out_filename, hdf5_out_group)));
        con.register_reader(GADGET_MESSAGE_ISMRMRD_IMAGE_REAL_FLOAT, boost::shared_ptr<GadgetronClientMessageReader>(new GadgetronClientImageMessageReader<float>(out_filename, hdf5_out_group)));
        con.register_reader(GADGET_MESSAGE_ISMRMRD_IMAGE_CPLX_FLOAT, boost::shared_ptr<GadgetronClientMessageReader>(new GadgetronClientImageMessageReader< std::complex<float> >(out_filename, hdf5_out_group)));

        //Image with attributes 
        con.register_reader(GADGET_MESSAGE_ISMRMRD_IMAGEWITHATTRIB_REAL_USHORT, boost::shared_ptr<GadgetronClientMessageReader>(new GadgetronClientAttribImageMessageReader<uint16_t>(out_filename, hdf5_out_group)));
        con.register_reader(GADGET_MESSAGE_ISMRMRD_IMAGEWITHATTRIB_REAL_FLOAT, boost::shared_ptr<GadgetronClientMessageReader>(new GadgetronClientAttribImageMessageReader<float>(out_filename, hdf5_out_group)));
        con.register_reader(GADGET_MESSAGE_ISMRMRD_IMAGEWITHATTRIB_CPLX_FLOAT, boost::shared_ptr<GadgetronClientMessageReader>(new GadgetronClientAttribImageMessageReader< std::complex<float> >(out_filename, hdf5_out_group)));
    }

    con.register_reader(GADGET_MESSAGE_DICOM, boost::shared_ptr<GadgetronClientMessageReader>(new GadgetronClientBlobMessageReader(std::string(hdf5_out_group), std::string("dcm"))));
    con.register_reader(GADGET_MESSAGE_DICOM_WITHNAME, boost::shared_ptr<GadgetronClientMessageReader>(new GadgetronClientBlobAttribMessageReader(std::string(), std::string("dcm"))));

    try {
        con.connect(host_name,port);
        if (vm.count("config-local")) {
            con.send_gadgetron_configuration_script(config_xml_local);
        } else {
            con.send_gadgetron_configuration_file(config_file);
        }
        con.send_gadgetron_parameters(xml_config);

        uint32_t acquisitions = 0;
        {
            mtx.lock();
            acquisitions = ismrmrd_dataset.getNumberOfAcquisitions();
            mtx.unlock();
        }

        ISMRMRD::Acquisition acq_tmp;
        for (uint32_t i = 0; i < acquisitions; i++) {
            {
                {
                    boost::mutex::scoped_lock scoped_lock(mtx);
                    ismrmrd_dataset.readAcquisition(i, acq_tmp);
                }
                con.send_ismrmrd_acquisition(acq_tmp);
            }
        }

        con.send_gadgetron_close();
        con.wait();

    } catch (std::exception& ex) {
        std::cout << "Error caught: " << ex.what() << std::endl;
    }

    return 0;
}
