/** \file       gtPlusIOAnalyze.h
    \brief      Implement the suppor for the Analzye75 medical image format
    \author     Hui Xue

    The ISMRMRD dimensions are mapped to Analyze75 format.

    Ref to:
    http://eeg.sourceforge.net/ANALYZE75.pdf
    http://ismrmrd.sourceforge.net/
*/

#pragma once

#include "gtPlusIOBase.h"

// the file input/output utility functions for the Analyze format

// the following Analyze75 data structured is defined as this online document eeg.sourceforge.net/ANALYZE75.pdf‎

enum AnalyzeDataType
{
    DT_ANA_UNKNOWN=0,
    DT_BINARY=1, 
    DT_UNSIGNED_CHAR=2,
    DT_SIGNED_SHORT=4,
    DT_UNSIGNED_SHORT=5,
    DT_SIGNED_INT=8,
    DT_UNSIGNED_INT=9,
    DT_FLOAT=16,
    DT_COMPLEX=32,
    DT_DOUBLE=64,
    DT_DOUBLECOMPLEX=96, // this type is added to support complex doulbe
    DT_RGB=128,
    DT_ALL=255
};

// the official definition of Analyze 7.5 file format
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

// to suppor the ISMRMRD format
// [Ro E1 Cha Slice E2 Con Phase Rep Set Seg]

namespace Gadgetron { namespace gtPlus {

class EXPORTGTPLUS gtPlusIOAnalyze
{
public:

    gtPlusIOAnalyze() { pixelSize_.resize(10, 1.0); }
    gtPlusIOAnalyze(float px, float py);
    gtPlusIOAnalyze(float px, float py, float pz);
    gtPlusIOAnalyze(float px, float py, float pz, float pt);
    gtPlusIOAnalyze(float px, float py, float pz, float pt, float pr);
    gtPlusIOAnalyze(float px, float py, float pz, float pt, float pr, float ps);
    gtPlusIOAnalyze(float px, float py, float pz, float pt, float pr, float ps, float pp);
    gtPlusIOAnalyze(float px, float py, float pz, float pt, float pr, float ps, float pp, float pq);

    void setPixelSize(float px, float py, float pz=1.0f, float pt=1.0f, float pr=1.0f, float ps=1.0f, float pp=1.0f, float pq=1.0f);

    virtual ~gtPlusIOAnalyze() {}

public:

    void printInfo(std::ostream& os);

    // export/input for 2D/3D/4D array
    // filename should be given without .hdr extension
    // the .hdr and .img extension will be added internally

    template <typename T> bool exportArray(const hoNDArray<T>& a, const std::string& filename);
    template <typename T> bool importArray(hoNDArray<T>& a, const std::string& filename);

    template <typename T> bool exportArrayComplex(const hoNDArray<T>& a, const std::string& filename);
    template <typename T> bool importArrayComplex(hoNDArray<T>& a, const std::string& filename);

    template <typename T> bool importArrayComplex(hoNDArray<T>& a, const std::string& filename_real, const std::string& filename_imag);

    // 2D array is exported as a 2D image
    template <typename T> bool export2DArray(const hoNDArray<T>& a, const std::string& filename);
    template <typename T> bool import2DArray(hoNDArray<T>& a, const std::string& filename);

    template <typename T> bool export2DArrayComplex(const hoNDArray<T>& a, const std::string& filename);
    template <typename T> bool import2DArrayComplex(hoNDArray<T>& a, const std::string& filename);

    // 3D array is exported as a 3D volume
    template <typename T> bool export3DArray(const hoNDArray<T>& a, const std::string& filename);
    template <typename T> bool import3DArray(hoNDArray<T>& a, const std::string& filename);

    template <typename T> bool export3DArrayComplex(const hoNDArray<T>& a, const std::string& filename);
    template <typename T> bool import3DArrayComplex(hoNDArray<T>& a, const std::string& filename);

    // 4D array is exported as multiple 3D volume
    template <typename T> bool export4DArray(const hoNDArray<T>& a, const std::string& filename);
    template <typename T> bool export4DArrayComplex(const hoNDArray<T>& a, const std::string& filename);

protected:

    std::vector<float> pixelSize_;

    template <typename T> bool array2Analyze(const hoNDArray<T>& a, dsr& header);
    template <typename T> bool analyze2Array(hoNDArray<T>& a, const dsr& header);

    // get the run-time type ID from analyze data type or vice versa
    std::string getRTTIFromAnalyzeDataType(AnalyzeDataType aDT);
    AnalyzeDataType getAnalyzeDataTypeFromRTTI(const std::string& name);

    // read/write the analyze header
    bool readAnalyzeHeader(const std::string& filename, dsr& header);
    bool writeAnalyzeHeader(const std::string& filename, const dsr& header);

    // read/write the analyze data file
    // len is the number of bytes
    template <typename T> bool readAnalyzeData(const std::string& filename, T* data, long long len);
    template <typename T> bool writeAnalyzeData(const std::string& filename, const T* data, long long len);
};

template <typename T> 
bool gtPlusIOAnalyze::exportArray(const hoNDArray<T>& a, const std::string& filename)
{
    try
    {
        dsr header;
        GADGET_CHECK_RETURN_FALSE(array2Analyze(a, header));
        GADGET_CHECK_RETURN_FALSE(writeAnalyzeHeader(filename, header));
        GADGET_CHECK_RETURN_FALSE(writeAnalyzeData(filename, a.begin(), a.get_number_of_bytes()));
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in gtPlusIOAnalyze::exportArray(const hoNDArray<T>& a, const std::string& filename) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusIOAnalyze::importArray(hoNDArray<T>& a, const std::string& filename)
{
    try
    {
        dsr header;
        GADGET_CHECK_RETURN_FALSE(readAnalyzeHeader(filename, header));
        GADGET_CHECK_RETURN_FALSE(analyze2Array(a, header));
        GADGET_CHECK_RETURN_FALSE(readAnalyzeData(filename, a.begin(), a.get_number_of_bytes()));
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in gtPlusIOAnalyze::importArray(const hoNDArray<T>& a, const std::string& filename) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusIOAnalyze::exportArrayComplex(const hoNDArray<T>& a, const std::string& filename)
{
    try
    {
        typedef typename Gadgetron::realType<T>::Type value_type;

        hoNDArray<value_type> buf;
        GADGET_CHECK_RETURN_FALSE(Gadgetron::complex_to_real(a, buf));

        std::string filenameReal = filename;
        filenameReal.append("_REAL");
        GADGET_CHECK_RETURN_FALSE(exportArray(buf, filenameReal));

        GADGET_CHECK_RETURN_FALSE(Gadgetron::complex_to_imag(a, buf));
        std::string filenameImag = filename;
        filenameImag.append("_IMAG");
        GADGET_CHECK_RETURN_FALSE(exportArray(buf, filenameImag));

        GADGET_CHECK_RETURN_FALSE(Gadgetron::absolute(a, buf));
        std::string filenameMag = filename;
        filenameMag.append("_MAG");
        GADGET_CHECK_RETURN_FALSE(exportArray(buf, filenameMag));

        GADGET_CHECK_RETURN_FALSE(Gadgetron::argument(a, buf));
        std::string filenamePhase = filename;
        filenamePhase.append("_PHASE");
        GADGET_CHECK_RETURN_FALSE(exportArray(buf, filenamePhase));
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in gtPlusIOAnalyze::exportArrayComplex(const hoNDArray<T>& a, const std::string& filename) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusIOAnalyze::importArrayComplex(hoNDArray<T>& a, const std::string& filename)
{
    try
    {
        typedef typename T::value_type value_type;
        hoNDArray<value_type> real, imag;

        std::string filenameReal = filename;
        filenameReal.append("_REAL");
        GADGET_CHECK_RETURN_FALSE(importArray(real, filenameReal));

        std::string filenameImag = filename;
        filenameImag.append("_IMAG");
        GADGET_CHECK_RETURN_FALSE(importArray(imag, filenameImag));

        GADGET_CHECK_RETURN_FALSE(Gadgetron::real_imag_to_complex(real, imag, a));
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in gtPlusIOAnalyze::importArrayComplex(const hoNDArray<T>& a, const std::string& filename) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusIOAnalyze::importArrayComplex(hoNDArray<T>& a, const std::string& filename_real, const std::string& filename_imag)
{
    try
    {
        typedef typename realType<T>::Type value_type;
        hoNDArray<value_type> real, imag;

        GADGET_CHECK_RETURN_FALSE(importArray(real, filename_real));
        GADGET_CHECK_RETURN_FALSE(importArray(imag, filename_imag));

        GADGET_CHECK_RETURN_FALSE(Gadgetron::real_imag_to_complex(real, imag, a));
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in gtPlusIOAnalyze::importArrayComplex(hoNDArray<T>& a, const std::string& filename_real, const std::string& filename_imag) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusIOAnalyze::export2DArray(const hoNDArray<T>& a, const std::string& filename)
{
    return exportArray(a, filename);
}

template <typename T> 
bool gtPlusIOAnalyze::import2DArray(hoNDArray<T>& a, const std::string& filename)
{
    return importArray(a, filename);
}

template <typename T> 
bool gtPlusIOAnalyze::export2DArrayComplex(const hoNDArray<T>& a, const std::string& filename)
{
    return exportArrayComplex(a, filename);
}

template <typename T> 
bool gtPlusIOAnalyze::import2DArrayComplex(hoNDArray<T>& a, const std::string& filename)
{
    return importArrayComplex(a, filename);
}

template <typename T> 
bool gtPlusIOAnalyze::export3DArray(const hoNDArray<T>& a, const std::string& filename)
{
    return exportArray(a, filename);
}

template <typename T> 
bool gtPlusIOAnalyze::import3DArray(hoNDArray<T>& a, const std::string& filename)
{
    return importArray(a, filename);
}

template <typename T> 
bool gtPlusIOAnalyze::export3DArrayComplex(const hoNDArray<T>& a, const std::string& filename)
{
    return exportArrayComplex(a, filename);
}

template <typename T> 
bool gtPlusIOAnalyze::import3DArrayComplex(hoNDArray<T>& a, const std::string& filename)
{
    return importArrayComplex(a, filename);
}

template <typename T> 
bool gtPlusIOAnalyze::export4DArray(const hoNDArray<T>& a, const std::string& filename)
{
    try
    {
        unsigned long long RO     = a.get_size(0);
        unsigned long long E1     = a.get_size(1);
        unsigned long long CHA    = a.get_size(2);
        unsigned long long N      = a.get_size(3);

        unsigned long long ii;
        for (ii=0; ii<N; ii++ )
        {
            std::vector<unsigned long long> dim(3);
            dim[0] = RO;
            dim[1] = E1;
            dim[2] = CHA;

            boost::shared_ptr< std::vector<unsigned long long> > sDim(&dim);
            hoNDArray<T> a3D(sDim, const_cast<T*>(a.begin()+ii*RO*E1*CHA), false);

            std::ostringstream ostr;
            ostr << filename << "_" << ii << std::ends;
            GADGET_CHECK_RETURN_FALSE(export3DArray(a3D, ostr.str()));
        }
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in gtPlusIOAnalyze::export4DArray(const hoNDArray<T>& a, const std::string& filename) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusIOAnalyze::export4DArrayComplex(const hoNDArray<T>& a, const std::string& filename)
{
    try
    {
        unsigned long long RO     = a.get_size(0);
        unsigned long long E1     = a.get_size(1);
        unsigned long long CHA    = a.get_size(2);
        unsigned long long N      = a.get_size(3);

        unsigned long long ii;
        for (ii=0; ii<N; ii++ )
        {
            std::vector<unsigned long long> dim(3);
            dim[0] = RO;
            dim[1] = E1;
            dim[2] = CHA;

            boost::shared_ptr< std::vector<unsigned long long> > sDim(&dim);
            hoNDArray<T> a3D(sDim, const_cast<T*>(a.begin()+ii*RO*E1*CHA), false);

            std::ostringstream ostr;
            ostr << filename << "_" << ii << std::ends;
            GADGET_CHECK_RETURN_FALSE(export3DArrayComplex(a3D, ostr.str()));
        }
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in gtPlusIOAnalyze::export4DArrayComplex(const hoNDArray<T>& a, const std::string& filename) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusIOAnalyze::array2Analyze(const hoNDArray<T>& a, dsr& header)
{
    try
    {
        // set everything to zero
        memset(&header, 0, sizeof(dsr));

        // header_key
        header.hk.sizeof_hdr = 348;
        unsigned long long i;
        for (i=0; i<10; i++ ) header.hk.data_type[i] = 0;
        for (i=0; i<18; i++ ) header.hk.db_name[i] = 0;
        header.hk.extents = 16384;
        header.hk.session_error = 0;
        header.hk.regular = 'r';
        header.hk.hkey_un0 = 0;

        // image_dimension
        unsigned long long NDim = a.get_number_of_dimensions();

        header.dime.dim[0] = (short)(NDim);
        header.dime.dim[1] = (short)(a.get_size(0));

        if ( NDim > 1 )
            header.dime.dim[2] = (short)(a.get_size(1));
        else
            header.dime.dim[2] = 1;

        if ( NDim > 2 )
            header.dime.dim[3] = (short)(a.get_size(2));
        else
            header.dime.dim[3] = 1;

        if ( NDim > 3 )
            header.dime.dim[4] = (short)(a.get_size(3));
        else
            header.dime.dim[4] = 1;

        if ( NDim > 4 )
            header.dime.dim[5] = (short)(a.get_size(4));
        else
            header.dime.dim[5] = 1;

        if ( NDim > 5 )
            header.dime.dim[6] = (short)(a.get_size(5));
        else
            header.dime.dim[6] = 1;

        if ( NDim > 6 )
            header.dime.dim[7] = (short)(a.get_size(6));
        else
            header.dime.dim[7] = 1;

        if ( NDim > 7 )
            header.dime.unused8 = (short)(a.get_size(7));
        else
            header.dime.unused8 = 1;

        if ( NDim > 8 )
            header.dime.unused9 = (short)(a.get_size(8));
        else
            header.dime.unused9 = 1;

        if ( NDim > 9 )
            header.dime.unused10 = (short)(a.get_size(9));
        else
            header.dime.unused10 = 1;

        header.dime.unused11 = 0;
        header.dime.unused12 = 0;
        header.dime.unused13 = 0;
        header.dime.unused14 = 0;

        std::string rttiID = std::string(typeid(T).name());
        header.dime.datatype = (short)getAnalyzeDataTypeFromRTTI(rttiID);
        header.dime.bitpix = (short)(8*sizeof(T));
        header.dime.dim_un0 = 0;

        // since the NDArray does not carry the pixel spacing
        header.dime.pixdim[0] = 0;
        if ( pixelSize_.size() > 1 )
            header.dime.pixdim[1] = pixelSize_[0];
        if ( pixelSize_.size() > 2 )
            header.dime.pixdim[2] = pixelSize_[1];
        if ( pixelSize_.size() > 3 )
            header.dime.pixdim[3] = pixelSize_[2];
        if ( pixelSize_.size() > 4 )
            header.dime.pixdim[4] = pixelSize_[3];
        if ( pixelSize_.size() > 5 )
            header.dime.pixdim[5] = pixelSize_[4];
        if ( pixelSize_.size() > 6 )
            header.dime.pixdim[6] = pixelSize_[5];
        if ( pixelSize_.size() > 7 )
            header.dime.pixdim[7] = pixelSize_[6];

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
        GADGET_ERROR_MSG("Errors in gtPlusIOAnalyze::array2Analyze(const hoNDArray<T>& a, dsr& header) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusIOAnalyze::analyze2Array(hoNDArray<T>& a, const dsr& header)
{
    try
    {
        std::string rttiID = std::string(typeid(T).name());
        GADGET_CHECK_RETURN_FALSE(rttiID==getRTTIFromAnalyzeDataType( (AnalyzeDataType)header.dime.datatype));

        std::vector<unsigned long long> dim(header.dime.dim[0]);
        unsigned long long ii;
        for ( ii=0; ii<dim.size(); ii++ )
        {
            if ( ii == 7 )
            {
                dim[ii] = header.dime.unused8;
            }
            else if ( ii == 8 )
            {
                dim[ii] = header.dime.unused9;
            }
            else if ( ii == 9 ) 
            {
                dim[ii] = header.dime.unused10;
            }
            else
            {
                dim[ii] = header.dime.dim[ii+1];
            }
        }

        pixelSize_.resize(dim.size());
        for ( ii=0; ii<dim.size(); ii++ )
        {
            if ( ii < 7 )
            {
                pixelSize_[ii] = header.dime.pixdim[ii+1];
            }
        }

        a.create(&dim);
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in gtPlusIOAnalyze::analyze2Array(hoNDArray<T>& a, const dsr& header) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusIOAnalyze::readAnalyzeData(const std::string& filename, T* data, long long len)
{
    try
    {
        std::string filenameData = filename;
        filenameData.append(".img");
        gtPlusIOWorker ioworker(filenameData, true);

        GADGET_CHECK_RETURN_FALSE(ioworker.open());
        GADGET_CHECK_RETURN_FALSE(ioworker.read(reinterpret_cast<char*>(data), len));
        GADGET_CHECK_RETURN_FALSE(ioworker.close());
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in gtPlusIOAnalyze::readAnalyzeData(const std::string& filename, T* data, long long len) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusIOAnalyze::writeAnalyzeData(const std::string& filename, const T* data, long long len)
{
    try
    {
        std::string filenameData = filename;
        filenameData.append(".img");
        gtPlusIOWorker ioworker(filenameData, false);

        GADGET_CHECK_RETURN_FALSE(ioworker.open());
        GADGET_CHECK_RETURN_FALSE(ioworker.write(reinterpret_cast<const char*>(data), len));
        GADGET_CHECK_RETURN_FALSE(ioworker.close());
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in gtPlusIOAnalyze::writeAnalyzeData(const std::string& filename, const T* data, long long len) ... ");
        return false;
    }

    return true;
}

}}
