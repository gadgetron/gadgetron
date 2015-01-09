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

class EXPORTGTPLUSIO gtPlusIOAnalyze : public gtPlusIOBase<dsr>
{
public:

    typedef gtPlusIOBase<dsr> BaseClass;
    typedef BaseClass::THeaderType HeaderType;

    gtPlusIOAnalyze();
    gtPlusIOAnalyze(float px, float py);
    gtPlusIOAnalyze(float px, float py, float pz);
    gtPlusIOAnalyze(float px, float py, float pz, float pt);
    gtPlusIOAnalyze(float px, float py, float pz, float pt, float pr);
    gtPlusIOAnalyze(float px, float py, float pz, float pt, float pr, float ps);
    gtPlusIOAnalyze(float px, float py, float pz, float pt, float pr, float ps, float pp);
    gtPlusIOAnalyze(float px, float py, float pz, float pt, float pr, float ps, float pp, float pq);

    virtual ~gtPlusIOAnalyze() {}

    virtual void printInfo(std::ostream& os);

    virtual bool exportArray(const hoNDArray<short>& a, const std::string& filename) { return this->exportArrayImpl(a, filename); }

    virtual bool exportArray(const hoNDArray<unsigned short>& a, const std::string& filename) { return this->exportArrayImpl(a, filename); }
    virtual bool exportArray(const hoNDArray<int>& a, const std::string& filename) { return this->exportArrayImpl(a, filename); }
    virtual bool exportArray(const hoNDArray<unsigned int>& a, const std::string& filename) { return this->exportArrayImpl(a, filename); }
    virtual bool exportArray(const hoNDArray<size_t>& a, const std::string& filename) { return this->exportArrayImpl(a, filename); }
    virtual bool exportArray(const hoNDArray<float>& a, const std::string& filename) { return this->exportArrayImpl(a, filename); }
    virtual bool exportArray(const hoNDArray<double>& a, const std::string& filename) { return this->exportArrayImpl(a, filename); }
    virtual bool exportArray(const hoNDArray< std::complex<float> >& a, const std::string& filename) { return this->exportArrayImpl(a, filename); }
    virtual bool exportArray(const hoNDArray< std::complex<double> >& a, const std::string& filename) { return this->exportArrayImpl(a, filename); }

    virtual bool importArray(hoNDArray<short>& a, const std::string& filename) { return this->importArrayImpl(a, filename); }
    virtual bool importArray(hoNDArray<unsigned short>& a, const std::string& filename) { return this->importArrayImpl(a, filename); }
    virtual bool importArray(hoNDArray<int>& a, const std::string& filename) { return this->importArrayImpl(a, filename); }
    virtual bool importArray(hoNDArray<unsigned int>& a, const std::string& filename) { return this->importArrayImpl(a, filename); }
    virtual bool importArray(hoNDArray<float>& a, const std::string& filename) { return this->importArrayImpl(a, filename); }
    virtual bool importArray(hoNDArray<double>& a, const std::string& filename) { return this->importArrayImpl(a, filename); }
    virtual bool importArray(hoNDArray< std::complex<float> >& a, const std::string& filename) { return this->importArrayImpl(a, filename); }
    virtual bool importArray(hoNDArray< std::complex<double> >& a, const std::string& filename) { return this->importArrayImpl(a, filename); }

    template <typename T> 
    bool exportArrayImpl(const hoNDArray<T>& a, const std::string& filename)
    {
        try
        {
            HeaderType header;
            GADGET_CHECK_RETURN_FALSE(this->array2Header(a, header));
            GADGET_CHECK_RETURN_FALSE(this->writeHeader(filename, header));

            std::string filenameData = filename;
            filenameData.append(".img");
            GADGET_CHECK_RETURN_FALSE(this->writeData(filenameData, a.begin(), a.get_number_of_bytes()));
        }
        catch(...)
        {
            GERROR_STREAM("Errors in gtPlusIOAnalyze::exportArrayImpl(const hoNDArray<T>& a, const std::string& filename) ... ");
            return false;
        }

        return true;
    }

    template <typename T> 
    bool importArrayImpl(hoNDArray<T>& a, const std::string& filename)
    {
        try
        {
            HeaderType header;
            GADGET_CHECK_RETURN_FALSE(this->readHeader(filename, header));
            GADGET_CHECK_RETURN_FALSE(this->header2Array(a, header));

            std::string filenameData = filename;
            filenameData.append(".img");
            GADGET_CHECK_RETURN_FALSE(this->readData(filenameData, a.begin(), a.get_number_of_bytes()));
        }
        catch(...)
        {
            GERROR_STREAM("Errors in gtPlusIOAnalyze::importArrayImpl(const hoNDArray<T>& a, const std::string& filename) ... ");
            return false;
        }

        return true;
    }

    template <typename T, unsigned int D> 
    bool exportImage(const hoNDImage<T,D>& a, const std::string& filename)
    {
        try
        {
            HeaderType header;
            GADGET_CHECK_RETURN_FALSE(this->image2Header(a, header));
            GADGET_CHECK_RETURN_FALSE(this->writeHeader(filename, header));

            std::string filenameData = filename;
            filenameData.append(".img");
            GADGET_CHECK_RETURN_FALSE(this->writeData(filenameData, a.begin(), a.get_number_of_bytes()));
        }
        catch(...)
        {
            GERROR_STREAM("Errors in exportImage(const hoNDImage<T,D>& a, const std::string& filename) ... ");
            return false;
        }

        return true;
    }

    template <typename T, unsigned int D> 
    bool importImage(hoNDImage<T,D>& a, const std::string& filename)
    {
        try
        {
            HeaderType header;
            GADGET_CHECK_RETURN_FALSE(this->readHeader(filename, header));
            GADGET_CHECK_RETURN_FALSE(this->header2Image(a, header));

            std::string filenameData = filename;
            filenameData.append(".img");
            GADGET_CHECK_RETURN_FALSE(this->readData(filenameData, a.begin(), a.get_number_of_bytes()));
        }
        catch(...)
        {
            GERROR_STREAM("Errors in importImage(const hoNDImage<T,D>& a, const std::string& filename) ... ");
            return false;
        }

        return true;
    }

/// image functions

    template <typename T, unsigned int D> 
    bool exportImageComplex(const hoNDImage<T,D>& a, const std::string& filename)
    {
        try
        {
            typedef typename Gadgetron::realType<T>::Type value_type;

            //hoNDImage<value_type, D> buf;
            //GADGET_CHECK_RETURN_FALSE(Gadgetron::complex_to_real(a, buf));

            //std::string filenameReal = filename;
            //filenameReal.append("_REAL");
            //GADGET_CHECK_RETURN_FALSE(exportImage(buf, filenameReal));

            //GADGET_CHECK_RETURN_FALSE(Gadgetron::complex_to_imag(a, buf));
            //std::string filenameImag = filename;
            //filenameImag.append("_IMAG");
            //GADGET_CHECK_RETURN_FALSE(exportImage(buf, filenameImag));

            //GADGET_CHECK_RETURN_FALSE(Gadgetron::abs(a, buf));
            //std::string filenameMag = filename;
            //filenameMag.append("_MAG");
            //GADGET_CHECK_RETURN_FALSE(exportImage(buf, filenameMag));

            //GADGET_CHECK_RETURN_FALSE(Gadgetron::argument(a, buf));
            //std::string filenamePhase = filename;
            //filenamePhase.append("_PHASE");
            //GADGET_CHECK_RETURN_FALSE(exportImage(buf, filenamePhase));

            long long num = (long long)a.get_number_of_elements();

            long long n;

            hoNDImage<value_type, D> rpart, ipart, mag, phs;
            rpart.create( *a.get_dimensions() );
            ipart.create( *a.get_dimensions() );
            mag.create( *a.get_dimensions() );
            phs.create( *a.get_dimensions() );

            const T* pA = a.begin();

            #pragma omp parallel for default(none) private(n) shared(num, pA, rpart, ipart, mag, phs)
            for ( n=0; n<num; n++ )
            {
                rpart(n) = pA[n].real();
                ipart(n) = pA[n].imag();
                mag(n) = std::abs( pA[n] );
                phs(n) = std::arg( pA[n] );
            }

            std::string filenameReal = filename;
            filenameReal.append("_REAL");
            GADGET_CHECK_RETURN_FALSE(exportImage(rpart, filenameReal));

            std::string filenameImag = filename;
            filenameImag.append("_IMAG");
            GADGET_CHECK_RETURN_FALSE(exportImage(ipart, filenameImag));

            std::string filenameMag = filename;
            filenameMag.append("_MAG");
            GADGET_CHECK_RETURN_FALSE(exportImage(mag, filenameMag));

            std::string filenamePhase = filename;
            filenamePhase.append("_PHASE");
            GADGET_CHECK_RETURN_FALSE(exportImage(phs, filenamePhase));
        }
        catch(...)
        {
            GERROR_STREAM("Errors in exportImageComplex(const hoNDImage<T,D>& a, const std::string& filename) ... ");
            return false;
        }

        return true;
    }

    template <typename T, unsigned int D> 
    bool importImageComplex(hoNDImage<T,D>& a, const std::string& filename)
    {
        try
        {
            typedef typename T::value_type value_type;
            hoNDImage<value_type, D> real, imag;

            std::string filenameReal = filename;
            filenameReal.append("_REAL");
            GADGET_CHECK_RETURN_FALSE(importImage(real, filenameReal));

            std::string filenameImag = filename;
            filenameImag.append("_IMAG");
            GADGET_CHECK_RETURN_FALSE(importImage(imag, filenameImag));

            a.create(real.get_dimensions());
            long long num = (long long)a.get_number_of_elements();

            long long n;
            T* pA = a.begin();

            #pragma omp parallel for default(none) private(n) shared(num, pA, real, imag)
            for ( n=0; n<num; n++ )
            {
                pA[n] = T( real(n), imag(n) );
            }
        }
        catch(...)
        {
            GERROR_STREAM("Errors in importImageComplex(const hoNDImage<T,D>& a, const std::string& filename) ... ");
            return false;
        }

        return true;
    }

    template <typename T, unsigned int D> 
    bool importImageComplex(hoNDImage<T,D>& a, const std::string& filename_real, const std::string& filename_imag)
    {
        try
        {
            typedef typename realType<T>::Type value_type;
            hoNDImage<value_type, D> real, imag;

            GADGET_CHECK_RETURN_FALSE(importImage(real, filename_real));
            GADGET_CHECK_RETURN_FALSE(importImage(imag, filename_imag));

            a.create(real.get_dimensions());
            long long num = (long long)a.get_number_of_elements();

            long long n;
            T* pA = a.begin();

            #pragma omp parallel for default(none) private(n) shared(num, pA, real, imag)
            for ( n=0; n<num; n++ )
            {
                pA[n] = T( real(n), imag(n) );
            }
        }
        catch(...)
        {
            GERROR_STREAM("Errors in importImageComplex(hoNDImage<T,D>& a, const std::string& filename_real, const std::string& filename_imag) ... ");
            return false;
        }

        return true;
    }

    template <typename T> 
    bool export2DImage(const hoNDImage<T,2>& a, const std::string& filename)
    {
        return exportImage(a, filename);
    }

    template <typename T> 
    bool import2DImage(hoNDImage<T,2>& a, const std::string& filename)
    {
        return importImage(a, filename);
    }

    template <typename T> 
    bool export2DImageComplex(const hoNDImage<T,2>& a, const std::string& filename)
    {
        return exportImageComplex(a, filename);
    }

    template <typename T> 
    bool import2DImageComplex(hoNDImage<T,2>& a, const std::string& filename)
    {
        return importImageComplex(a, filename);
    }

    template <typename T> 
    bool export3DImage(const hoNDImage<T,3>& a, const std::string& filename)
    {
        return exportImage(a, filename);
    }

    template <typename T> 
    bool import3DImage(hoNDImage<T,3>& a, const std::string& filename)
    {
        return importImage(a, filename);
    }

    template <typename T> 
    bool export3DImageComplex(const hoNDImage<T,3>& a, const std::string& filename)
    {
        return exportImageComplex(a, filename);
    }

    template <typename T> 
    bool import3DImageComplex(hoNDImage<T,3>& a, const std::string& filename)
    {
        return importImageComplex(a, filename);
    }

    template <typename T> 
    bool export4DImage(const hoNDImage<T,4>& a, const std::string& filename)
    {
        try
        {
            size_t RO     = a.get_size(0);
            size_t E1     = a.get_size(1);
            size_t CHA    = a.get_size(2);
            size_t N      = a.get_size(3);

            size_t ii;
            for (ii=0; ii<N; ii++ )
            {
                std::vector<size_t> dim(3);
                dim[0] = RO;
                dim[1] = E1;
                dim[2] = CHA;

                hoNDImage<T, 3> a3D(dim, const_cast<T*>(a.begin()+ii*RO*E1*CHA), false);

                std::ostringstream ostr;
                ostr << filename << "_" << ii << std::ends;
                GADGET_CHECK_RETURN_FALSE(export3DImage(a3D, ostr.str()));
            }
        }
        catch(...)
        {
            GERROR_STREAM("Errors in export4DImage(const hoNDImage<T>& a, const std::string& filename) ... ");
            return false;
        }

        return true;
    }

    template <typename T> 
    bool export4DImageComplex(const hoNDImage<T,4>& a, const std::string& filename)
    {
        try
        {
            size_t RO     = a.get_size(0);
            size_t E1     = a.get_size(1);
            size_t CHA    = a.get_size(2);
            size_t N      = a.get_size(3);

            size_t ii;
            for (ii=0; ii<N; ii++ )
            {
                std::vector<size_t> dim(3);
                dim[0] = RO;
                dim[1] = E1;
                dim[2] = CHA;

                hoNDImage<T, 3> a3D(dim, const_cast<T*>(a.begin()+ii*RO*E1*CHA), false);

                std::ostringstream ostr;
                ostr << filename << "_" << ii << std::ends;
                GADGET_CHECK_RETURN_FALSE(export3DImageComplex(a3D, ostr.str()));
            }
        }
        catch(...)
        {
            GERROR_STREAM("Errors in export4DImageComplex(const hoNDImage<T>& a, const std::string& filename) ... ");
            return false;
        }

        return true;
    }

protected:

    template <typename T> bool array2Header(const hoNDArray<T>& a, HeaderType& header);
    template <typename T> bool header2Array(hoNDArray<T>& a, const HeaderType& header);

    template <typename T, unsigned int D> bool image2Header(const hoNDImage<T, D>& a, HeaderType& header);
    template <typename T, unsigned int D> bool header2Image(hoNDImage<T, D>& a, const HeaderType& header);

    // read/write the analyze header
    bool readHeader(const std::string& filename, HeaderType& header);
    bool writeHeader(const std::string& filename, const HeaderType& header);
};

template <typename T> 
bool gtPlusIOAnalyze::array2Header(const hoNDArray<T>& a, HeaderType& header)
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
        size_t NDim = a.get_number_of_dimensions();

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
        header.dime.datatype = (short)getDataTypeFromRTTI(rttiID);
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
        GERROR_STREAM("Errors in gtPlusIOAnalyze::array2Analyze(const hoNDArray<T>& a, dsr& header) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusIOAnalyze::header2Array(hoNDArray<T>& a, const HeaderType& header)
{
    try
    {
        std::string rttiID = std::string(typeid(T).name());
        GADGET_CHECK_RETURN_FALSE(rttiID==getRTTIFromDataType( (GtDataType)header.dime.datatype));

        std::vector<size_t> dim(header.dime.dim[0]);
        size_t ii;
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
        GERROR_STREAM("Errors in gtPlusIOAnalyze::analyze2Array(hoNDArray<T>& a, const dsr& header) ... ");
        return false;
    }

    return true;
}

template <typename T, unsigned int D> 
bool gtPlusIOAnalyze::image2Header(const hoNDImage<T,D>& a, HeaderType& header)
{
    try
    {
        typedef typename hoNDImage<T,D>::coord_type coord_type;

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
        size_t NDim = D;

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
        header.dime.datatype = (short)getDataTypeFromRTTI(rttiID);
        header.dime.bitpix = (short)(8*sizeof(T));
        header.dime.dim_un0 = 0;

        header.dime.pixdim[0] = 0;
        header.dime.pixdim[1] = a.get_pixel_size(0);
        header.dime.pixdim[2] = 1;
        header.dime.pixdim[3] = 1;
        if ( NDim > 1 )
            header.dime.pixdim[2] = a.get_pixel_size(1);
        if ( NDim > 2 )
            header.dime.pixdim[3] = a.get_pixel_size(2);
        if ( NDim > 3 )
            header.dime.pixdim[4] = a.get_pixel_size(3);
        if ( NDim > 4 )
            header.dime.pixdim[5] = a.get_pixel_size(4);
        if ( NDim > 5 )
            header.dime.pixdim[6] = a.get_pixel_size(5);
        if ( NDim > 6 )
            header.dime.pixdim[7] = a.get_pixel_size(6);

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

        // store image origin and axis
        // total number of bytes are
        size_t numOfBytes = sizeof(float)*(D+D*D);
        if ( numOfBytes <= sizeof(data_history) )
        {
            std::vector<float> buf(D+D*D, 0);

            unsigned int ii;
            for ( ii=0; ii<D; ii++ )
            {
                buf[ii] = (float)a.get_origin(ii);
            }

            unsigned int jj;
            for ( ii=0; ii<D; ii++ )
            {
                for ( jj=0; jj<D; jj++ )
                {
                    buf[D+ii*D+jj] = (float)a.get_axis(ii, jj);
                }
            }

            memcpy(&header.hist, &buf[0], numOfBytes);
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusIOAnalyze::image2Analyze(const hoNDImage<T>& a, dsr& header) ... ");
        return false;
    }

    return true;
}

template <typename T, unsigned int D> 
bool gtPlusIOAnalyze::header2Image(hoNDImage<T,D>& a, const HeaderType& header)
{
    try
    {
        std::string rttiID = std::string(typeid(T).name());
        GADGET_CHECK_RETURN_FALSE(rttiID==getRTTIFromDataType( (GtDataType)header.dime.datatype));

        std::vector<size_t> dim(header.dime.dim[0]);

        if ( D > dim.size() ) return false;

        size_t ii;
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

        a.create(dim);

        for ( ii=0; ii<dim.size(); ii++ )
        {
            if ( ii < 7 )
            {
                a.set_pixel_size(ii, header.dime.pixdim[ii+1]);
            }
        }

        // get origin and axis
        size_t numOfBytes = sizeof(float)*(D+D*D);
        if ( numOfBytes <= sizeof(data_history) )
        {
            std::vector<float> buf(D+D*D);
            memcpy(&buf[0], &header.hist, numOfBytes);

            unsigned int ii;
            for ( ii=0; ii<D; ii++ )
            {
                a.set_origin(ii, buf[ii]);
            }

            unsigned int jj;
            for ( ii=0; ii<D; ii++ )
            {
                typename hoNDImage<T,D>::coord_type v(0);
                typename hoNDImage<T,D>::coord_type mag(0);

                for ( jj=0; jj<D; jj++ )
                {
                    v = buf[D+ii*D+jj];
                    mag += v*v;
                    a.set_axis(ii, jj, v);
                }

                if ( mag < FLT_EPSILON )
                {
                    for ( jj=0; jj<D; jj++ )
                    {
                        if ( ii != jj )
                        {
                            a.set_axis(ii, jj, 0);
                        }
                        else
                        {
                            a.set_axis(ii, jj, (typename hoNDImage<T,D>::coord_type)(1.0) );
                        }
                    }
                }
            }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusIOAnalyze::analyze2Image(hoNDImage<T,D>& a, const dsr& header) ... ");
        return false;
    }

    return true;
}

}}
