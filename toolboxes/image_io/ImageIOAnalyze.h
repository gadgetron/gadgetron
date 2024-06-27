/** \file       ImageIOAnalyze.h
    \brief      Implement the suppor for the Analzye75 medical image format
    \author     Hui Xue

    The ISMRMRD dimensions are mapped to Analyze75 format.

    Ref to:
    http://eeg.sourceforge.net/ANALYZE75.pdf
    http://ismrmrd.sourceforge.net/
*/

#pragma once

#include "ImageIOBase.h"
#include <filesystem>

// the file input/output utility functions for the Analyze format

// the following Analyze75 data structured is defined as this online document eeg.sourceforge.net/ANALYZE75.pdf

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

namespace Gadgetron { 

class EXPORTIMAGEIO ImageIOAnalyze : public ImageIOBase<dsr>
{
public:

    typedef ImageIOBase<dsr> BaseClass;
    typedef BaseClass::THeaderType HeaderType;

    ImageIOAnalyze();
    ImageIOAnalyze(float px, float py);
    ImageIOAnalyze(float px, float py, float pz);
    ImageIOAnalyze(float px, float py, float pz, float pt);
    ImageIOAnalyze(float px, float py, float pz, float pt, float pr);
    ImageIOAnalyze(float px, float py, float pz, float pt, float pr, float ps);
    ImageIOAnalyze(float px, float py, float pz, float pt, float pr, float ps, float pp);
    ImageIOAnalyze(float px, float py, float pz, float pt, float pr, float ps, float pp, float pq);

    virtual ~ImageIOAnalyze() {}

    virtual void export_array(const hoNDArray<short>& a, const std::string& filename) { this->export_array_impl(a, filename); }

    virtual void export_array(const hoNDArray<unsigned short>& a, const std::string& filename) { this->export_array_impl(a, filename); }
    virtual void export_array(const hoNDArray<int>& a, const std::string& filename) { this->export_array_impl(a, filename); }
    virtual void export_array(const hoNDArray<unsigned int>& a, const std::string& filename) { this->export_array_impl(a, filename); }
    virtual void export_array(const hoNDArray<size_t>& a, const std::string& filename) { this->export_array_impl(a, filename); }
    virtual void export_array(const hoNDArray<float>& a, const std::string& filename) { this->export_array_impl(a, filename); }
    virtual void export_array(const hoNDArray<double>& a, const std::string& filename) { this->export_array_impl(a, filename); }
    virtual void export_array(const hoNDArray< std::complex<float> >& a, const std::string& filename) { this->export_array_impl(a, filename); }
    virtual void export_array(const hoNDArray< std::complex<double> >& a, const std::string& filename) { this->export_array_impl(a, filename); }

    virtual void import_array(hoNDArray<short>& a, const std::string& filename) { this->import_array_impl(a, filename); }
    virtual void import_array(hoNDArray<unsigned short>& a, const std::string& filename) { this->import_array_impl(a, filename); }
    virtual void import_array(hoNDArray<int>& a, const std::string& filename) { this->import_array_impl(a, filename); }
    virtual void import_array(hoNDArray<unsigned int>& a, const std::string& filename) { this->import_array_impl(a, filename); }
    virtual void import_array(hoNDArray<float>& a, const std::string& filename) { this->import_array_impl(a, filename); }
    virtual void import_array(hoNDArray<double>& a, const std::string& filename) { this->import_array_impl(a, filename); }
    virtual void import_array(hoNDArray< std::complex<float> >& a, const std::string& filename) { this->import_array_impl(a, filename); }
    virtual void import_array(hoNDArray< std::complex<double> >& a, const std::string& filename) { this->import_array_impl(a, filename); }

    template <typename T> 
    void export_array_impl(const hoNDArray<T>& a, const std::string& filename)
    {
        try
        {
            // Check if folder exists
            std::filesystem::path folder_path(filename);
            folder_path.remove_filename();
            if ( !std::filesystem::is_directory(folder_path) )
            {
                GWARN_STREAM("Failed to write " << filename << " because parent folder does not exist");
                return;
            }

            HeaderType header;
            GADGET_CHECK_THROW(this->array_to_header(a, header));
            GADGET_CHECK_THROW(this->write_header(filename, header));

            std::string filenameData = filename;
            filenameData.append(".img");
            this->write_data(filenameData, a.begin(), a.get_number_of_bytes());
        }
        catch(...)
        {
            GADGET_THROW("Errors in ImageIOAnalyze::export_array_impl(const hoNDArray<T>& a, const std::string& filename) ... ");
        }
    }

    template <typename T> 
    void import_array_impl(hoNDArray<T>& a, const std::string& filename)
    {
        try
        {
            HeaderType header;
            GADGET_CHECK_THROW(this->read_header(filename, header));
            GADGET_CHECK_THROW(this->header_to_array(a, header));

            std::string filenameData = filename;
            filenameData.append(".img");
            this->read_data(filenameData, a.begin(), a.get_number_of_bytes());
        }
        catch(...)
        {
            GADGET_THROW("Errors in ImageIOAnalyze::import_array_impl(const hoNDArray<T>& a, const std::string& filename) ... ");
        }
    }

    template <typename T, unsigned int D> 
    void export_image(const hoNDImage<T,D>& a, const std::string& filename)
    {
        try
        {
            HeaderType header;
            GADGET_CHECK_THROW(this->image_to_header(a, header));
            GADGET_CHECK_THROW(this->write_header(filename, header));

            std::string filenameData = filename;
            filenameData.append(".img");
            this->write_data(filenameData, a.begin(), a.get_number_of_bytes());
        }
        catch(...)
        {
            GADGET_THROW("Errors in export_image(const hoNDImage<T,D>& a, const std::string& filename) ... ");
        }
    }

    template <typename T, unsigned int D> 
    void import_image(hoNDImage<T,D>& a, const std::string& filename)
    {
        try
        {
            HeaderType header;
            GADGET_CHECK_THROW(this->read_header(filename, header));
            GADGET_CHECK_THROW(this->header_to_image(a, header));

            std::string filenameData = filename;
            filenameData.append(".img");
            this->read_data(filenameData, a.begin(), a.get_number_of_bytes());
        }
        catch(...)
        {
            GADGET_THROW("Errors in import_image(const hoNDImage<T,D>& a, const std::string& filename) ... ");
        }
    }

/// image functions

    template <typename T, unsigned int D> 
    void export_image_complex(const hoNDImage<T,D>& a, const std::string& filename)
    {
        try
        {
            typedef typename Gadgetron::realType<T>::Type value_type;

            long long num = (long long)a.get_number_of_elements();

            long long n;

            hoNDImage<value_type, D> rpart, ipart, mag, phs;
            rpart.create( a.get_dimensions() );
            ipart.create( a.get_dimensions() );
            mag.create( a.get_dimensions() );
            phs.create( a.get_dimensions() );

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
            GADGET_CATCH_THROW(export_image(rpart, filenameReal));

            std::string filenameImag = filename;
            filenameImag.append("_IMAG");
            GADGET_CATCH_THROW(export_image(ipart, filenameImag));

            std::string filenameMag = filename;
            filenameMag.append("_MAG");
            GADGET_CATCH_THROW(export_image(mag, filenameMag));

            std::string filenamePhase = filename;
            filenamePhase.append("_PHASE");
            GADGET_CATCH_THROW(export_image(phs, filenamePhase));
        }
        catch(...)
        {
            GADGET_THROW("Errors in export_image_complex(const hoNDImage<T,D>& a, const std::string& filename) ... ");
        }
    }

    template <typename T, unsigned int D> 
    void import_image_complex(hoNDImage<T,D>& a, const std::string& filename)
    {
        try
        {
            typedef typename T::value_type value_type;
            hoNDImage<value_type, D> real, imag;

            std::string filenameReal = filename;
            filenameReal.append("_REAL");
            GADGET_CATCH_THROW(import_image(real, filenameReal));

            std::string filenameImag = filename;
            filenameImag.append("_IMAG");
            GADGET_CATCH_THROW(import_image(imag, filenameImag));

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
            GADGET_THROW("Errors in import_image_complex(const hoNDImage<T,D>& a, const std::string& filename) ... ");
        }
    }

    template <typename T, unsigned int D> 
    void import_image_complex(hoNDImage<T,D>& a, const std::string& filename_real, const std::string& filename_imag)
    {
        try
        {
            typedef typename realType<T>::Type value_type;
            hoNDImage<value_type, D> real, imag;

            GADGET_CATCH_THROW(import_image(real, filename_real));
            GADGET_CATCH_THROW(import_image(imag, filename_imag));

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
            GADGET_THROW("Errors in import_image_complex(hoNDImage<T,D>& a, const std::string& filename_real, const std::string& filename_imag) ... ");
        }
    }

    template <typename T> 
    void export_2d_image(const hoNDImage<T,2>& a, const std::string& filename)
    {
        export_image(a, filename);
    }

    template <typename T> 
    void import_2d_image(hoNDImage<T,2>& a, const std::string& filename)
    {
        import_image(a, filename);
    }

    template <typename T> 
    void export_2d_image_complex(const hoNDImage<T,2>& a, const std::string& filename)
    {
        export_image_complex(a, filename);
    }

    template <typename T> 
    void import_2d_image_complex(hoNDImage<T,2>& a, const std::string& filename)
    {
        import_image_complex(a, filename);
    }

    template <typename T> 
    void export_3d_image(const hoNDImage<T,3>& a, const std::string& filename)
    {
        export_image(a, filename);
    }

    template <typename T> 
    void import_3d_image(hoNDImage<T,3>& a, const std::string& filename)
    {
        import_image(a, filename);
    }

    template <typename T> 
    void export_3d_image_complex(const hoNDImage<T,3>& a, const std::string& filename)
    {
        export_image_complex(a, filename);
    }

    template <typename T> 
    void import_3d_image_complex(hoNDImage<T,3>& a, const std::string& filename)
    {
        import_image_complex(a, filename);
    }

    template <typename T> 
    void export_4d_image(const hoNDImage<T,4>& a, const std::string& filename)
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
                GADGET_CATCH_THROW(export_3d_image(a3D, ostr.str()));
            }
        }
        catch(...)
        {
            GADGET_THROW("Errors in export_4d_image(const hoNDImage<T>& a, const std::string& filename) ... ");
        }
    }

    template <typename T> 
    void export_4d_image_complex(const hoNDImage<T,4>& a, const std::string& filename)
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
                GADGET_CATCH_THROW(export_3d_image_complex(a3D, ostr.str()));
            }
        }
        catch(...)
        {
            GADGET_THROW("Errors in export_4d_image_complex(const hoNDImage<T>& a, const std::string& filename) ... ");
        }
    }

protected:

    template <typename T> bool array_to_header(const hoNDArray<T>& a, HeaderType& header);
    template <typename T> bool header_to_array(hoNDArray<T>& a, const HeaderType& header);

    template <typename T, unsigned int D> bool image_to_header(const hoNDImage<T, D>& a, HeaderType& header);
    template <typename T, unsigned int D> bool header_to_image(hoNDImage<T, D>& a, const HeaderType& header);

    // read/write the analyze header
    bool read_header(const std::string& filename, HeaderType& header);
    bool write_header(const std::string& filename, const HeaderType& header);
};

template <typename T> 
bool ImageIOAnalyze::array_to_header(const hoNDArray<T>& a, HeaderType& header)
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
        GERROR_STREAM("Errors in ImageIOAnalyze::array2Analyze(const hoNDArray<T>& a, dsr& header) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool ImageIOAnalyze::header_to_array(hoNDArray<T>& a, const HeaderType& header)
{
    try
    {
        std::string rttiID = std::string(typeid(T).name());
        GADGET_CHECK_THROW(rttiID==getRTTIFromDataType( (ImageIODataType)header.dime.datatype));

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

        a.create(dim);
    }
    catch(...)
    {
        GERROR_STREAM("Errors in ImageIOAnalyze::analyze2Array(hoNDArray<T>& a, const dsr& header) ... ");
        return false;
    }

    return true;
}

template <typename T, unsigned int D> 
bool ImageIOAnalyze::image_to_header(const hoNDImage<T,D>& a, HeaderType& header)
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
        GERROR_STREAM("Errors in ImageIOAnalyze::image2Analyze(const hoNDImage<T>& a, dsr& header) ... ");
        return false;
    }

    return true;
}

template <typename T, unsigned int D> 
bool ImageIOAnalyze::header_to_image(hoNDImage<T,D>& a, const HeaderType& header)
{
    try
    {
        std::string rttiID = std::string(typeid(T).name());
        GADGET_CHECK_THROW(rttiID==getRTTIFromDataType( (ImageIODataType)header.dime.datatype));

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
        GERROR_STREAM("Errors in ImageIOAnalyze::analyze2Image(hoNDImage<T,D>& a, const dsr& header) ... ");
        return false;
    }

    return true;
}

}
