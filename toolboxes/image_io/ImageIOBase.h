/** \file       ImageIOBase.h
    \brief      Define the base IO funcatinality for image_analyze_io toolbox
    \author     Hui Xue
*/

#pragma once

#include <iostream>
#include <typeinfo>

#include "NDArray.h"
#include "complext.h"

#include "hoNDArray.h"
#include "hoNDImage.h"

#include "hoNDArray_fileio.h"

namespace Gadgetron {

struct rgb_type { unsigned char r,g,b; };
struct rgba_type { unsigned char r,g,b,a; };

class ImageIOWorker
{
public:

    ImageIOWorker(const std::string& ioTag, bool readFlag=true);
    virtual ~ImageIOWorker();

    // open the file stream
    // readFlag: true, read mode; false, write mode
    virtual bool open();

    // close the file stream
    virtual bool close();

    // the current file offset
    long tell();

    // set the file offset
    bool seek(long long offset);

    // reset the file to the beginning
    bool reset() { return (this->seek(0)); }

    // check the status of i/o operations
    bool IOinError();

    // read/write
    // len: number of bytes in data
    bool read(char* data, long long len);
    bool write(const char* data, long long len);

protected:

    std::string ioTag_;
    std::fstream fid_;
    bool readFlag_;
};

#ifdef DT_UNKNOWN
    #undef DT_UNKNOWN
#endif // DT_UNKNOWN

enum ImageIODataType
{
    DT_ANA_UNKNOWN             =0,
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

template <typename HeaderType>
class ImageIOBase
{
public:

    typedef HeaderType THeaderType;

    ImageIOBase()
    {
        pixelSize_.resize(10, 1.0);
    }

    ImageIOBase(float px, float py)
    {
        pixelSize_.resize(2);
        pixelSize_[0] = px;
        pixelSize_[1] = py;
    }

    ImageIOBase(float px, float py, float pz)
    {
        pixelSize_.resize(3);
        pixelSize_[0] = px;
        pixelSize_[1] = py;
        pixelSize_[2] = pz;
    }

    ImageIOBase(float px, float py, float pz, float pt)
    {
        pixelSize_.resize(4);
        pixelSize_[0] = px;
        pixelSize_[1] = py;
        pixelSize_[2] = pz;
        pixelSize_[3] = pt;
    }

    ImageIOBase(float px, float py, float pz, float pt, float pr)
    {
        pixelSize_.resize(5);
        pixelSize_[0] = px;
        pixelSize_[1] = py;
        pixelSize_[2] = pz;
        pixelSize_[3] = pt;
        pixelSize_[4] = pr;
    }

    ImageIOBase(float px, float py, float pz, float pt, float pr, float ps)
    {
        pixelSize_.resize(6);
        pixelSize_[0] = px;
        pixelSize_[1] = py;
        pixelSize_[2] = pz;
        pixelSize_[3] = pt;
        pixelSize_[4] = pr;
        pixelSize_[5] = ps;
    }

    ImageIOBase(float px, float py, float pz, float pt, float pr, float ps, float pp)
    {
        pixelSize_.resize(7);
        pixelSize_[0] = px;
        pixelSize_[1] = py;
        pixelSize_[2] = pz;
        pixelSize_[3] = pt;
        pixelSize_[4] = pr;
        pixelSize_[5] = ps;
        pixelSize_[6] = pp;
    }

    ImageIOBase(float px, float py, float pz, float pt, float pr, float ps, float pp, float pq)
    {
        pixelSize_.resize(8);
        pixelSize_[0] = px;
        pixelSize_[1] = py;
        pixelSize_[2] = pz;
        pixelSize_[3] = pt;
        pixelSize_[4] = pr;
        pixelSize_[5] = ps;
        pixelSize_[6] = pp;
        pixelSize_[7] = pq;
    }

    void set_pixel_size(float px, float py, float pz=1.0f, float pt=1.0f, float pr=1.0f, float ps=1.0f, float pp=1.0f, float pq=1.0f)
    {
        pixelSize_.resize(8);
        pixelSize_[0] = px;
        pixelSize_[1] = py;
        pixelSize_[2] = pz;
        pixelSize_[3] = pt;
        pixelSize_[4] = pr;
        pixelSize_[5] = ps;
        pixelSize_[6] = pp;
        pixelSize_[7] = pq;
    }

    void set_pixel_size(double px, double py, double pz=1.0, double pt=1.0, double pr=1.0, double ps=1.0, double pp=1.0, double pq=1.0)
    {
        pixelSize_.resize(8);
        pixelSize_[0] = (float)px;
        pixelSize_[1] = (float)py;
        pixelSize_[2] = (float)pz;
        pixelSize_[3] = (float)pt;
        pixelSize_[4] = (float)pr;
        pixelSize_[5] = (float)ps;
        pixelSize_[6] = (float)pp;
        pixelSize_[7] = (float)pq;
    }

    virtual ~ImageIOBase()
    {
    }

    /// export/input for 2D/3D/4D array
    /// filename should be given without extension

    virtual void export_array(const hoNDArray<short>& a, const std::string& filename) = 0;
    virtual void export_array(const hoNDArray<unsigned short>& a, const std::string& filename) = 0;
    virtual void export_array(const hoNDArray<int>& a, const std::string& filename) = 0;
    virtual void export_array(const hoNDArray<unsigned int>& a, const std::string& filename) = 0;
    virtual void export_array(const hoNDArray<float>& a, const std::string& filename) = 0;
    virtual void export_array(const hoNDArray<double>& a, const std::string& filename) = 0;
    virtual void export_array(const hoNDArray< std::complex<float> >& a, const std::string& filename) = 0;
    virtual void export_array(const hoNDArray< std::complex<double> >& a, const std::string& filename) = 0;

    virtual void import_array(hoNDArray<short>& a, const std::string& filename) = 0;
    virtual void import_array(hoNDArray<unsigned short>& a, const std::string& filename) = 0;
    virtual void import_array(hoNDArray<int>& a, const std::string& filename) = 0;
    virtual void import_array(hoNDArray<unsigned int>& a, const std::string& filename) = 0;
    virtual void import_array(hoNDArray<float>& a, const std::string& filename) = 0;
    virtual void import_array(hoNDArray<double>& a, const std::string& filename) = 0;
    virtual void import_array(hoNDArray< std::complex<float> >& a, const std::string& filename) = 0;
    virtual void import_array(hoNDArray< std::complex<double> >& a, const std::string& filename) = 0;

    template <typename T>
    void export_array_complex_real_imag(const hoNDArray<T>& a, const std::string& filename)
    {
        try
        {
            typedef typename Gadgetron::realType<T>::Type value_type;

            hoNDArray<value_type> buf(a.dimensions());

            long long num = (long long)a.get_number_of_elements();

            long long n;
            #pragma omp parallel for default(none) private(n) shared(num, a, buf)
            for ( n=0; n<num; n++ )
            {
                buf(n) = a(n).real();
            }

            std::string filenameReal = filename;
            filenameReal.append("_REAL");
            export_array(buf, filenameReal);

            #pragma omp parallel for default(none) private(n) shared(num, a, buf)
            for ( n=0; n<num; n++ )
            {
                buf(n) = a(n).imag();
            }

            std::string filenameImag = filename;
            filenameImag.append("_IMAG");
            export_array(buf, filenameImag);
        }
        catch(...)
        {
            GADGET_THROW("Errors in export_array_complex_real_imag(const hoNDArray<T>& a, const std::string& filename) ... ");
        }
    }

    template <typename T>
    void export_array_complex(const hoNDArray<T>& a, const std::string& filename)
    {
        try
        {
            typedef typename Gadgetron::realType<T>::Type value_type;

            hoNDArray<value_type> rpart, ipart, mag, phs;
            rpart.create(a.dimensions());
            ipart.create(a.dimensions());
            mag.create(a.dimensions());
            phs.create(a.dimensions());

            long long num = (long long)a.get_number_of_elements();

            long long n;

            #pragma omp parallel for default(none) private(n) shared(num, a, rpart, ipart, mag, phs)
            for ( n=0; n<num; n++ )
            {
                rpart(n) = a(n).real();
                ipart(n) = a(n).imag();
                mag(n) = std::abs( a(n) );
                phs(n) = std::arg( a(n) );
            }

            std::string filenameReal = filename;
            filenameReal.append("_REAL");
            GADGET_CATCH_THROW(export_array(rpart, filenameReal));

            std::string filenameImag = filename;
            filenameImag.append("_IMAG");
            GADGET_CATCH_THROW(export_array(ipart, filenameImag));

            std::string filenameMag = filename;
            filenameMag.append("_MAG");
            GADGET_CATCH_THROW(export_array(mag, filenameMag));

            std::string filenamePhase = filename;
            filenamePhase.append("_PHASE");
            GADGET_CATCH_THROW(export_array(phs, filenamePhase));
        }
        catch(...)
        {
            GADGET_THROW("Errors in export_array_complex(const hoNDArray<T>& a, const std::string& filename) ... ");
        }
    }

    template <typename T>
    void import_array_complex(hoNDArray<T>& a, const std::string& filename)
    {
        try
        {
            typedef typename T::value_type value_type;
            hoNDArray<value_type> real, imag;

            std::string filenameReal = filename;
            filenameReal.append("_REAL");
            GADGET_CATCH_THROW(import_array(real, filenameReal));

            std::string filenameImag = filename;
            filenameImag.append("_IMAG");
            GADGET_CATCH_THROW(import_array(imag, filenameImag));

            a.create(real.dimensions());
            long long num = (long long)real.get_number_of_elements();

            long long n;
            #pragma omp parallel for private(n) shared(num, a, real, imag)
            for ( n=0; n<num; n++ )
            {
                a(n) = T(real(n), imag(n));
            }
        }
        catch(...)
        {
            GADGET_THROW("Errors in import_array_complex(const hoNDArray<T>& a, const std::string& filename) ... ");
        }
    }

    template <typename T>
    void import_array_complex(hoNDArray<T>& a, const std::string& filename_real, const std::string& filename_imag)
    {
        try
        {
            typedef typename realType<T>::Type value_type;
            hoNDArray<value_type> real, imag;

            GADGET_CATCH_THROW(import_array(real, filename_real));
            GADGET_CATCH_THROW(import_array(imag, filename_imag));

            a.create(real.dimensions());
            long long num = (long long)real.get_number_of_elements();

            long long n;
            #pragma omp parallel for private(n) shared(num, a, real, imag)
            for ( n=0; n<num; n++ )
            {
                a(n) = T(real(n), imag(n));
            }
        }
        catch(...)
        {
            GADGET_THROW("Errors in import_array_complex(hoNDArray<T>& a, const std::string& filename_real, const std::string& filename_imag) ... ");
        }
    }

    template <typename T>
    void export_2d_array(const hoNDArray<T>& a, const std::string& filename)
    {
        export_array(a, filename);
    }

    template <typename T>
    void import_2d_array(hoNDArray<T>& a, const std::string& filename)
    {
        import_array(a, filename);
    }

    template <typename T>
    void export_2d_array_complex(const hoNDArray<T>& a, const std::string& filename)
    {
        export_array_complex(a, filename);
    }

    template <typename T>
    void import_2d_array_complex(hoNDArray<T>& a, const std::string& filename)
    {
        import_array_complex(a, filename);
    }

    template <typename T>
    void export_3d_array(const hoNDArray<T>& a, const std::string& filename)
    {
        export_array(a, filename);
    }

    template <typename T>
    void import_3d_array(hoNDArray<T>& a, const std::string& filename)
    {
        import_array(a, filename);
    }

    template <typename T>
    void export_3d_array_complex(const hoNDArray<T>& a, const std::string& filename)
    {
        export_array_complex(a, filename);
    }

    template <typename T>
    void import_3d_array_complex(hoNDArray<T>& a, const std::string& filename)
    {
        import_array_complex(a, filename);
    }

    template <typename T>
    void export_4d_array(const hoNDArray<T>& a, const std::string& filename)
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

                hoNDArray<T> a3D(dim, const_cast<T*>(a.begin()+ii*RO*E1*CHA), false);

                std::ostringstream ostr;
                ostr << filename << "_" << ii << std::ends;
                GADGET_CATCH_THROW(export3DArray(a3D, ostr.str()));
            }
        }
        catch(...)
        {
            GADGET_THROW("Errors in export_4d_array(const hoNDArray<T>& a, const std::string& filename) ... ");
        }
    }

    template <typename T>
    void export_4d_array_complex(const hoNDArray<T>& a, const std::string& filename)
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

                hoNDArray<T> a3D(dim, const_cast<T*>(a.begin()+ii*RO*E1*CHA), false);

                std::ostringstream ostr;
                ostr << filename << "_" << ii << std::ends;
                GADGET_CATCH_THROW(export3DArrayComplex(a3D, ostr.str()));
            }
        }
        catch(...)
        {
            GADGET_THROW("Errors in export_4d_array_complex(const hoNDArray<T>& a, const std::string& filename) ... ");
        }
    }

    static void readFromFile(const std::string& filename, char*& data, long long& length)
    {
        try
        {
            if (data!=NULL) delete [] data;

            ImageIOWorker ioworker_(filename, true);

            GADGET_CHECK_THROW(ioworker_.open());

            // read the total length
            long long totalLen;
            GADGET_CHECK_THROW(ioworker_.read(reinterpret_cast<char*>(&totalLen), sizeof(long long)));

            length = totalLen - sizeof(long long);

            data = new char[length];
            GADGET_CHECK_THROW(data!=NULL);

            GADGET_CHECK_THROW(ioworker_.read(data, length));

            GADGET_CHECK_THROW(ioworker_.close());
        }
        catch (...)
        {
            GADGET_THROW("Errors in ImageIOBase::readFromFile(const std::string& filename, char*& data, long long& length) ... ");
        }
    }

    static void writeToFile(const std::string& filename, char* data, long long length)
    {
        try
        {
            if ( length == 0 ) return;

            GADGET_CHECK_THROW(data!=NULL);

            ImageIOWorker ioworker_(filename, false);

            GADGET_CHECK_THROW(ioworker_.open());

            // write the total lengh
            const long long totalLen = length+sizeof(long long);
            GADGET_CHECK_THROW(ioworker_.write(reinterpret_cast<const char*>(&totalLen), sizeof(long long)));

            // write the data
            GADGET_CHECK_THROW(ioworker_.write(data, length));

            // close the file
            GADGET_CHECK_THROW(ioworker_.close());
        }
        catch (...)
        {
            GADGET_THROW("Errors in ImageIOBase::writeToFile(const std::string& filename, char* data, long long length) ... ");
        }
    }

protected:

    std::vector<float> pixelSize_;

    // get the run-time type ID from analyze data type or vice versa
    std::string getRTTIFromDataType(ImageIODataType aDT)
    {
        std::string rttiID;

        switch (aDT)
        {
        case DT_INT8 :
            rttiID = typeid(char).name();
            break;

        case DT_UNSIGNED_CHAR :
            rttiID = typeid(unsigned char).name();
            break;

        case DT_SIGNED_SHORT :
            rttiID = typeid(short).name();
            break;

        case DT_UNSIGNED_SHORT :
        case DT_UINT16 :
            rttiID = typeid(unsigned short).name();
            break;

        case DT_SIGNED_INT :
            rttiID = typeid(int).name();
            break;

        case DT_UINT32 :
            rttiID = typeid(unsigned int).name();
            break;

        case DT_INT64 :
            rttiID = typeid(long long).name();
            break;

        case DT_UINT64 :
            rttiID = typeid(unsigned long long).name();
            break;

        case DT_FLOAT :
            rttiID = typeid(float).name();
            break;

        case DT_DOUBLE :
            rttiID = typeid(double).name();
            break;

        case DT_FLOAT128 :
            rttiID = typeid(long double).name();
            break;

        case DT_COMPLEX :
            rttiID = typeid( std::complex<float> ).name();
            break;

        case DT_COMPLEX128 :
            rttiID = typeid( std::complex<double> ).name();
            break;

        case DT_COMPLEX256 :
            rttiID = typeid( std::complex<long double> ).name();
            break;

        case DT_RGB :
            rttiID = typeid( Gadgetron::rgb_type ).name();
            break;

        case DT_RGBA32 :
            rttiID = typeid( Gadgetron::rgba_type ).name();
            break;

        default:
            rttiID = "UNKOWN TYPE";
        }

        return rttiID;
    }

    ImageIODataType getDataTypeFromRTTI(const std::string& name)
    {
        ImageIODataType analyzeDT = DT_ANA_UNKNOWN;

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

        if ( name == typeid( std::complex<float> ).name() )
        {
            analyzeDT = DT_COMPLEX;
        }

        if ( name == typeid( std::complex<double> ).name() )
        {
            analyzeDT = DT_COMPLEX128;
        }

        if ( name == typeid(std::complex<long double>).name() )
        {
            analyzeDT = DT_COMPLEX256;
        }

        if ( name == typeid(Gadgetron::rgb_type).name() )
        {
            analyzeDT = DT_RGB;
        }

        if ( name == typeid(Gadgetron::rgba_type).name() )
        {
            analyzeDT = DT_RGBA32;
        }

        return analyzeDT;
    }

    template <typename T>
    void read_data(const std::string& filename, T* data, long long len)
    {
        try
        {
            ImageIOWorker ioworker(filename, true);

            GADGET_CHECK_THROW(ioworker.open());
            GADGET_CHECK_THROW(ioworker.read(reinterpret_cast<char*>(data), len));
            GADGET_CHECK_THROW(ioworker.close());
        }
        catch(...)
        {
            GADGET_THROW("Errors in read_data(const std::string& filename, T* data, long long len) ... ");
        }
    }

    template <typename T>
    void write_data(const std::string& filename, const T* data, long long len)
    {
        try
        {
            ImageIOWorker ioworker(filename, false);

            GADGET_CHECK_THROW(ioworker.open());
            GADGET_CHECK_THROW(ioworker.write(reinterpret_cast<const char*>(data), len));
            GADGET_CHECK_THROW(ioworker.close());
        }
        catch(...)
        {
            GADGET_THROW("Errors in write_data(const std::string& filename, const T* data, long long len) ... ");
        }
    }
};

}
