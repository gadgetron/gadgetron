/** \file       gtPlusIONifti.h
    \brief      Implement the suppor for the Nifti medical image format
    \author     Hui Xue

    The ISMRMRD dimensions are mapped to Nifti format.

    Ref to:
    http://eeg.sourceforge.net/ANALYZE75.pdf
    http://ismrmrd.sourceforge.net/
*/

#pragma once

#include <gtPlusIOBase.h>

// ------------------------------------------------------------------------------------------------ //

// the file input/output utility functions for the Nfiti format

// the following Nifti data structured is defined as this online document
// http://nifti.nimh.nih.gov/pub/dist/src/niftilib/nifti1.h

// the official definition of Nfiti file format
struct nifti_1_header 
{ /* NIFTI-1 usage         */  /* ANALYZE 7.5 field(s) */
    /*************************/  /************************/

    /*--- was header_key substruct ---*/
    int   sizeof_hdr;    /*!< MUST be 348           */  /* int sizeof_hdr;      */
    char  data_type[10]; /*!< ++UNUSED++            */  /* char data_type[10];  */
    char  db_name[18];   /*!< ++UNUSED++            */  /* char db_name[18];    */
    int   extents;       /*!< ++UNUSED++            */  /* int extents;         */
    short session_error; /*!< ++UNUSED++            */  /* short session_error; */
    char  regular;       /*!< ++UNUSED++            */  /* char regular;        */
    char  dim_info;      /*!< MRI slice ordering.   */  /* char hkey_un0;       */

    /*--- was image_dimension substruct ---*/
    short dim[8];        /*!< Data array dimensions.*/  /* short dim[8];        */
    float intent_p1 ;    /*!< 1st intent parameter. */  /* short unused8;       */
    /* short unused9;       */
    float intent_p2 ;    /*!< 2nd intent parameter. */  /* short unused10;      */
    /* short unused11;      */
    float intent_p3 ;    /*!< 3rd intent parameter. */  /* short unused12;      */
    /* short unused13;      */
    short intent_code ;  /*!< NIFTI_INTENT_* code.  */  /* short unused14;      */
    short datatype;      /*!< Defines data type!    */  /* short datatype;      */
    short bitpix;        /*!< Number bits/voxel.    */  /* short bitpix;        */
    short slice_start;   /*!< First slice index.    */  /* short dim_un0;       */
    float pixdim[8];     /*!< Grid spacings.        */  /* float pixdim[8];     */
    float vox_offset;    /*!< Offset into .nii file */  /* float vox_offset;    */
    float scl_slope ;    /*!< Data scaling: slope.  */  /* float funused1;      */
    float scl_inter ;    /*!< Data scaling: offset. */  /* float funused2;      */
    short slice_end;     /*!< Last slice index.     */  /* float funused3;      */
    char  slice_code ;   /*!< Slice timing order.   */
    char  xyzt_units ;   /*!< Units of pixdim[1..4] */
    float cal_max;       /*!< Max display intensity */  /* float cal_max;       */
    float cal_min;       /*!< Min display intensity */  /* float cal_min;       */
    float slice_duration;/*!< Time for 1 slice.     */  /* float compressed;    */
    float toffset;       /*!< Time axis shift.      */  /* float verified;      */
    int   glmax;         /*!< ++UNUSED++            */  /* int glmax;           */
    int   glmin;         /*!< ++UNUSED++            */  /* int glmin;           */

    /*--- was data_history substruct ---*/
    char  descrip[80];   /*!< any text you like.    */  /* char descrip[80];    */
    char  aux_file[24];  /*!< auxiliary filename.   */  /* char aux_file[24];   */

    short qform_code ;   /*!< NIFTI_XFORM_* code.   */  /*-- all ANALYZE 7.5 ---*/
    short sform_code ;   /*!< NIFTI_XFORM_* code.   */  /*   fields below here  */
    /*   are replaced       */
    float quatern_b ;    /*!< Quaternion b param.   */
    float quatern_c ;    /*!< Quaternion c param.   */
    float quatern_d ;    /*!< Quaternion d param.   */
    float qoffset_x ;    /*!< Quaternion x shift.   */
    float qoffset_y ;    /*!< Quaternion y shift.   */
    float qoffset_z ;    /*!< Quaternion z shift.   */

    float srow_x[4] ;    /*!< 1st row affine transform.   */
    float srow_y[4] ;    /*!< 2nd row affine transform.   */
    float srow_z[4] ;    /*!< 3rd row affine transform.   */

    char intent_name[16];/*!< 'name' or meaning of data.  */

    char magic[4] ;      /*!< MUST be "ni1\0" or "n+1\0". */

    char extension[4];
} ;                   /**** 348+4 bytes total ****/

typedef struct nifti_1_header nifti_1_header ;

/*! \defgroup NIFTI1_SLICE_ORDER
    \brief nifti1 slice order codes, describing the acquisition order
           of the slices
    @{
 */
#define NIFTI_SLICE_UNKNOWN   0
#define NIFTI_SLICE_SEQ_INC   1
#define NIFTI_SLICE_SEQ_DEC   2
#define NIFTI_SLICE_ALT_INC   3
#define NIFTI_SLICE_ALT_DEC   4
#define NIFTI_SLICE_ALT_INC2  5  /* 05 May 2005: RWCox */
#define NIFTI_SLICE_ALT_DEC2  6  /* 05 May 2005: RWCox */
/* @} */

#define NIFTI_UNITS_UNKNOWN 0 // NIFTI code for unspecified units. /
                              // Space codes are multiples of 1. /
#define NIFTI_UNITS_METER   1 // NIFTI code for meters. /
#define NIFTI_UNITS_MM      2 // NIFTI code for millimeters. /
#define NIFTI_UNITS_MICRON  3 // NIFTI code for micrometers. /


                               // Time codes are multiples of 8. /
#define NIFTI_UNITS_SEC     8 // NIFTI code for seconds. /
#define NIFTI_UNITS_MSEC   16 // NIFTI code for milliseconds. /
#define NIFTI_UNITS_USEC   24 // NIFTI code for microseconds. /


                               // These units are for spectral data: /
#define NIFTI_UNITS_HZ     32 // NIFTI code for Hertz. /
#define NIFTI_UNITS_PPM    40 // NIFTI code for ppm. /
#define NIFTI_UNITS_RADS   48 // NIFTI code for radians per second. */

#define XYZT_TO_SPACE(xyzt)       ( (xyzt) & 0x07 )
#define XYZT_TO_TIME(xyzt)        ( (xyzt) & 0x38 )
#define SPACE_TIME_TO_XYZT(ss,tt) (  (((char)(ss)) & 0x07)   \
                                   | (((char)(tt)) & 0x38) )

/* [qs]form_code value:  /      / x,y,z coordinate system refers to:    /
/-----------------------/      /---------------------------------------*/

#define NIFTI_XFORM_UNKNOWN      0 /* Arbitrary coordinates (Method 1). */
#define NIFTI_XFORM_SCANNER_ANAT 1 /* Scanner-based anatomical coordinates */
#define NIFTI_XFORM_ALIGNED_ANAT 2 /* Coordinates aligned to another file's, or to anatomical "truth".           */
#define NIFTI_XFORM_TALAIRACH    3 /* Coordinates aligned to Talairach-Tournoux Atlas; (0,0,0)=AC, etc.*/
#define NIFTI_XFORM_MNI_152      4 /* MNI 152 normalized coordinates. */

/*! \struct nifti1_extender
    \brief This structure represents a 4-byte string that should follow the
           binary nifti_1_header data in a NIFTI-1 header file.  If the char
           values are {1,0,0,0}, the file is expected to contain extensions,
           values of {0,0,0,0} imply the file does not contain extensions.
           Other sequences of values are not currently defined.
 */
struct nifti1_extender { char extension[4] ; } ;
typedef struct nifti1_extender nifti1_extender ;

/*! \struct nifti1_extension
    \brief Data structure defining the fields of a header extension.
 */
struct nifti1_extension {
   int    esize ; /*!< size of extension, in bytes (must be multiple of 16) */
   int    ecode ; /*!< extension code, one of the NIFTI_ECODE_ values       */
   char * edata ; /*!< raw data, with no byte swapping (length is esize-8)  */

   nifti1_extension()
   {
        esize = sizeof(int) + sizeof(int);
        ecode = 0;
        edata = NULL;
   }

   ~nifti1_extension()
   {
        if ( !edata ) delete [] edata;
   }
} ;
typedef struct nifti1_extension nifti1_extension ;

/*.................*/
/*! Given a nifti_1_header struct, check if it has a good magic number.
    Returns NIFTI version number (1..9) if magic is good, 0 if it is not. */

#define NIFTI_VERSION(h)                               \
 ( ( (h).magic[0]=='n' && (h).magic[3]=='\0'    &&     \
     ( (h).magic[1]=='i' || (h).magic[1]=='+' ) &&     \
     ( (h).magic[2]>='1' && (h).magic[2]<='9' )   )    \
 ? (h).magic[2]-'0' : 0 )

/*.................*/
/*! Check if a nifti_1_header struct says if the data is stored in the
    same file or in a separate file.  Returns 1 if the data is in the same
    file as the header, 0 if it is not.                                   */

#define NIFTI_ONEFILE(h) ( (h).magic[1] == '+' )

/*.................*/
/*! Check if a nifti_1_header struct needs to be byte swapped.
    Returns 1 if it needs to be swapped, 0 if it does not.     */

#define NIFTI_NEEDS_SWAP(h) ( (h).dim[0] < 0 || (h).dim[0] > 7 )

// ------------------------------------------------------------------------------------------------ //
// to suppor the ISMRMRD format
// [Ro E1 Cha Slice E2 Con Phase Rep Set Seg]

namespace Gadgetron { namespace gtPlus {

    class EXPORTGTPLUS gtPlusIONifti : public gtPlusIOBase<nifti_1_header>
    {
    public:

        typedef gtPlusIOBase<nifti_1_header> BaseClass;
        typedef BaseClass::THeaderType HeaderType;

        gtPlusIONifti();
        gtPlusIONifti(float px, float py);
        gtPlusIONifti(float px, float py, float pz);
        gtPlusIONifti(float px, float py, float pz, float pt);
        gtPlusIONifti(float px, float py, float pz, float pt, float pr);
        gtPlusIONifti(float px, float py, float pz, float pt, float pr, float ps);
        gtPlusIONifti(float px, float py, float pz, float pt, float pr, float ps, float pp);
        gtPlusIONifti(float px, float py, float pz, float pt, float pr, float ps, float pp, float pq);

        virtual ~gtPlusIONifti() {}

        virtual void printInfo(std::ostream& os);

        virtual bool exportArray(const hoNDArray<short>& a, const std::string& filename) { return this->exportArrayImpl(a, filename); }
        virtual bool exportArray(const hoNDArray<unsigned short>& a, const std::string& filename) { return this->exportArrayImpl(a, filename); }
        virtual bool exportArray(const hoNDArray<int>& a, const std::string& filename) { return this->exportArrayImpl(a, filename); }
        virtual bool exportArray(const hoNDArray<unsigned int>& a, const std::string& filename) { return this->exportArrayImpl(a, filename); }
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
                std::string filenameData = filename;
                filenameData.append(".nii");

                HeaderType header;
                GADGET_CHECK_RETURN_FALSE(this->array2Header(a, header));

                gtPlusIOWorker ioworker(filenameData, false);

                GADGET_CHECK_RETURN_FALSE(ioworker.open());
                GADGET_CHECK_RETURN_FALSE(ioworker.write(reinterpret_cast<const char*>(&header), sizeof(HeaderType)));
                GADGET_CHECK_RETURN_FALSE(ioworker.write(reinterpret_cast<const char*>(a.begin()), a.get_number_of_bytes()));
                GADGET_CHECK_RETURN_FALSE(ioworker.close());
            }
            catch(...)
            {
                GADGET_ERROR_MSG("Errors in gtPlusIONifti::exportArrayImpl(const hoNDArray<T>& a, const std::string& filename) ... ");
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

                std::string filenameData = filename;
                filenameData.append(".nii");

                gtPlusIOWorker ioworker(filenameData, true);
                GADGET_CHECK_RETURN_FALSE(ioworker.open());
                GADGET_CHECK_RETURN_FALSE(ioworker.read(reinterpret_cast<char*>(&header), sizeof(HeaderType)));

                GADGET_CHECK_RETURN_FALSE(this->header2Array(a, header));

                GADGET_CHECK_RETURN_FALSE(ioworker.read(reinterpret_cast<char*>(a.begin()), a.get_number_of_bytes()));
                GADGET_CHECK_RETURN_FALSE(ioworker.close());
            }
            catch(...)
            {
                GADGET_ERROR_MSG("Errors in gtPlusIONifti::importArrayImpl(const hoNDArray<T>& a, const std::string& filename) ... ");
                return false;
            }

            return true;
        }

        template <typename T, unsigned int D> 
        bool exportImage(const hoNDImage<T,D>& a, const std::string& filename)
        {
            try
            {
                std::string filenameData = filename;
                filenameData.append(".nii");

                HeaderType header;
                std::vector<nifti1_extension> extenInfo;
                GADGET_CHECK_RETURN_FALSE(this->image2Header(a, header, extenInfo));

                gtPlusIOWorker ioworker(filenameData, false);

                GADGET_CHECK_RETURN_FALSE(ioworker.open());
                GADGET_CHECK_RETURN_FALSE(ioworker.write(reinterpret_cast<const char*>(&header), sizeof(HeaderType)));

                unsigned int ii;
                for ( ii=0; ii<extenInfo.size(); ii++ )
                {
                    GADGET_CHECK_RETURN_FALSE(ioworker.write(reinterpret_cast<const char*>(&extenInfo[ii]), extenInfo[ii].esize));
                }

                GADGET_CHECK_RETURN_FALSE(ioworker.write(reinterpret_cast<const char*>(a.begin()), a.get_number_of_bytes()));
                GADGET_CHECK_RETURN_FALSE(ioworker.close());
            }
            catch(...)
            {
                GADGET_ERROR_MSG("Errors in gtPlusIONifti::exportImage(const hoNDImage<T,D>& a, const std::string& filename) ... ");
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
                std::vector<nifti1_extension> extenInfo;

                std::string filenameData = filename;
                filenameData.append(".nii");

                gtPlusIOWorker ioworker(filenameData, true);
                GADGET_CHECK_RETURN_FALSE(ioworker.open());

                // read header
                GADGET_CHECK_RETURN_FALSE(ioworker.read(reinterpret_cast<char*>(&header), sizeof(HeaderType)));

                size_t ind = sizeof(HeaderType);

                if ( header.extension[0] > 0 )
                {
                    while ( ind < header.vox_offset )
                    {
                        nifti1_extension extenInfoObj;

                        // read extension info
                        GADGET_CHECK_RETURN_FALSE(ioworker.read(reinterpret_cast<char*>(&(extenInfoObj.esize)), sizeof(int)));
                        GADGET_CHECK_RETURN_FALSE(ioworker.read(reinterpret_cast<char*>(&(extenInfoObj.ecode)), sizeof(int)));

                        extenInfoObj.edata = new char[ extenInfoObj.esize - 2*sizeof(int) ];
                        GADGET_CHECK_RETURN_FALSE(ioworker.read(extenInfoObj.edata, extenInfoObj.esize-2*sizeof(int)));

                        extenInfo.push_back( extenInfoObj );
                        extenInfoObj.edata = NULL;

                        ind += extenInfoObj.esize;
                    }
                }

                GADGET_CHECK_RETURN_FALSE(this->header2Image(a, header, extenInfo));

                GADGET_CHECK_RETURN_FALSE(ioworker.seek(header.vox_offset));
                GADGET_CHECK_RETURN_FALSE(ioworker.read(reinterpret_cast<char*>(a.begin()), a.get_number_of_bytes()));
                GADGET_CHECK_RETURN_FALSE(ioworker.close());
            }
            catch(...)
            {
                GADGET_ERROR_MSG("Errors in gtPlusIONifti::importArray(hoNDImage<T,D>& a, const std::string& filename) ... ");
                return false;
            }

            return true;
        }

    protected:

        template <typename T> bool array2Header(const hoNDArray<T>& a, HeaderType& header);
        template <typename T> bool header2Array(hoNDArray<T>& a, const HeaderType& header);

        template <typename T, unsigned int D> bool image2Header(const hoNDImage<T, D>& a, HeaderType& header, std::vector<nifti1_extension>& extenInfo);
        template <typename T, unsigned int D> bool header2Image(hoNDImage<T, D>& a, const HeaderType& header, const std::vector<nifti1_extension>& extenInfo);
    };

    template <typename T> 
    bool gtPlusIONifti::array2Header(const hoNDArray<T>& a, HeaderType& header)
    {
        try
        {
            header.sizeof_hdr = 348;
            memset(header.data_type, '\0', 10);
            memset(header.db_name, '\0', 18);
            header.extents = 0;
            header.session_error = 0;
            header.regular = 0;

            header.dim_info = NIFTI_SLICE_UNKNOWN;

            // dimension
            size_t NDim = a.get_number_of_dimensions();

            header.dim[0] = (short)(NDim);
            header.dim[1] = (short)(a.get_size(0));

            if ( NDim > 1 )
                header.dim[2] = (short)(a.get_size(1));
            else
                header.dim[2] = 1;

            if ( NDim > 2 )
                header.dim[3] = (short)(a.get_size(2));
            else
                header.dim[3] = 1;

            if ( NDim > 3 )
                header.dim[4] = (short)(a.get_size(3));
            else
                header.dim[4] = 1;

            if ( NDim > 4 )
                header.dim[5] = (short)(a.get_size(4));
            else
                header.dim[5] = 1;

            if ( NDim > 5 )
                header.dim[6] = (short)(a.get_size(5));
            else
                header.dim[6] = 1;

            if ( NDim > 6 )
                header.dim[7] = (short)(a.get_size(6));
            else
                header.dim[7] = 1;

            if ( NDim > 7 )
                header.intent_p1 = (short)(a.get_size(7));
            else
                header.intent_p1 = 1;

            if ( NDim > 8 )
                header.intent_p2 = (short)(a.get_size(8));
            else
                header.intent_p2 = 1;

            if ( NDim > 9 )
                header.intent_p3 = (short)(a.get_size(9));
            else
                header.intent_p3 = 1;

            header.intent_code = 0;

            std::string rttiID = std::string(typeid(T).name());
            header.datatype = (short)getDataTypeFromRTTI(rttiID);
            header.bitpix = (short)(8*sizeof(T));

            header.slice_start = 0;

            // since the NDArray does not carry the pixel spacing
            header.pixdim[0] = 0;
            if ( pixelSize_.size() > 1 )
                header.pixdim[1] = pixelSize_[0];
            if ( pixelSize_.size() > 2 )
                header.pixdim[2] = pixelSize_[1];
            if ( pixelSize_.size() > 3 )
                header.pixdim[3] = pixelSize_[2];
            if ( pixelSize_.size() > 4 )
                header.pixdim[4] = pixelSize_[3];
            if ( pixelSize_.size() > 5 )
                header.pixdim[5] = pixelSize_[4];
            if ( pixelSize_.size() > 6 )
                header.pixdim[6] = pixelSize_[5];
            if ( pixelSize_.size() > 7 )
                header.pixdim[7] = pixelSize_[6];

            header.vox_offset = sizeof(nifti_1_header);

            header.scl_slope = 1;
            header.scl_inter = 0;

            header.slice_end = header.dim[3]-1;

            header.slice_code = NIFTI_SLICE_UNKNOWN;

            header.xyzt_units = SPACE_TIME_TO_XYZT(NIFTI_UNITS_MM, NIFTI_UNITS_MSEC);

            header.cal_max = 0;
            header.cal_min = 0;

            header.slice_duration = 0;

            header.toffset = 0;

            header.glmax = 0;
            header.glmin = 0;

            memset(header.descrip, '\0', 80);
            std::sprintf(header.descrip, "%s", " Gadgetron NDArray v2.5 ");

            memset(header.aux_file, '\0', 24);

            header.qform_code = NIFTI_XFORM_UNKNOWN; // since it is a NDArray
            header.sform_code = NIFTI_XFORM_UNKNOWN;

            header.quatern_b = 0;
            header.quatern_c = 0;
            header.quatern_d = 0;
            header.qoffset_x = 0;
            header.qoffset_y = 0;
            header.qoffset_z = 0;

            header.srow_x[0] = 1; header.srow_x[1] = 0; header.srow_x[2] = 0; header.srow_x[3] = 0;
            header.srow_y[0] = 0; header.srow_y[1] = 1; header.srow_y[2] = 0; header.srow_y[3] = 0;
            header.srow_z[0] = 0; header.srow_z[1] = 0; header.srow_z[2] = 1; header.srow_z[3] = 0;

            memset(header.intent_name, '\0', 16);
            std::sprintf(header.intent_name, "%s", "Pixel Values");

            header.magic[0] = 'n';
            header.magic[1] = '+';
            header.magic[2] = '1';
            header.magic[3] = '\0';

            header.extension[0] = 0; // for a ND array, no extension info is appended
            header.extension[1] = 0;
            header.extension[2] = 0;
            header.extension[3] = 0;
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in gtPlusIONifti::array2Header(const hoNDArray<T>& a, HeaderType& header) ... ");
            return false;
        }

        return true;
    }

    template <typename T> 
    bool gtPlusIONifti::header2Array(hoNDArray<T>& a, const HeaderType& header)
    {
        try
        {
            std::string rttiID = std::string(typeid(T).name());
            GADGET_CHECK_RETURN_FALSE(rttiID==getRTTIFromDataType( (GtDataType)header.datatype));

            std::vector<size_t> dim(header.dim[0]);
            size_t ii;
            for ( ii=0; ii<dim.size(); ii++ )
            {
                if ( ii == 7 )
                {
                    dim[ii] = (size_t)header.intent_p1;
                }
                else if ( ii == 8 )
                {
                    dim[ii] = (size_t)header.intent_p2;
                }
                else if ( ii == 9 ) 
                {
                    dim[ii] = (size_t)header.intent_p3;
                }
                else
                {
                    dim[ii] = header.dim[ii+1];
                }
            }

            pixelSize_.resize(dim.size());
            for ( ii=0; ii<dim.size(); ii++ )
            {
                if ( ii < 7 )
                {
                    pixelSize_[ii] = header.pixdim[ii+1];
                }
            }

            a.create(&dim);
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in gtPlusIONifti::header2Array(hoNDArray<T>& a, const HeaderType& header) ... ");
            return false;
        }

        return true;
    }

    template <typename T, unsigned int D> 
    bool gtPlusIONifti::image2Header(const hoNDImage<T,D>& a, HeaderType& header, std::vector<nifti1_extension>& extenInfo)
    {
        typedef typename hoNDImage<T,D>::coord_type coord_type;

        char* buf = NULL;
        coord_type* originBuf = NULL;
        coord_type* axisBuf = NULL;
        try
        {
            extenInfo.clear();
            extenInfo.resize(3);

            // serialize the attributes
            buf = NULL;
            size_t len(0);
            a.attrib_.serialize(buf, len);

            extenInfo[0].esize = sizeof(int) + sizeof(int) + len;
            extenInfo[0].ecode = 0; // user defined type
            extenInfo[0].edata = buf;

            // origin
            originBuf = new coord_type[D];
            unsigned int ii;
            for ( ii=0; ii<D; ii++ )
            {
                originBuf[ii] = a.get_origin(ii);
            }

            extenInfo[1].esize = sizeof(int) + sizeof(int) + sizeof(coord_type)*D;
            extenInfo[1].ecode = 0; // user defined type
            extenInfo[1].edata = reinterpret_cast<char*>(originBuf);

            // axis
            axisBuf = new coord_type[D*D];
            unsigned int jj;
            for ( ii=0; ii<D; ii++ )
            {
                for ( jj=0; jj<D; jj++ )
                {
                    axisBuf[ii*D+jj] = a.get_axis(ii, jj);
                }
            }

            extenInfo[2].esize = sizeof(int) + sizeof(int) + sizeof(coord_type)*D*D;
            extenInfo[2].ecode = 0; // user defined type
            extenInfo[2].edata = reinterpret_cast<char*>(axisBuf);

            // header info
            header.sizeof_hdr = 348;
            memset(header.data_type, '\0', 10);
            memset(header.db_name, '\0', 18);
            header.extents = 0;
            header.session_error = 0;
            header.regular = 0;

            header.dim_info = NIFTI_SLICE_UNKNOWN;

            // dimension
            size_t NDim = a.get_number_of_dimensions();

            header.dim[0] = (short)(NDim);
            header.dim[1] = (short)(a.get_size(0));

            if ( NDim > 1 )
                header.dim[2] = (short)(a.get_size(1));
            else
                header.dim[2] = 1;

            if ( NDim > 2 )
                header.dim[3] = (short)(a.get_size(2));
            else
                header.dim[3] = 1;

            if ( NDim > 3 )
                header.dim[4] = (short)(a.get_size(3));
            else
                header.dim[4] = 1;

            if ( NDim > 4 )
                header.dim[5] = (short)(a.get_size(4));
            else
                header.dim[5] = 1;

            if ( NDim > 5 )
                header.dim[6] = (short)(a.get_size(5));
            else
                header.dim[6] = 1;

            if ( NDim > 6 )
                header.dim[7] = (short)(a.get_size(6));
            else
                header.dim[7] = 1;

            if ( NDim > 7 )
                header.intent_p1 = (short)(a.get_size(7));
            else
                header.intent_p1 = 1;

            if ( NDim > 8 )
                header.intent_p2 = (short)(a.get_size(8));
            else
                header.intent_p2 = 1;

            if ( NDim > 9 )
                header.intent_p3 = (short)(a.get_size(9));
            else
                header.intent_p3 = 1;

            header.intent_code = 0;

            std::string rttiID = std::string(typeid(T).name());
            header.datatype = (short)getDataTypeFromRTTI(rttiID);
            header.bitpix = (short)(8*sizeof(T));

            header.slice_start = 0;

            // since the NDArray does not carry the pixel spacing
            header.pixdim[0] = 0;
            if ( pixelSize_.size() > 1 )
                header.pixdim[1] = pixelSize_[0];
            if ( pixelSize_.size() > 2 )
                header.pixdim[2] = pixelSize_[1];
            if ( pixelSize_.size() > 3 )
                header.pixdim[3] = pixelSize_[2];
            if ( pixelSize_.size() > 4 )
                header.pixdim[4] = pixelSize_[3];
            if ( pixelSize_.size() > 5 )
                header.pixdim[5] = pixelSize_[4];
            if ( pixelSize_.size() > 6 )
                header.pixdim[6] = pixelSize_[5];
            if ( pixelSize_.size() > 7 )
                header.pixdim[7] = pixelSize_[6];

            header.vox_offset = sizeof(nifti_1_header) + extenInfo[0].esize + extenInfo[1].esize + extenInfo[2].esize;

            header.scl_slope = 1;
            header.scl_inter = 0;

            header.slice_end = header.dim[3]-1;

            header.slice_code = NIFTI_SLICE_UNKNOWN;

            header.xyzt_units = SPACE_TIME_TO_XYZT(NIFTI_UNITS_MM, NIFTI_UNITS_MSEC);

            header.cal_max = 0;
            header.cal_min = 0;

            header.slice_duration = 0;

            header.toffset = 0;

            header.glmax = 0;
            header.glmin = 0;

            memset(header.descrip, '\0', 80);
            std::sprintf(header.descrip, "%s", " Gadgetron NDImage v2.5 ");

            memset(header.aux_file, '\0', 24);

            header.qform_code = NIFTI_XFORM_SCANNER_ANAT;
            header.sform_code = NIFTI_XFORM_UNKNOWN;

            float quat[4];
            a.get_image_orientation(quat);

            header.quatern_b = quat[1];
            header.quatern_c = quat[2];
            header.quatern_d = quat[3];
            header.qoffset_x = a.get_image_position(0);
            header.qoffset_y = a.get_image_position(1);
            header.qoffset_z = a.get_image_position(2);

            header.srow_x[0] = 1; header.srow_x[1] = 0; header.srow_x[2] = 0; header.srow_x[3] = 0;
            header.srow_y[0] = 0; header.srow_y[1] = 1; header.srow_y[2] = 0; header.srow_y[3] = 0;
            header.srow_z[0] = 0; header.srow_z[1] = 0; header.srow_z[2] = 1; header.srow_z[3] = 0;

            memset(header.intent_name, '\0', 16);
            std::sprintf(header.intent_name, "%s", "Pixel Values");

            header.magic[0] = 'n';
            header.magic[1] = '+';
            header.magic[2] = '1';
            header.magic[3] = '\0';

            header.extension[0] = 1;
            header.extension[1] = 0;
            header.extension[2] = 0;
            header.extension[3] = 0;
        }
        catch(...)
        {
            delete [] buf;
            delete [] originBuf;
            delete [] axisBuf;

            GADGET_ERROR_MSG("Errors in gtPlusIONifti::image2Header(const hoNDImage<T>& a, HeaderType& header, std::vector<nifti1_extension>& extenInfo) ... ");
            return false;
        }

        return true;
    }

    template <typename T, unsigned int D> 
    bool gtPlusIONifti::header2Image(hoNDImage<T,D>& a, const HeaderType& header, const std::vector<nifti1_extension>& extenInfo)
    {
        try
        {
            std::string rttiID = std::string(typeid(T).name());
            GADGET_CHECK_RETURN_FALSE(rttiID==getRTTIFromDataType( (GtDataType)header.datatype));

            std::vector<size_t> dim(header.dim[0]);

            if ( D != dim.size() ) return false;

            size_t ii;
            for ( ii=0; ii<dim.size(); ii++ )
            {
                if ( ii == 7 )
                {
                    dim[ii] = header.intent_p1;
                }
                else if ( ii == 8 )
                {
                    dim[ii] = header.intent_p2;
                }
                else if ( ii == 9 ) 
                {
                    dim[ii] = header.intent_p3;
                }
                else
                {
                    dim[ii] = header.dim[ii+1];
                }
            }

            a.create(dim);

            for ( ii=0; ii<dim.size(); ii++ )
            {
                if ( ii < 7 )
                {
                    a.set_pixel_size(ii, header.pixdim[ii+1]);
                }
            }

            // set the image position and orientation patient
            a.set_image_position(0, header.qoffset_x);
            a.set_image_position(1, header.qoffset_y);
            a.set_image_position(2, header.qoffset_z);

            float quat[4];

            quat[1] = header.quatern_b;
            quat[2] = header.quatern_c;
            quat[3] = header.quatern_d;
            float v = header.quatern_b*header.quatern_b - header.quatern_c*header.quatern_c -header.quatern_d*header.quatern_d;
            if ( v < 1 )
            {
                quat[0] = (float)std::sqrt(1.0 - v);
            }
            else
            {
                quat[0] = 0;
            }

            a.set_image_orientation(quat);

            // get origin
            std::vector<float> buf(D);
            memcpy(&buf[0], extenInfo[1].edata, D*sizeof(float));

            for ( ii=0; ii<D; ii++ )
            {
                a.set_origin(ii, buf[ii]);
            }

            // get axis
            buf.resize(D*D);
            memcpy(&buf[0], extenInfo[2].edata, D*D*sizeof(float));

            unsigned int jj;
            for ( ii=0; ii<D; ii++ )
            {
                typename hoNDImage<T,D>::coord_type v(0);
                typename hoNDImage<T,D>::coord_type mag(0);

                for ( jj=0; jj<D; jj++ )
                {
                    v = buf[ii*D+jj];
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

            // get attributes
            a.attrib_.deserialize(extenInfo[2].edata, extenInfo[2].esize-2*sizeof(int));
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in gtPlusIONifti::header2Image(hoNDImage<T,D>& a, const HeaderType& header, const std::vector<nifti1_extension>& extenInfo) ... ");
            return false;
        }

        return true;
    }

}}
