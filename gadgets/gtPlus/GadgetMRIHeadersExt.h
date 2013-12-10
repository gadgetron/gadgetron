#ifndef GADGETMRIHEADERSEXT_H
#define GADGETMRIHEADERSEXT_H

#include "gadgetronMrRecon_export.h"
#include "GadgetronMrReconCommon.h"
#include "GadgetMRIHeaders.h"
#include "ismrmrd.h"
#include "core/basic/Common.h"
#include <algorithm/MrRecon/basic/MrReconNDArray.h>

#include <vector>

/** @name OS and compiler version */
//@{
#ifdef _WIN32
    // assume microsft visual c++ compiler if on windows
    #define GADGETRON_FTK_VISUAL_CPP
#elif defined WIN32
    #define GADGETRON_FTK_VISUAL_CPP
#elif defined WINDOWS
    #define GADGETRON_FTK_VISUAL_CPP
#else
    // not the visual studio, maybe gcc
    #define NOT_WIN32
    #define GADGETRON_FTK_DEPRECATED
#endif

#ifdef GADGETRON_FTK_VISUAL_CPP
    #if _MSC_VER >= 1300 // vc 7 or higher, only vc6 does not support template very well
        #define GADGETRON_FTK_TEMPLATE_SUPPORT
    #else
        #ifndef GADGETRON_FTK_OLD_VC_FLAG
            #define GADGETRON_FTK_OLD_VC_FLAG // vc 6 flag
        #endif
    #endif
#elif defined NOT_WIN32 // gcc or others
    #define GADGETRON_FTK_TEMPLATE_SUPPORT
#endif

// settings specific for microsoft compiler
#ifdef GADGETRON_FTK_VISUAL_CPP
    // disable warnings on 255 char debug symbols
    #pragma warning (disable : 4786)

    // disable warnings on exporting classes in DLL which has STL members
    #pragma warning (disable : 4251)

    // disable warnings on using 'this' in initializer list
    #pragma warning (disable : 4355)

    // disable warnings when specifying functions with a throw specifier
    #pragma warning( disable : 4290 )

    // disable warnings for implicit conversions
    //#pragma warning( disable : 4244 )

    // disable warnings for unknown pragma
    #pragma warning( disable : 4068 )
    
    // disable warnings for unsafe functions
    #pragma warning( disable : 4996 )

    // disable warnings for warning C4275: non dll-interface class 
    // 'std::_Complex_base<float>' used as base for dll-interface 
    //class 'std::complex<float>'
    #pragma warning( disable : 4275 )

    /// disable warning for constant conditional expression
    #pragma warning( disable : 4127)

    /// disable warning for unreachable code
    #pragma warning( disable : 4702)

    /// 'identifier' : decorated name length exceeded, name was truncated
    /// The decorated name was longer than the maximum the compiler allows (247), 
    /// and was truncated. To avoid this warning and the truncation, reduce the number of arguments or name length of identifiers used.
    #pragma warning( disable : 4503)

    #pragma warning( disable : 4267)
    #pragma warning( disable : 4244)
    #pragma warning( disable : 4996)

    // debug functionality
    // #include <crtdbg.h>

    // make code portable between VSS 6.0 and .NET
    #if _MSC_VER >= 1300 // check for .NET
    #define GADGETRON_FTK_DEPRECATED __declspec(deprecated)
    #else
    #define GADGETRON_FTK_DEPRECATED
    #endif

#endif
//@}

// -----------------------------------------------------------------
// info zone

enum PATRefScanMode
{
    PAT_REF_SCAN_UNDEFINED      = 0x01, // e.g. if no PAT is selected
    PAT_REF_SCAN_INPLACE        = 0x02, // sequence supplies inplace reference lines
    PAT_REF_SCAN_EXTRA          = 0x04, // sequence supplies extra reference lines
    PAT_REF_SCAN_PRESCAN        = 0x08, // sequence does not supply reference lines, the data must have been acquired with a previous measurement
    PAT_REF_SCAN_INTRINSIC_AVE  = 0x10, // The sequence contains intrinsic ref.lines due to sharing e.g. in the averages dimension
    PAT_REF_SCAN_INTRINSIC_REP  = 0x20, // The sequence contains intrinsic ref.lines due to sharing e.g. in the repetition or phases dimension (i.e., TSENSE)
    PAT_REF_SCAN_INTRINSIC_PHS  = 0x40, // The sequence contains intrinsic ref.lines due to sharing e.g. in the repetition or phases dimension (i.e., TSENSE)
    PAT_REF_SCAN_INPLACE_LET    = 0x80  // A single (L)ong (E)cho (T)rain acquires reference lines and imaging lines
};

struct LoopCounters
{
    ACE_UINT16 line;
    ACE_UINT16 acquisition;
    ACE_UINT16 slice;
    ACE_UINT16 partition;
    ACE_UINT16 echo;
    ACE_UINT16 phase;
    ACE_UINT16 repetition;
    ACE_UINT16 set;
    ACE_UINT16 segment;
    ACE_UINT16 channel;
};

#define MDH_FREEHDRPARA         4
#define MDH_FREEHDRPARAOFFSET   4

// aushIceProgramPara

// in the user_int
#define WIP_INDEX_TR                                MDH_FREEHDRPARAOFFSET+0
#define WIP_INDEX_TE                                MDH_FREEHDRPARAOFFSET+1
#define WIP_INDEX_FOV                               MDH_FREEHDRPARAOFFSET+2
#define WIP_INDEX_SliceThickness                    MDH_FREEHDRPARAOFFSET+3

// in the user_float
#define WIP_INDEX_BaseResolution                    0
#define WIP_INDEX_KernelSelection                   1
#define WIP_INDEX_MoCoRecon                         2
#define WIP_INDEX_AcceFactor                        3
#define WIP_INDEX_NumRepForMoCo                     4
#define WIP_INDEX_NumOfLine                         5

// -----------------------------------------------------------------

// [Col Line Cha Slice Partition Echo Phase Rep Set Seg]
#ifdef VDimMrRecon
    #undef VDimMrRecon
#endif // VDimMrRecon
#define  VDimMrRecon 10

#ifdef BufferLengthMrRecon
    #undef BufferLengthMrRecon
#endif // BufferLengthMrRecon
#define  BufferLengthMrRecon 2048

struct  EXPORTGADGETSMRRECON GadgetMessageImageExt : public ISMRMRD::ImageHeader
{
    // fields added to store the time_stamp and pmu_time_stamp for every incoming read-out line
    // if one line is not acquried, the corresponding time is -1
    std::vector<int>     time_stamps;
    std::vector<int>     pmu_time_stamps;

    GadgetMessageImageExt();
    ~GadgetMessageImageExt();

    void copy(GadgetMessageImageExt& aMessageImage);
    void set_matrix_size(unsigned int index, ACE_UINT16 size);
    void dump();
}; 

// [Col Line Cha Slice Partition Echo Phase Rep Set Seg]
//   0   1    2   3     4         5    6     7   8   9
// store a scan with 10 dimensions
struct  EXPORTGADGETSMRRECON GadgetMessageImageArray
{
    // size of the image array
    ACE_UINT16 matrix_size[10];

    // kspace center column number
    ACE_UINT16 kSpace_centre_col_no;
    // kspace max acquired col number
    ACE_UINT16 kSpace_max_acquired_col_no;

    // kspace center line number
    ACE_UINT16 kSpace_centre_line_no;
    // kspace max acquired line number
    ACE_UINT16 kSpace_max_acquired_line_no;

    // kspace center partition number
    ACE_UINT16 kSpace_centre_partition_no;
    // kspace max acquired partition number
    ACE_UINT16 kSpace_max_acquired_partition_no;

    // message information for every 2D image [Slice Partition Echo Phase Rep Set Seg]
    GadgetMessageImageExt* imageArray_;

    GadgetMessageImageArray();
    GadgetMessageImageArray(int aSize[10]);
    ~GadgetMessageImageArray();

    void resize(int aSize[10]);
    void copy(GadgetMessageImageArray& imageArray);
    int get_offset(int slc, int par, int eco, int phs, int rep, int set, int seg);
    void extractMessageImageArrayForSLC(int slc, GadgetMessageImageArray& imageArray);
    void extractMessageImageArrayForREP(int rep, GadgetMessageImageArray& imageArray);

    void dump();
};

struct EXPORTGADGETSMRRECON KSpaceBuffer
{
    typedef FTK_NAMESPACE_NAME::MrReconNDArray< std::complex<float>, VDimMrRecon > MrReconBufferType;

    // kspace data
    MrReconBufferType buffer_;

    // reference ACS data
    MrReconBufferType ref_;

    // other data, e.g. AIF data
    MrReconBufferType other_;

    // whether it is ipat or pat with seperate ref
    bool isIPAT;

    KSpaceBuffer();
    ~KSpaceBuffer();
};

#endif  //GADGETMRIHEADERSEXT_H
