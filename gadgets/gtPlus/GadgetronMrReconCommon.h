#ifndef GADGETRONMRRECONCOMMON_H
#define GADGETRONMRRECONCOMMON_H

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

    // warning C4305: 'argument' : truncation
    #pragma warning( disable : 4305)

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

#endif  // GADGETRONMRRECONCOMMON_H
