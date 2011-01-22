#ifndef __DEFINE_DLLEXPORT__
#define __DEFINE_DLLEXPORT__

#if defined (WIN32) || defined (_WIN32)

#ifdef GADGETS_BUILD_DLL
#define DLLEXPORT __declspec(dllexport)
#else
#define DLLEXPORT __declspec(dllimport)
#endif

#else
/* In Linux, DLLEXPORT is ignored */
#define DLLEXPORT
#endif

#endif
