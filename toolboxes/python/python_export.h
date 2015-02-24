#ifndef PYTHON_EXPORT_H_
#define PYTHON_EXPORT_H_

#if defined (WIN32)
    #if defined (__BUILD_GADGETRON_PYTHON__) || defined (gadgetron_toolbox_python_EXPORTS)
        #define EXPORTPYTHON __declspec(dllexport)
    #else
        #define EXPORTPYTHON __declspec(dllimport)
    #endif
#else
    #define EXPORTPYTHON
#endif

#endif
