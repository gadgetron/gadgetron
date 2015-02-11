#ifndef PYTHON_MATH_EXPORT_H_
#define PYTHON_MATH_EXPORT_H_

#if defined (WIN32)
    #if defined (__BUILD_GADGETRON_PYTHON_MATH__) || defined (gadgetron_toolbox_python_math_EXPORTS)
        #define EXPORTPYTHONMATH __declspec(dllexport)
    #else
        #define EXPORTPYTHONMATH __declspec(dllimport)
    #endif
#else
    #define EXPORTPYTHONMATH
#endif

#endif
