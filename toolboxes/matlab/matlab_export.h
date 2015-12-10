/** \file matlab_export.h
    \brief Required definitions for Windows, importing/exporting dll symbols 
*/

#ifndef MATLAB_EXPORT_H_
#define MATLAB_EXPORT_H_

#if defined (WIN32)
    #if defined (__BUILD_GADGETRON_TOOLBOX_MATLAB__) || defined (gadgetron_toolbox_matlab_EXPORTS)
        #define EXPORTMATLAB __declspec(dllexport)
    #else
        #define EXPORTMATLAB __declspec(dllimport)
    #endif
#else
    #define EXPORTMATLAB
#endif

#endif /* MATLAB_EXPORT_H_ */
