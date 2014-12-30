#ifndef _MRI_CORE_EXPORT_H_
#define _MRI_CORE_EXPORT_H_

#if defined (WIN32)
    #ifdef BUILD_TOOLBOX_STATIC
        #define EXPORTMRICORE
    #else
        #if defined (__BUILD_GADGETRON_MRI_CORE__) || defined (mri_core_EXPORTS)
            #define EXPORTMRICORE __declspec(dllexport)
        #else
            #define EXPORTMRICORE __declspec(dllimport)
        #endif
    #endif
#else
    #define EXPORTMRICORE
#endif

#endif /* _MRI_CORE_EXPORT_H_ */
