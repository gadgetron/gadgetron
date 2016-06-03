#ifndef _MRI_CORE_EXPORT_H_
#define _MRI_CORE_EXPORT_H_

#if defined (WIN32)
    #if defined (__BUILD_GADGETRON_MRI_CORE__) || defined (gadgetron_toolbox_mri_core_EXPORTS)
        #define EXPORTMRICORE __declspec(dllexport)
    #else
        #define EXPORTMRICORE __declspec(dllimport)
    #endif
#else
    #define EXPORTMRICORE
#endif

#endif /* _MRI_CORE_EXPORT_H_ */
