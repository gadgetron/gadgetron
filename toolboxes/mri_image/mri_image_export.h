#ifndef _MRI_IMAGE_EXPORT_H_
#define _MRI_IMAGE_EXPORT_H_

#if defined (WIN32)
    #if defined (__BUILD_GADGETRON_MRI_IMAGE__) || defined (gadgetron_toolbox_mri_image_EXPORTS)
        #define EXPORTMRIIMAGE __declspec(dllexport)
    #else
        #define EXPORTMRIIMAGE __declspec(dllimport)
    #endif
#else
    #define EXPORTMRIIMAGE
#endif

#endif /* _MRI_IMAGE_EXPORT_H_ */
