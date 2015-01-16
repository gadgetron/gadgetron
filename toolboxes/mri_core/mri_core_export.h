/** \file       mri_core_export.h
    \brief      Implement windows export/import for mri_core toolbox
    \author     Hui Xue
*/

#pragma once

#if defined (WIN32)
    #if defined (__BUILD_GADGETRON_TOOLBOX_MRI_CORE__) || defined (gadgetron_toolbox_mri_core_EXPORTS)
        #define EXPORTMRICORE __declspec(dllexport)
    #else
        #define EXPORTMRICORE __declspec(dllimport)
    #endif
#else
    #define EXPORTMRICORE
#endif
