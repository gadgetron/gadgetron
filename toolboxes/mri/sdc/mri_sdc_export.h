/**
 * \file mri_sdc_export.h
 * \brief Export definitions for sampling density compensation.
 */

#ifndef MRI_SDC_EXPORT_H_
#define MRI_SDC_EXPORT_H_

#if defined (WIN32)
    #if defined (__BUILD_MRI_NONCARTESIAN__)
        #define EXPORTSDC __declspec(dllexport)
    #else
        #define EXPORTSDC __declspec(dllimport)
    #endif
#else
    #define EXPORTSDC
#endif

#endif  // MRI_SDC_EXPORT_H_ 