/** \file       EPIExport.h
    \brief      Implement windows export/import for EPI toolbox
    \author     Souheil Inati
*/

#pragma once

#if defined (WIN32)
    #if defined (__BUILD_GADGETRON_EPI__)
        #define EXPORTEPI __declspec(dllexport)
    #else
        #define EXPORTEPI __declspec(dllimport)
    #endif
#else
    #define EXPORTEPI
#endif
