/** \file       GtPlusExport.h
    \brief      Implement windows export/import for GtPlus toolbox
    \author     Hui Xue
*/

#pragma once

#if defined (WIN32)
    #if defined (__BUILD_GADGETRON_PLUS__) || defined (gtplus_EXPORTS)
        #define EXPORTGTPLUS __declspec(dllexport)
    #else
        #define EXPORTGTPLUS __declspec(dllimport)
    #endif
#else
    #define EXPORTGTPLUS
#endif
