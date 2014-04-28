/** \file       GtPlusIOExport.h
    \brief      Implement export/import for GtPlus toolbox
    \author     Hui Xue
*/

#pragma once

#if defined (WIN32)
    #ifdef BUILD_TOOLBOX_STATIC
        #define EXPORTGTPLUSIO 
    #else
        #if defined (__BUILD_GADGETRON_PLUS__) || defined (gtplus_io_EXPORTS)
            #define EXPORTGTPLUSIO __declspec(dllexport)
        #else
            #define EXPORTGTPLUSIO __declspec(dllimport)
        #endif
    #endif
#else
    #define EXPORTGTPLUSIO
#endif
