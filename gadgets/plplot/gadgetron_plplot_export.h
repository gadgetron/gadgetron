/** \file   gadgetron_plplot_export.h
    \brief  The windows export/import definition for the gadgetron plplot related gadgets
    \author Hui Xue
*/

#pragma once

#if defined (WIN32)
    #if defined (__BUILD_GADGETRON_PLPLOT__) || defined (gadgetron_plplot_EXPORTS)
        #define EXPORTPLPLOTGADGET __declspec(dllexport)
    #else
        #define EXPORTPLPLOTGADGET __declspec(dllimport)
    #endif
#else
    #define EXPORTPLPLOTGADGET
#endif
