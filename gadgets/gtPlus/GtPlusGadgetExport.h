/** \file   GtPlusGadgetExport.h
    \brief  The windows export/import definition for the GtPlus reconstruction gadget
    \author Hui Xue
*/

#pragma once

#if defined (WIN32)
    #if defined (__BUILD_GADGETS__) || defined (gadgetronPlus_EXPORTS)
        #define EXPORTGTPLUSGADGET __declspec(dllexport)
    #else
        #define EXPORTGTPLUSGADGET __declspec(dllimport)
    #endif
#else
    #define EXPORTGTPLUSGADGET
#endif
