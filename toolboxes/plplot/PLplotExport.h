/** \file       PLplotExport.h
    \brief      Implement windows export/import for PLplot gadgetron utility toolbox
    \author     Hui Xue
*/

#pragma once

#if defined (WIN32)
    #if defined (__BUILD_GADGETRON_TOOLBOX_PLPLOT__) || defined (gadgetron_toolbox_plplot_EXPORTS)
        #define EXPORTGTPLPLOT __declspec(dllexport)
    #else
        #define EXPORTGTPLPLOT __declspec(dllimport)
    #endif
#else
    #define EXPORTGTPLPLOT
#endif
