/** \file       cpuOperatorExport.h
    \brief      Implement windows export/import for gadgetron cpu operator toolbox
    \author     Hui Xue
*/

#pragma once

#if defined (WIN32)
    #if defined (__BUILD_GADGETRON_CPUOPERATORS__) || defined (gadgetron_toolbox_cpuoperator_EXPORTS)
        #define EXPORTCPUOPERATOR __declspec(dllexport)
    #else
        #define EXPORTCPUOPERATOR __declspec(dllimport)
    #endif
#else
    #define EXPORTCPUOPERATOR
#endif
