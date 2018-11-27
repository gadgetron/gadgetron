#pragma once

#if defined (WIN32)
#if defined (__BUILD_GADGETRON_GADGETCORE__) || defined (gadgetron_gadgetcore_EXPORTS)
#define EXPORTGADGETCORE __declspec(dllexport)
#else
#define EXPORTGADGETCORE __declspec(dllimport)
#endif
#else
#define EXPORTGADGETCORE
#endif


