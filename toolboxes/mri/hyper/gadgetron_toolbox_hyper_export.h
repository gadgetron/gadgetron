#pragma once

#if defined (WIN32)
#if defined (__BUILD_GADGETRON_TOOLBOX_HYPER__)
#define EXPORTHYPER __declspec(dllexport)
#else
#define EXPORTHYPER __declspec(dllimport)
#endif
#else
#define EXPORTHYPER
#endif