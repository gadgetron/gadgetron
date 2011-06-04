#pragma once

#if defined (WIN32)
#ifdef __BUILD_GADGETRON_GPUCORE__
#define EXPORTGPUCORE __declspec(dllexport)
#else
#define EXPORTGPUCORE __declspec(dllimport)
#endif
#else
#define EXPORTGPUCORE
#endif

#if defined (WIN32)
#ifdef __BUILD_GADGETRON_GPUCG__
#define EXPORTGPUCG __declspec(dllexport)
#else
#define EXPORTGPUCG __declspec(dllimport)
#endif
#else
#define EXPORTGPUCG
#endif

#if defined (WIN32)
#ifdef __BUILD_GADGETRON_GPUPMRI__
#define EXPORTGPUPMRI __declspec(dllexport)
#else
#define EXPORTGPUPMRI __declspec(dllimport)
#endif
#else
#define EXPORTGPUPMRI
#endif

#if defined (WIN32)
#ifdef __BUILD_GADGETRON_GPUNFFT__
#define EXPORTGPUNFFT __declspec(dllexport)
#else
#define EXPORTGPUNFFT __declspec(dllimport)
#endif
#else
#define EXPORTGPUNFFT
#endif

#if defined (WIN32)
#ifdef __BUILD_GADGETRON_HOSTUTILS__
#define EXPORTHOSTUTILS __declspec(dllexport)
#else
#define EXPORTHOSTUTILS __declspec(dllimport)
#endif
#else
#define EXPORTHOSTUTILS
#endif
