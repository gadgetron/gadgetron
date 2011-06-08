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

#if defined (WIN32)
#ifdef __BUILD_GADGETRON_CORE__
#define EXPORTGADGETSCORE __declspec(dllexport)
#else
#define EXPORTGADGETSCORE __declspec(dllimport)
#endif
#else
#define EXPORTGADGETSCORE
#endif

#if defined (WIN32)
#ifdef __BUILD_GADGETRON_GRAPPA__
#define EXPORTGADGETSGRAPPA __declspec(dllexport)
#else
#define EXPORTGADGETSGRAPPA __declspec(dllimport)
#endif
#else
#define EXPORTGADGETSGRAPPA
#endif

#if defined (WIN32)
#ifdef __BUILD_GADGETRON_CGSENSE__
#define EXPORTGADGETSCGSENSE __declspec(dllexport)
#else
#define EXPORTGADGETSCGSENSE __declspec(dllimport)
#endif
#else
#define EXPORTGADGETSCGSENSE
#endif
