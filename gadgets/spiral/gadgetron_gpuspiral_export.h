#pragma once

#if defined (WIN32)
#ifdef __BUILD_GADGETRON_GPUSPIRAL__
#define EXPORTGADGETS_GPUSPIRAL __declspec(dllexport)
#else
#define EXPORTGADGETS_GPUSPIRAL __declspec(dllimport)
#endif
#else
#define EXPORTGADGETS_GPUSPIRAL
#endif
