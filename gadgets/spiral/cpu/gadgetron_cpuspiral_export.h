#pragma once

#if defined (WIN32)
#ifdef __BUILD_GADGETRON_CPUSPIRAL__
#define EXPORTGADGETS_CPUSPIRAL __declspec(dllexport)
#else
#define EXPORTGADGETS_CPUSPIRAL __declspec(dllimport)
#endif
#else
#define EXPORTGADGETS_CPUSPIRAL
#endif
