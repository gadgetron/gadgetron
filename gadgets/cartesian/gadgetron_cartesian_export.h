#pragma once

#if defined (WIN32)
#ifdef __BUILD_GADGETRON_CARTESIAN__
#define EXPORTGADGETS_CARTESIAN __declspec(dllexport)
#else
#define EXPORTGADGETS_CARTESIAN __declspec(dllimport)
#endif
#else
#define EXPORTGADGETS_CARTESIAN
#endif
