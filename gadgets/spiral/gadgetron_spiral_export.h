#pragma once

#if defined (WIN32)
#ifdef __BUILD_GADGETRON_SPIRAL__
#define EXPORTGADGETS_SPIRAL __declspec(dllexport)
#else
#define EXPORTGADGETS_SPIRAL __declspec(dllimport)
#endif
#else
#define EXPORTGADGETS_SPIRAL
#endif
