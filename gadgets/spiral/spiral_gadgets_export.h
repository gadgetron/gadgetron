#pragma once

#if defined (WIN32)
#ifdef __BUILD_GADGETRON_SPIRAL__
#define EXPORTGADGETSSPIRAL __declspec(dllexport)
#else
#define EXPORTGADGETSSPIRAL __declspec(dllimport)
#endif
#else
#define EXPORTGADGETSSPIRAL
#endif
