#pragma once

#if defined (WIN32)
#ifdef __BUILD_GADGETRON_TINYXML__
#define EXPORTTINYXML __declspec(dllexport)
#else
#define EXPORTTINYXML __declspec(dllimport)
#endif
#else
#define EXPORTTINYXML
#endif