#ifndef GADGETRONEXPORT_H
#define GADGETRONEXPORT_H
#pragma once

#if defined (WIN32)
#ifdef __BUILD_GADGETS__
#define GADGETEXPORT __declspec(dllexport)
#else
#define GADGETEXPORT __declspec(dllimport)
#endif
#else
#define GADGETEXPORT
#endif

//In CPP file add this macro add the end
#define GADGETRON_LOADABLE_FACTORY_DECLARE(CLASS, COMPONENT)	\
extern "C" GADGETEXPORT CLASS * make_##COMPONENT (void);        \
CLASS * make_##COMPONENT (void)       				\
{							       	\
  return new COMPONENT;                                         \
}                                                              


#endif
