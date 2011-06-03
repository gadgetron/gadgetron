#ifndef GADGETRONEXPORT_H
#define GADGETRONEXPORT_H

#if defined (WIN32) || defined (_WIN32)

#ifdef GADGETS_BUILD_DLL
#define DLLEXPORT __declspec(dllexport)
#else
#define DLLEXPORT __declspec(dllimport)
#endif

#else
/* In Linux, DLLEXPORT is ignored */
#define DLLEXPORT
#endif


//In header file add this macro
#define GADGETRON_LOADABLE_DECLARE(COMPONENT)                   \
  void *operator new (size_t bytes);                            \
  void operator delete (void *ptr);                             \
  void *operator new(size_t s, void * p) { return p; }

//In CPP file add this macro add the end
#define GADGETRON_LOADABLE_FACTORY_DECLARE(CLASS, COMPONENT)	\
extern "C" DLLEXPORT CLASS * make_##COMPONENT (void);           \
CLASS * make_##COMPONENT (void)       				\
{							       	\
  return new COMPONENT;                                         \
}                                                               \
void * COMPONENT ::operator new (size_t bytes)                  \
{                                                               \
  return ::new char[bytes];                                     \
}                                                               \
void COMPONENT ::operator delete (void *ptr)                    \
{                                                               \
  delete [] static_cast <char *> (ptr);                         \
} 


#endif
