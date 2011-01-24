#ifndef GADGETRONRUNTIMELINKING_H
#define GADGETRONRUNTIMELINKING_H

#include "ace/OS_NS_stdlib.h"
#include "ace/OS_NS_string.h"
#include "ace/OS_NS_stdio.h"
#include "ace/DLL_Manager.h"

#include "ace/Log_Msg.h"
#include "Gadgetron.h"


//In header file add this macro
#define GADGETRON_LOADABLE_DECLARE(COMPONENT)                   \
  void *operator new (size_t bytes);                            \
  void operator delete (void *ptr);                             \

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



#endif //GADGETRONRUNTIMELINKING_H
