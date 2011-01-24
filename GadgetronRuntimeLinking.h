#ifndef GADGETRONRUNTIMELINKING_H
#define GADGETRONRUNTIMELINKING_H

#include "ace/OS_NS_stdlib.h"
#include "ace/OS_NS_string.h"
#include "ace/OS_NS_stdio.h"
#include "ace/DLL.h"
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

template <class T> inline T* GadgetronLoadComponent(const char* DLL, const char* component_name)
{
  ACE_DLL dll;

  ACE_TCHAR dllname[1024];
  ACE_OS::sprintf(dllname, "%s%s",ACE_DLL_PREFIX, DLL);

  ACE_TCHAR factoryname[1024];
  ACE_OS::sprintf(factoryname, "make_%s", component_name);

  int retval = dll.open (dllname);

  if (retval != 0)
    ACE_ERROR_RETURN ((LM_ERROR,
                       "%p, ---%s---, %s\n",
                       "dll.open", dllname, dll.error()),
                      0);

  //Function pointer
  typedef T* (*ComponentCreator) (void);
  
  
  void *void_ptr = dll.symbol (factoryname);
  ptrdiff_t tmp = reinterpret_cast<ptrdiff_t> (void_ptr);
  ComponentCreator cc = reinterpret_cast<ComponentCreator> (tmp);

  if (cc == 0) {
    ACE_ERROR_RETURN ((LM_ERROR,
		       "%p,  ---%s---, %s\n",
		       "dll.symbol", factoryname, dll.error()),
		      0);
  }


  T* c = cc();
  
  if (!c) {
    GADGET_DEBUG1("Failed to create component using factory\n");
    return 0;
  }

  dll.close ();

  return c;
}


#endif //GADGETRONRUNTIMELINKING_H
