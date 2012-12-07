#ifndef GADGETRON_H
#define GADGETRON_H

#include "ace/Log_Msg.h"

//#include "Gadget.h"
//#include "GadgetContainerMessage.h"

//Return messages
#define GADGET_FAIL -1
#define GADGET_OK    0


//MACROS FOR LOGGING
#define GADGET_DEBUG1(_fmt) \
  ACE_DEBUG( (LM_DEBUG, \
	      ACE_TEXT("[file %N, line %l] " _fmt)) ) 

#define GADGET_DEBUG2(_fmt, ...) \
  ACE_DEBUG( (LM_DEBUG, \
	      ACE_TEXT("[file %N, line %l] " _fmt),	\
	      __VA_ARGS__) )
//MACROS FOR LOGGING
#define GADGET_DEBUG3(err, message); \
	{std::string gdb ("[file %N, line %l] "); \
	gdb += message; \
	gdb += diagnostic_information(err); \
  ACE_DEBUG( (LM_DEBUG, \
	      ACE_TEXT(gdb.c_str() )));}

#endif  //GADGETRON_H
