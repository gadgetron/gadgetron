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
#define GADGET_DEBUG_EXCEPTION(err, message); \
	{std::string gdb ("[file %N, line %l] "); \
	gdb += message; \
	gdb += err.what(); \
  ACE_DEBUG( (LM_DEBUG, \
	      ACE_TEXT(gdb.c_str() )));}

//MACROS FOR LOGGING
#define GADGET_MSG(message) { std::cout << message << std::endl; }
#define GADGET_ERROR_MSG(message) { std::cout << " (" << __FILE__ << ", " << __LINE__ << ") -> error happend: " << message << std::endl; }
#define GADGET_WARN_MSG(message) { std::cout << " (" << __FILE__ << ", " << __LINE__ << ") -> warning released: " << message << std::endl; }

#define GADGET_THROW(msg) { GADGET_ERROR_MSG(msg); BOOST_THROW_EXCEPTION( runtime_error(msg)); }
#define GADGET_CHECK_THROW(con) { if ( !(con) ) { GADGET_ERROR_MSG(#con); BOOST_THROW_EXCEPTION( runtime_error(#con)); } }

#define GADGET_CHECK_RETURN(con, value) { if ( ! (con) ) { GADGET_ERROR_MSG("Returning '" << value << "' due to failed check: '" << #con << "'"); return (value); } }
#define GADGET_CHECK_RETURN_FALSE(con) { if ( ! (con) ) { GADGET_ERROR_MSG("Returning false due to failed check: '" << #con << "'"); return false; } }

#define GADGET_CHECK_PERFORM(con, action) { if ( con ) { action; } }

#define GADGET_EXPORT_ARRAY(debugFolder, exporter, a, filename) { if ( !debugFolder.empty() ) { exporter.exportArray(a, debugFolder+filename); } }
#define GADGET_EXPORT_ARRAY_COMPLEX(debugFolder, exporter, a, filename) { if ( !debugFolder.empty() ) { exporter.exportArrayComplex(a, debugFolder+filename); } }

#ifdef GADGET_DEBUG_MODE
    #define GADGET_DEBUG_CHECK_THROW(con) GADGET_CHECK_THROW(con)
    #define GADGET_DEBUG_CHECK_RETURN(con, value) GADGET_CHECK_RETURN(con, value)
    #define GADGET_DEBUG_CHECK_RETURN_FALSE(con) GADGET_CHECK_RETURN_FALSE(con)
#else
    #define GADGET_DEBUG_CHECK_THROW(con)
    #define GADGET_DEBUG_CHECK_RETURN(con, value)
    #define GADGET_DEBUG_CHECK_RETURN_FALSE(con)
#endif // GADGET_DEBUG_MODE

// MACROS FOR UTILITY
#define GT_MIN(a,b)    (((a)<(b))?(a):(b))
#define GT_MAX(a,b)    (((a)>(b))?(a):(b))
#define GT_ABS(a)      (((a)>=0)?(a):(-(a)))
#define GT_SGN(a)      (((a)>=0)?(1):(-1))
#define GT_PI          3.141592653589793238462

#endif  //GADGETRON_H
