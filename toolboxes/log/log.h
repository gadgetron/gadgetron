#ifndef GADGETRON_LOG_H
#define GADGETRON_LOG_H

#include "log_export.h"

#include <stdint.h>

#include <sstream> //For deprecated macros

//Log levels
#define GADGETRON_LOG_LEVEL_DEBUG   (1<<0)
#define GADGETRON_LOG_LEVEL_INFO    (1<<1)
#define GADGETRON_LOG_LEVEL_WARNING (1<<2)
#define GADGETRON_LOG_LEVEL_ERROR   (1<<3)

//Log output options
#define GADGETRON_LOG_PRINT_FILELOC  (1<<8)
#define GADGETRON_LOG_PRINT_LEVEL    (1<<9)
#define GADGETRON_LOG_PRINT_DATETIME (1<<10)

#define GADGETRON_LOG_MASK_ENVIRONMENT "GADGETRON_LOG_MASK"

namespace Gadgetron
{
  class EXPORTGADGETRONLOG GadgetronLogger
  {
  public:
    static GadgetronLogger* instance();
    void Log(uint16_t LEVEL, const char* filename, int lineno, const char* cformatting, ...);

  protected:
    GadgetronLogger();
    static GadgetronLogger* instance_;
    uint16_t log_mask_;
  };
}

#define GDEBUG(...) Gadgetron::GadgetronLogger::instance()->Log(GADGETRON_LOG_LEVEL_DEBUG,   __FILE__, __LINE__, __VA_ARGS__)
#define GINFO(...)  Gadgetron::GadgetronLogger::instance()->Log(GADGETRON_LOG_LEVEL_INFO,    __FILE__, __LINE__, __VA_ARGS__)
#define GWARN(...)  Gadgetron::GadgetronLogger::instance()->Log(GADGETRON_LOG_LEVEL_WARNING, __FILE__, __LINE__, __VA_ARGS__)
#define GERROR(...) Gadgetron::GadgetronLogger::instance()->Log(GADGETRON_LOG_LEVEL_ERROR,   __FILE__, __LINE__, __VA_ARGS__)

#define GEXCEPTION(err, message);	  \
  {					  \
    std::string gdb(message);		  \
    gdb += std::string(" --> ");	  \
    gdb += err.what();			  \
    GDEBUG(gdb.c_str());		  \
 }

//Deprecated logging macros kept here for backwards compatibility
#define GDEBUG_STREAM(message)				\
  {							\
    std::stringstream gadget_msg_dep_str;		\
    gadget_msg_dep_str  << message << std::endl;	\
    GDEBUG(gadget_msg_dep_str.str().c_str());		\
  }

#define GERROR_STREAM(message)					\
  {								\
    std::stringstream gadget_msg_dep_str;			\
    gadget_msg_dep_str  << message << std::endl;		\
    GERROR(gadget_msg_dep_str.str().c_str());			\
  }
     
#define GWARN_STREAM(message)					\
  {								\
    std::stringstream gadget_msg_dep_str;			\
    gadget_msg_dep_str  << message << std::endl;		\
    GERROR(gadget_msg_dep_str.str().c_str());			\
  }

//Older debugging macros
//TODO: Review and check that they are up to date
#define GDEBUG_CONDITION_STREAM(con, message) { if ( con ) GDEBUG_STREAM(message) }
#define GWARN_CONDITION_STREAM(con, message) { if ( con ) GWARN_STREAM(message) }
     
#define GADGET_THROW(msg) { GERROR_STREAM(msg); throw std::runtime_error(msg); }
#define GADGET_CHECK_THROW(con) { if ( !(con) ) { GERROR_STREAM(#con); throw std::runtime_error(#con); } }

#define GADGET_CATCH_THROW(con) { try { con; } catch(...) { GERROR_STREAM(#con); throw std::runtime_error(#con); } }

#define GADGET_CHECK_RETURN(con, value) { if ( ! (con) ) { GERROR_STREAM("Returning '" << value << "' due to failed check: '" << #con << "'"); return (value); } }
#define GADGET_CHECK_RETURN_FALSE(con) { if ( ! (con) ) { GERROR_STREAM("Returning false due to failed check: '" << #con << "'"); return false; } }

#define GADGET_CHECK_EXCEPTION_RETURN(con, value) { try { con; } catch(...) { GERROR_STREAM("Returning '" << value << "' due to failed check: '" << #con << "'"); return (value); } }
#define GADGET_CHECK_EXCEPTION_RETURN_FALSE(con) { try { con; } catch(...) { GERROR_STREAM("Returning false due to failed check: '" << #con << "'"); return false; } }

#ifdef GADGET_DEBUG_MODE
#define GADGET_DEBUG_CHECK_THROW(con) GADGET_CHECK_THROW(con)
#define GADGET_DEBUG_CHECK_RETURN(con, value) GADGET_CHECK_RETURN(con, value)
#define GADGET_DEBUG_CHECK_RETURN_FALSE(con) GADGET_CHECK_RETURN_FALSE(con)
#else
#define GADGET_DEBUG_CHECK_THROW(con)
#define GADGET_DEBUG_CHECK_RETURN(con, value)
#define GADGET_DEBUG_CHECK_RETURN_FALSE(con)
#endif // GADGET_DEBUG_MODE


#endif //GADGETRON_LOG_H
