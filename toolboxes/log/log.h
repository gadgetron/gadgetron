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

#define GEXCEPTION(err, message); \
  {					  \
    std::string gdb(message);		  \
    gdb += std::string(" --> ");	  \
    gdb += err.what();			  \
    GDEBUG(gdb.c_str());		  \
 }

//Deprecated logging macros kept here for backwards compatibility
#define GADGET_MSG_DEPRECATED(message)			\
  {							\
    std::stringstream gadget_msg_dep_str;		\
    gadget_msg_dep_str  << message << std::endl;	\
    GDEBUG(gadget_msg_dep_str.str().c_str());		\
  }

//Old macros defined but Hui Xue, but should no longer be used
//MACROS FOR LOGGING
//#define GADGET_MSG_DEPRECATED(message) { std::cout << message << std::endl; }
#define GADGET_ERROR_MSG(message) { std::cout << " (" << __FILE__ << ", " << __LINE__ << ") -> error happend: " << message << std::endl; }
#define GADGET_WARN_MSG(message) { std::cout << " (" << __FILE__ << ", " << __LINE__ << ") -> warning released: " << message << std::endl; }

#define GADGET_CONDITION_MSG(con, message) { if ( con ) GADGET_MSG_DEPRECATED(message) }
#define GADGET_CONDITION_WARN_MSG(con, message) { if ( con ) GADGET_WARN_MSG(message) }

#define GADGET_THROW(msg) { GADGET_ERROR_MSG(msg); throw std::runtime_error(msg); }
#define GADGET_CHECK_THROW(con) { if ( !(con) ) { GADGET_ERROR_MSG(#con); throw std::runtime_error(#con); } }

#define GADGET_CATCH_THROW(con) { try { con; } catch(...) { GADGET_ERROR_MSG(#con); throw std::runtime_error(#con); } }

#define GADGET_CHECK_RETURN(con, value) { if ( ! (con) ) { GADGET_ERROR_MSG("Returning '" << value << "' due to failed check: '" << #con << "'"); return (value); } }
#define GADGET_CHECK_RETURN_FALSE(con) { if ( ! (con) ) { GADGET_ERROR_MSG("Returning false due to failed check: '" << #con << "'"); return false; } }

#define GADGET_CHECK_EXCEPTION_RETURN(con, value) { try { con; } catch(...) { GADGET_ERROR_MSG("Returning '" << value << "' due to failed check: '" << #con << "'"); return (value); } }
#define GADGET_CHECK_EXCEPTION_RETURN_FALSE(con) { try { con; } catch(...) { GADGET_ERROR_MSG("Returning false due to failed check: '" << #con << "'"); return false; } }

#ifdef GADGET_DEBUG_MODE
#define GADGET_DEBUG_CHECK_THROW(con) GADGET_CHECK_THROW(con)
#define GADGET_DEBUG_CHECK_RETURN(con, value) GADGET_CHECK_RETURN(con, value)
#define GADGET_DEBUG_CHECK_RETURN_FALSE(con) GADGET_CHECK_RETURN_FALSE(con)
#else
#define GADGET_DEBUG_CHECK_THROW(con)
#define GADGET_DEBUG_CHECK_RETURN(con, value)
#define GADGET_DEBUG_CHECK_RETURN_FALSE(con)
#endif // GADGET_DEBUG_MODE

// Macros FOR TIMING
#define GADGET_START_TIMING(timer, oper) { timer.start(#oper); } 
#define GADGET_STOP_TIMING(timer) { timer.stop(); }

#define GADGET_START_TIMING_CONDITION(timer, oper, con) { if ( con ) { timer.start(#oper); } } 
#define GADGET_STOP_TIMING_CONDITION(timer, con) { if ( con ) { timer.stop(); } }

// MACROS FOR PRINTING
#define GADGET_OSTREAM_PRINT(os, content) { os << #content << " is " << content << std::endl; }

#define GADGET_CHECK_PERFORM(con, action) { if ( con ) { action; } }

// MACROS for EXPORTING
#define GADGET_EXPORT_ARRAY(debugFolder, exporter, a, filename) { if ( !debugFolder.empty() ) { exporter.exportArray(a, debugFolder+filename); } }
#define GADGET_EXPORT_ARRAY_COMPLEX(debugFolder, exporter, a, filename) { if ( !debugFolder.empty() ) { exporter.exportArrayComplex(a, debugFolder+filename); } }
#define GADGET_EXPORT_ARRAY_COMPLEX_REAL_IMAG(debugFolder, exporter, a, filename) { if ( !debugFolder.empty() ) { exporter.exportArrayComplexRealImag(a, debugFolder+filename); } }

#define GADGET_EXPORT_IMAGE(debugFolder, exporter, a, filename) { if ( !debugFolder.empty() ) { exporter.exportImage(a, debugFolder+filename); } }
#define GADGET_EXPORT_IMAGE_COMPLEX(debugFolder, exporter, a, filename) { if ( !debugFolder.empty() ) { exporter.exportImageComplex(a, debugFolder+filename); } }


#endif //GADGETRON_LOG_H
