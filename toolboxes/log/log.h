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
#define GADGETRON_LOG_PRINT_FILELOC (1<<8)

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
#define GADGET_MSG_DEPRECATED(message) \
  {				       \
     std::stringstream str;	       \
     str << message << std::endl;      \
     GDEBUG(str.str().c_str());	       \
  }

#endif //GADGETRON_LOG_H
