#include "log.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string>

namespace Gadgetron
{
  GadgetronLogger* GadgetronLogger::instance()
  {
    if (!instance_) instance_ = new GadgetronLogger();
    return instance_;
  }
  
  GadgetronLogger* GadgetronLogger::instance_ = NULL;
  
  GadgetronLogger::GadgetronLogger()
    : log_mask_(0)
  {
    char* log_mask = getenv(GADGETRON_LOG_MASK_ENVIRONMENT);
    if ( log_mask != NULL) {
      log_mask_ = static_cast<uint16_t>(atoi(log_mask));
    } else {
      log_mask_ |= GADGETRON_LOG_LEVEL_DEBUG;
      log_mask_ |= GADGETRON_LOG_LEVEL_INFO;
      log_mask_ |= GADGETRON_LOG_LEVEL_WARNING;
      log_mask_ |= GADGETRON_LOG_LEVEL_ERROR;
      log_mask_ |= GADGETRON_LOG_PRINT_FILELOC;
    }
  }


  void GadgetronLogger::Log(uint16_t LEVEL, const char* filename, int lineno, const char* cformatting, ...)
  {
    //Check if we should log this message
    if (!(LEVEL & log_mask_)) return;

    const char* fmt = cformatting;

    std::string fmt_str;
    if (log_mask_ & GADGETRON_LOG_PRINT_FILELOC) {
      char linenostr[8];sprintf(linenostr, "%d", lineno);
      fmt_str += std::string("[") + std::string(filename);
      fmt_str += std::string(":") + std::string(linenostr);
      fmt_str += std::string("] ");
      fmt_str += std::string(cformatting); 
      fmt = fmt_str.c_str();
    }

    va_list args;
    va_start (args, cformatting);
    vprintf(fmt, args);
    va_end (args);
    fflush(stdout);
  }
}
