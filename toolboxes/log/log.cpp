#include "log.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string>
#include <time.h>

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
      log_mask_ |= GADGETRON_LOG_PRINT_LEVEL;
      log_mask_ |= GADGETRON_LOG_PRINT_DATETIME;
    }
  }


  void GadgetronLogger::Log(uint16_t LEVEL, const char* filename, int lineno, const char* cformatting, ...)
  {
    //Check if we should log this message
    if (!(LEVEL & log_mask_)) return;

    const char* fmt = cformatting;
    std::string fmt_str;

    if (log_mask_ & GADGETRON_LOG_PRINT_DATETIME) {
      time_t rawtime;
      struct tm * timeinfo;

      time ( &rawtime );
      timeinfo = localtime ( &rawtime );
      
      //Time the format YYYY-MM-DD HH:MM:SS
      char timestr[22];sprintf(timestr, "%d-%02d-%02d %02d:%02d:%02d ",
			       timeinfo->tm_year+1900, timeinfo->tm_mon+1, timeinfo->tm_mday,
			       timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec);

      fmt_str += std::string(timestr);
      fmt = fmt_str.c_str();
    }

    if (log_mask_ & GADGETRON_LOG_PRINT_LEVEL) {
      switch (LEVEL) {
      case GADGETRON_LOG_LEVEL_DEBUG:
	fmt_str += "DEBUG ";
	break;
      case GADGETRON_LOG_LEVEL_INFO:
	fmt_str += "INFO ";
	break;
      case GADGETRON_LOG_LEVEL_WARNING:
	fmt_str += "WARNING ";
	break;
      case GADGETRON_LOG_LEVEL_ERROR:
	fmt_str += "ERROR ";
	break;
      default:
	;
      }
      fmt = fmt_str.c_str();
    }

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
