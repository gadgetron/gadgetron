#include "log.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string>
#include <time.h>
#include <cstring>

namespace Gadgetron
{
  GadgetronLogger* GadgetronLogger::instance()
  {
    if (!instance_) instance_ = new GadgetronLogger();
    return instance_;
  }
  
  GadgetronLogger* GadgetronLogger::instance_ = NULL;
  
  GadgetronLogger::GadgetronLogger()
    : level_mask_(0)
    , print_mask_(0)
  {
    char* log_mask = getenv(GADGETRON_LOG_MASK_ENVIRONMENT);
    if ( log_mask != NULL) {
      
      std::string log_mask_str(log_mask);

      //Which log levels are enabled
      if (log_mask_str.find("ALL") != std::string::npos) {
	enableAllOutputOptions();
	enableAllLogLevels();
	return;
      }
      
      if (log_mask_str.find("LEVEL_DEBUG") != std::string::npos) 
	enableLogLevel(GADGETRON_LOG_LEVEL_DEBUG);

      if (log_mask_str.find("LEVEL_INFO") != std::string::npos) 
	enableLogLevel(GADGETRON_LOG_LEVEL_INFO);

      if (log_mask_str.find("LEVEL_WARNING") != std::string::npos) 
	enableLogLevel(GADGETRON_LOG_LEVEL_WARNING);

      if (log_mask_str.find("LEVEL_ERROR") != std::string::npos) 
	enableLogLevel(GADGETRON_LOG_LEVEL_ERROR);

      if (log_mask_str.find("PRINT_FILELOC") != std::string::npos) 
	enableOutputOption(GADGETRON_LOG_PRINT_FILELOC);

      if (log_mask_str.find("PRINT_LEVEL") != std::string::npos) 
	enableOutputOption(GADGETRON_LOG_PRINT_LEVEL);
      
      if (log_mask_str.find("PRINT_DATETIME") != std::string::npos) 
	enableOutputOption(GADGETRON_LOG_PRINT_DATETIME);
    } else {
      enableLogLevel(GADGETRON_LOG_LEVEL_DEBUG);
      enableLogLevel(GADGETRON_LOG_LEVEL_INFO);
      enableLogLevel(GADGETRON_LOG_LEVEL_WARNING);
      enableLogLevel(GADGETRON_LOG_LEVEL_ERROR);
      enableOutputOption(GADGETRON_LOG_PRINT_FILELOC);
      enableOutputOption(GADGETRON_LOG_PRINT_LEVEL);
      enableOutputOption(GADGETRON_LOG_PRINT_DATETIME);
    }
  }


  void GadgetronLogger::log(GadgetronLogLevel LEVEL, const char* filename, int lineno, const char* cformatting, ...)
  {
    //Check if we should log this message
    if (!isLevelEnabled(LEVEL)) return;

    const char* fmt = cformatting;
    std::string fmt_str;
    bool append_cformatting_needed = false; //Will be set to true if we add any additional labels

    if (isOutputOptionEnabled(GADGETRON_LOG_PRINT_DATETIME)) {
      time_t rawtime;
      struct tm * timeinfo;

      time ( &rawtime );
      timeinfo = localtime ( &rawtime );
      
      //Time the format YYYY-MM-DD HH:MM:SS
      char timestr[22];sprintf(timestr, "%d-%02d-%02d %02d:%02d:%02d ",
			       timeinfo->tm_year+1900, timeinfo->tm_mon+1, timeinfo->tm_mday,
			       timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec);

      fmt_str += std::string(timestr);
      append_cformatting_needed = true;
    }

    if (isOutputOptionEnabled(GADGETRON_LOG_PRINT_LEVEL)) {
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
      append_cformatting_needed = true;
    }

    if (isOutputOptionEnabled(GADGETRON_LOG_PRINT_FILELOC)) {
      if (!isOutputOptionEnabled(GADGETRON_LOG_PRINT_FOLDER)) {
	const char* base_start = strrchr(filename,'/');
	if (!base_start) {
	  base_start = strrchr(filename,'\\'); //Maybe using backslashes
	}
	if (base_start) {
	  base_start++;
	  fmt_str += std::string("[") + std::string(base_start);
	} else {
	  std::string("[") + std::string(filename);
	}
      } else {
	fmt_str += std::string("[") + std::string(filename);
      }
      char linenostr[8];sprintf(linenostr, "%d", lineno);
      fmt_str += std::string(":") + std::string(linenostr);
      fmt_str += std::string("] ");
      append_cformatting_needed = true;
    }

    if (append_cformatting_needed) {
      fmt_str += std::string(cformatting); 
      fmt = fmt_str.c_str();      
    }

    va_list args;
    va_start (args, cformatting);
    vprintf(fmt, args);
    va_end (args);
    fflush(stdout);
  }

  void GadgetronLogger::enableLogLevel(GadgetronLogLevel LEVEL)
  {
    if (LEVEL < GADGETRON_LOG_LEVEL_MAX) {
      level_mask_ |= (1<<LEVEL);
    }
  }

  void GadgetronLogger::disableLogLevel(GadgetronLogLevel LEVEL)
  {
    if (LEVEL < GADGETRON_LOG_LEVEL_MAX) {
      level_mask_ &= ~(1<<LEVEL);
    }
  }
  
  bool GadgetronLogger::isLevelEnabled(GadgetronLogLevel LEVEL)
  {
    if (LEVEL >= GADGETRON_LOG_LEVEL_MAX) return false;

    return ((1<<LEVEL) & level_mask_)>0;
  }
  
  void GadgetronLogger::enableAllLogLevels()
  {
    uint64_t null_mask = 0;
    level_mask_ |= ~null_mask;
  }

  void GadgetronLogger::disableAllLogLevels()
  {
    level_mask_ = 0;
  }

  void GadgetronLogger::enableOutputOption(GadgetronLogOutput OUTPUT) 
  {
    if (OUTPUT < GADGETRON_LOG_PRINT_MAX) {
      print_mask_ |= (1<<OUTPUT);
    }    
  }

  void GadgetronLogger::disableOutputOption(GadgetronLogOutput OUTPUT) 
  {
    if (OUTPUT < GADGETRON_LOG_PRINT_MAX) {
      print_mask_ &= ~(1<<OUTPUT);
    }    
  }
 
  bool GadgetronLogger::isOutputOptionEnabled(GadgetronLogOutput OUTPUT)
  {
    if (OUTPUT >= GADGETRON_LOG_PRINT_MAX) return false;
    return ((1<<OUTPUT) & print_mask_)>0;
  }
  
  void GadgetronLogger::enableAllOutputOptions() 
  {
    uint64_t null_mask = 0;
    print_mask_ |= ~null_mask;
  }

  void GadgetronLogger::disableAllOutputOptions() 
  {
    print_mask_ = 0;
  }

}
