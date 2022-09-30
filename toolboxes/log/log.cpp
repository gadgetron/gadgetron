#include "log.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string>
#include <time.h>
#include <cstring>
#include <chrono>


namespace Gadgetron
{
  GadgetronLogger* GadgetronLogger::instance()
  {
    if (!instance_) instance_ = new GadgetronLogger();
    return instance_;
  }

  GadgetronLogger* GadgetronLogger::instance_ = NULL;

  GadgetronLogger::GadgetronLogger()
    : level_mask_(GADGETRON_LOG_LEVEL_MAX,false)
    , print_mask_(GADGETRON_LOG_PRINT_MAX, false)
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

    // Redirect stderr to a log file
    char *log_file = getenv(GADGETRON_LOG_FILE_ENVIRONMENT);
    if (log_file != NULL) {
       fflush(stderr);
       FILE *newStdOut = freopen(log_file, "a", stderr);
       if (newStdOut == NULL) {
         fprintf(stderr, "Unable to redirect stderr to %s\n", log_file);
         fflush(stderr);
       }
    }
  }


  void GadgetronLogger::log(GadgetronLogLevel LEVEL, const char* filename, int lineno, const char* cformatting, ...)
  {
    //Check if we should log this message
    if (!isLevelEnabled(LEVEL)) return;

    std::unique_lock<std::mutex> lock(m);
    const char* fmt = cformatting;
    std::string fmt_str;
    bool append_cformatting_needed = false; //Will be set to true if we add any additional labels

    if (isOutputOptionEnabled(GADGETRON_LOG_PRINT_DATETIME)) {
      time_t rawtime;
      struct tm * timeinfo;

      auto curtime= std::chrono::system_clock::now();
      rawtime = std::chrono::system_clock::to_time_t(curtime);
      timeinfo = localtime ( &rawtime );

      auto duration = curtime.time_since_epoch();
      int micros = std::chrono::duration_cast<std::chrono::microseconds>(duration).count() % 1000000;

      //Time the format MM-DD HH:MM:SS.uuu
      char timestr[66];sprintf(timestr, "%02d-%02d %02d:%02d:%02d.%03d ",
			       timeinfo->tm_mon+1, timeinfo->tm_mday,
			       timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec, micros/1000);

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
	  fmt_str += std::string("[") + std::string(filename);
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
    vfprintf(stderr, fmt, args);
    va_end (args);
    fflush(stderr);
  }

  void GadgetronLogger::enableLogLevel(GadgetronLogLevel LEVEL)
  {
    if (LEVEL < level_mask_.size()) {
      level_mask_[LEVEL] = true;
    }
  }

  void GadgetronLogger::disableLogLevel(GadgetronLogLevel LEVEL)
  {
    if (LEVEL < level_mask_.size()) {
      level_mask_[LEVEL] = false;
    }
  }

  bool GadgetronLogger::isLevelEnabled(GadgetronLogLevel LEVEL)
  {
    if (LEVEL >= level_mask_.size()) return false;
    return level_mask_[LEVEL];
  }

  void GadgetronLogger::enableAllLogLevels()
  {
    level_mask_.assign(GADGETRON_LOG_LEVEL_MAX,true);
  }

  void GadgetronLogger::disableAllLogLevels()
  {
    level_mask_.assign(GADGETRON_LOG_LEVEL_MAX,false);
  }

  void GadgetronLogger::enableOutputOption(GadgetronLogOutput OUTPUT)
  {
    if (OUTPUT < print_mask_.size()) {
      print_mask_[OUTPUT] = true;
    }
  }

  void GadgetronLogger::disableOutputOption(GadgetronLogOutput OUTPUT)
  {
    if (OUTPUT < print_mask_.size()) {
      print_mask_[OUTPUT] = false;
    }
  }

  bool GadgetronLogger::isOutputOptionEnabled(GadgetronLogOutput OUTPUT)
  {
    if (OUTPUT < print_mask_.size()) {
      return print_mask_[OUTPUT];
    }
    return false;
  }

  void GadgetronLogger::enableAllOutputOptions()
  {
    print_mask_.assign(GADGETRON_LOG_PRINT_MAX, true);
  }

  void GadgetronLogger::disableAllOutputOptions()
  {
    print_mask_.assign(GADGETRON_LOG_PRINT_MAX, false);
  }
}
