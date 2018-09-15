#include "log.h"

#include "bart/bart_embed_api.h"

enum debug_levels { DP_ERROR, DP_WARN, DP_INFO, DP_DEBUG1, DP_DEBUG2, DP_DEBUG3, DP_DEBUG4, DP_TRACE, DP_ALL };

extern "C"
void vendor_log(int level,
        const char* /*func_name*/,
		const char* file,
		unsigned int line,
		const char* message)
{
     const auto fname(file + std::string(" (BART)"));
     const auto msg (message + std::string("\n"));
     
     if (-1 == debug_level) {
	  char* str = getenv("DEBUG_LEVEL");
	  debug_level = (NULL != str) ? atoi(str) : DP_INFO;
     }

     if (level <= debug_level) {
	  switch(level) {
	  case DP_ERROR:
	       Gadgetron::GadgetronLogger::instance()->log(Gadgetron::GADGETRON_LOG_LEVEL_ERROR, fname.c_str(), line, msg.c_str());
	       break;
	  case DP_WARN:
	       Gadgetron::GadgetronLogger::instance()->log(Gadgetron::GADGETRON_LOG_LEVEL_WARNING, fname.c_str(), line, msg.c_str());
	       break;
	  case DP_INFO:
	       Gadgetron::GadgetronLogger::instance()->log(Gadgetron::GADGETRON_LOG_LEVEL_INFO, fname.c_str(), line, msg.c_str());
	       break;
	  default:
	       Gadgetron::GadgetronLogger::instance()->log(Gadgetron::GADGETRON_LOG_LEVEL_DEBUG, fname.c_str(), line, msg.c_str());
	  }
     }
}
