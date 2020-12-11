#ifndef GADGETRON_LOG_H
#define GADGETRON_LOG_H

#include "log_export.h"

#include <vector> //For mask fields

#include <sstream> //For deprecated macros

#define GADGETRON_LOG_MASK_ENVIRONMENT "GADGETRON_LOG_MASK"
#define GADGETRON_LOG_FILE_ENVIRONMENT "GADGETRON_LOG_FILE"

namespace Gadgetron
{
  /**
     Gadgetron log levels
   */
  enum GadgetronLogLevel
  {
    GADGETRON_LOG_LEVEL_DEBUG = 0,     //!< Debug information
    GADGETRON_LOG_LEVEL_INFO,          //!< Regular application information
    GADGETRON_LOG_LEVEL_WARNING,       //!< Warnings about events that could lead to failues
    GADGETRON_LOG_LEVEL_ERROR,         //!< Errors after which application will be unable to continue
    GADGETRON_LOG_LEVEL_VERBOSE,       //!< Verbose information about algorithm parameters, etc. 
    GADGETRON_LOG_LEVEL_MAX            //!< All log levels must have values lower than this
  };

  /**
     Gadgetron output options. These options control what context information
     will be printed with the log statements.
   */
  enum GadgetronLogOutput
  {
    GADGETRON_LOG_PRINT_FILELOC = 0,  //!< Print filename and line in file
    GADGETRON_LOG_PRINT_FOLDER,       //!< Print the folder name too (full filename)
    GADGETRON_LOG_PRINT_LEVEL,        //!< Print Log Level
    GADGETRON_LOG_PRINT_DATETIME,     //!< Print date and time
    GADGETRON_LOG_PRINT_MAX           //!< All print options must have lower values than this
  };

  /**
     Main logging utility class for the Gadgetron and associated toolboxes. 

     This is a process wide singleton. 

     Logging/Debug messages should be done with the convenience macros:

     GDEBUG
     GINFO
     GWARN
     GERROR
     GVERBOSE

     These macros use a printf style syntax:

     GDEBUG("Here we are logging some values %f, %d\n", myFloat, myInt);

     The c++ std::cout is not recommended for logging as it does not add
     filename, log level, timing or other context information to the logging
     statements. Use of std::cout can also cause log lines from different threads
     to be interleaved. For people more comfortable with the std::cout style 
     syntax, we provide the macros:

     GDEBUG_STREAM
     GINFO_STREAM
     GWARN_STREAM
     GERROR_STREAM
     GVERBOSE_STREAM

     To use them:
     
     GDEBUG("Here we are logging some values " << myFloat << ", " << myInt << std::endl);

     It is possible to control which log levels are output using the @enableLogLevel, @disableLogLevel, 
     @enableOutputOption, and @disableOutputOption functions. 

     Log levels are defined in @GadgetronLogLevel
     Ouput options are defined in @GadgetronLogOutput
 
     The logger checks the environment variable GADGETRON_LOG_MASK. If it is set,
     it disables all log levels and outputs and only enables the ones in the mask. 
     It can be specified with (on unix system):

     export GADGETRON_LOG_MASK="LEVEL_INFO,LEVEL_DEBUG,PRINT_FILELOC,PRINT_DATETIME"
     
     Any (or no) seperator is allowed between the levels and ourput options.

   */
  class EXPORTGADGETRONLOG GadgetronLogger
  {
  public:
    ///Function for accessing the process wide singleton
    static GadgetronLogger* instance();

    ///Generic log function. Use the logging macros for easy access to this function
    void log(GadgetronLogLevel LEVEL, const char* filename, int lineno, const char* cformatting, ...);

    void enableLogLevel(GadgetronLogLevel LEVEL);
    void disableLogLevel(GadgetronLogLevel LEVEL);
    bool isLevelEnabled(GadgetronLogLevel LEVEL);
    void enableAllLogLevels();
    void disableAllLogLevels();

    void enableOutputOption(GadgetronLogOutput OUTPUT);
    void disableOutputOption(GadgetronLogOutput OUTPUT);
    bool isOutputOptionEnabled(GadgetronLogOutput OUTPUT);
    void enableAllOutputOptions();
    void disableAllOutputOptions();

  protected:
    GadgetronLogger();
    static GadgetronLogger* instance_;
    std::vector<bool> level_mask_;
    std::vector<bool> print_mask_;
  };
}

#define GDEBUG(...)   Gadgetron::GadgetronLogger::instance()->log(Gadgetron::GADGETRON_LOG_LEVEL_DEBUG,   __FILE__, __LINE__, __VA_ARGS__)
#define GINFO(...)    Gadgetron::GadgetronLogger::instance()->log(Gadgetron::GADGETRON_LOG_LEVEL_INFO,    __FILE__, __LINE__, __VA_ARGS__)
#define GWARN(...)    Gadgetron::GadgetronLogger::instance()->log(Gadgetron::GADGETRON_LOG_LEVEL_WARNING, __FILE__, __LINE__, __VA_ARGS__)
#define GERROR(...)   Gadgetron::GadgetronLogger::instance()->log(Gadgetron::GADGETRON_LOG_LEVEL_ERROR,   __FILE__, __LINE__, __VA_ARGS__)
#define GVERBOSE(...) Gadgetron::GadgetronLogger::instance()->log(Gadgetron::GADGETRON_LOG_LEVEL_VERBOSE,   __FILE__, __LINE__, __VA_ARGS__)

#define GEXCEPTION(err, message);	  \
  {					  \
    std::string gdb(message);		  \
    gdb += std::string(" --> ");	  \
    gdb += err.what();			  \
    GDEBUG(gdb.c_str());		  \
 }

//Stream syntax log level functions
#define GINFO_STREAM(message)				\
  {							\
    std::stringstream gadget_msg_dep_str;		\
    gadget_msg_dep_str  << message << std::endl;	\
    GINFO(gadget_msg_dep_str.str().c_str());		\
  }

#define GVERBOSE_STREAM(message)					\
  {								\
    std::stringstream gadget_msg_dep_str;			\
    gadget_msg_dep_str  << message << std::endl;		\
    GVERBOSE(gadget_msg_dep_str.str().c_str());			\
  }

#ifndef MATLAB_MEX_COMPILE

#define GDEBUG_STREAM(message)				\
{							\
    std::stringstream gadget_msg_dep_str;		\
    gadget_msg_dep_str  << message << std::endl;	\
    GDEBUG(gadget_msg_dep_str.str().c_str());		\
}

#define GWARN_STREAM(message)					\
  {								\
    std::stringstream gadget_msg_dep_str;			\
    gadget_msg_dep_str  << message << std::endl;		\
    GWARN(gadget_msg_dep_str.str().c_str());			\
  }

#define GERROR_STREAM(message)					\
  {								\
    std::stringstream gadget_msg_dep_str;			\
    gadget_msg_dep_str  << message << std::endl;		\
    GERROR(gadget_msg_dep_str.str().c_str());			\
  }

#else
    #pragma message ("Use matlab definition for GDEBUG stream ... ")

    #ifdef _DEBUG
        #define GDEBUG_STREAM(message) { std::ostrstream outs; outs << " (" << __FILE__ << ", " << __LINE__ << "): " << message << std::endl << '\0'; mexPrintf("%s", outs.str()); }
    #else
        #define GDEBUG_STREAM(message) { std::ostrstream outs; outs << message << std::endl << '\0'; mexPrintf("%s", outs.str()); }
    #endif // _DEBUG

    #ifdef _DEBUG
        #define GWARN_STREAM(message) { std::ostrstream outs; outs << " (" << __FILE__ << ", " << __LINE__ << "): " << message << std::endl << '\0'; mexWarnMsgTxt(outs.str()); }
    #else
        #define GWARN_STREAM(message) { std::ostrstream outs; outs << message << std::endl << '\0'; mexWarnMsgTxt(outs.str()); }
    #endif // _DEBUG

    #define GERROR_STREAM(message) GDEBUG_STREAM(message) 
#endif // MATLAB_MEX_COMPILE

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


//namespace boost
//{
//#ifdef BOOST_NO_EXCEPTIONS
//    inline void throw_exception(std::exception const & e)
//    {
//        throw e; // or whatever
//    };
//#endif
//
//}// namespace boost

#endif //GADGETRON_LOG_H
