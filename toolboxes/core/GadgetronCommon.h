#ifndef GADGETRONCOMMON_H
#define GADGETRONCOMMON_H

#ifndef _WIN32

#define GCC_VERSION (__GNUC__ * 10000           \
                     + __GNUC_MINOR__ * 1000    \
                     + __GNUC_PATCHLEVEL__)

#if GCC_VERSION < 42000
#pragma message ("GCC version is older than 4.2.0")
#define GCC_OLD_FLAG 1
#endif

#else

#endif // _WIN32

//MACROS FOR LOGGING
#define GADGET_MSG(message) { std::cout << message << std::endl; }
#define GADGET_ERROR_MSG(message) { std::cout << " (" << __FILE__ << ", " << __LINE__ << ") -> error happend: " << message << std::endl; }
#define GADGET_WARN_MSG(message) { std::cout << " (" << __FILE__ << ", " << __LINE__ << ") -> warning released: " << message << std::endl; }

#define GADGET_CONDITION_MSG(con, message) { if ( con ) GADGET_MSG(message) }
#define GADGET_CONDITION_WARN_MSG(con, message) { if ( con ) GADGET_WARN_MSG(message) }

#define GADGET_THROW(msg) { GADGET_ERROR_MSG(msg); BOOST_THROW_EXCEPTION( runtime_error(msg)); }
#define GADGET_CHECK_THROW(con) { if ( !(con) ) { GADGET_ERROR_MSG(#con); BOOST_THROW_EXCEPTION( runtime_error(#con)); } }

#define GADGET_CHECK_RETURN(con, value) { if ( ! (con) ) { GADGET_ERROR_MSG("Returning '" << value << "' due to failed check: '" << #con << "'"); return (value); } }
#define GADGET_CHECK_RETURN_FALSE(con) { if ( ! (con) ) { GADGET_ERROR_MSG("Returning false due to failed check: '" << #con << "'"); return false; } }

#ifdef GADGET_DEBUG_MODE
#define GADGET_DEBUG_CHECK_THROW(con) GADGET_CHECK_THROW(con)
#define GADGET_DEBUG_CHECK_RETURN(con, value) GADGET_CHECK_RETURN(con, value)
#define GADGET_DEBUG_CHECK_RETURN_FALSE(con) GADGET_CHECK_RETURN_FALSE(con)
#else
#define GADGET_DEBUG_CHECK_THROW(con)
#define GADGET_DEBUG_CHECK_RETURN(con, value)
#define GADGET_DEBUG_CHECK_RETURN_FALSE(con)
#endif // GADGET_DEBUG_MODE

// MACROS FOR TIMING
#define GADGET_START_TIMING(timer, oper) { timer.start(#oper); } 
#define GADGET_STOP_TIMING(timer) { timer.stop(); }

#define GADGET_START_TIMING_CONDITION(timer, oper, con) { if ( con ) { timer.start(#oper); } } 
#define GADGET_STOP_TIMING_CONDITION(timer, con) { if ( con ) { timer.stop(); } }

// MACROS FOR PRINTING
#define GADGET_OSTREAM_PRINT(os, content) { os << #content << " is " << content << std::endl; }

// MACROS FOR UTILITY
#define GT_MIN(a,b)    (((a)<(b))?(a):(b))
#define GT_MAX(a,b)    (((a)>(b))?(a):(b))
#define GT_ABS(a)      (((a)>=0)?(a):(-(a)))
#define GT_SGN(a)      (((a)>=0)?(1):(-1))
#define GT_PI          3.141592653589793238462
#define GT_IMAGING_GEOMETRY_DELTA 0.001

#endif  //GADGETRONCOMMON_H
