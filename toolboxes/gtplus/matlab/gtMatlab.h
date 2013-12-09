/********************************************************************
    created:    2013/10/03
    created:    3:10:2013   16:15
    author:     Hui Xue

    purpose:    Header to enable matlab print out info
*********************************************************************/

#pragma once 

#include <strstream>

#ifdef GADGET_MSG
    #undef GADGET_MSG
#endif // GADGET_MSG

#ifdef GADGET_ERROR_MSG
    #undef GADGET_ERROR_MSG
#endif // GADGET_ERROR_MSG

#ifdef GADGET_WARN_MSG
    #undef GADGET_WARN_MSG
#endif // GADGET_WARN_MSG

#ifdef _DEBUG
    #define GADGET_MSG(message) { std::ostrstream outs; outs << " (" << __FILE__ << ", " << __LINE__ << "): " << message << std::endl << '\0'; mexPrintf("%s", outs.str()); }
#else
    #define GADGET_MSG(message) { std::ostrstream outs; outs << message << std::endl << '\0'; mexPrintf("%s", outs.str()); }
#endif // _DEBUG

#ifdef _DEBUG
    #define GADGET_WARN_MSG(message) { std::ostrstream outs; outs << " (" << __FILE__ << ", " << __LINE__ << "): " << message << std::endl << '\0'; mexWarnMsgTxt(outs.str()); }
#else
    #define GADGET_WARN_MSG(message) { std::ostrstream outs; outs << message << std::endl << '\0'; mexWarnMsgTxt(outs.str()); }
#endif // _DEBUG

#define GADGET_ERROR_MSG(message) GADGET_MSG(message) 
