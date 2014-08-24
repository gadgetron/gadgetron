/********************************************************************
    created:    2013/10/03
    created:    3:10:2013   16:15
    author:     Hui Xue

    purpose:    Header to enable matlab print out info
*********************************************************************/

#pragma once 

#include <sstream>
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

#ifdef GADGET_CHECK_RETURN_FALSE
    #undef GADGET_CHECK_RETURN_FALSE
#endif // GADGET_CHECK_RETURN_FALSE
#define GADGET_CHECK_RETURN_FALSE(con) { if ( ! (con) ) { return false; } }

#ifdef GADGET_DEBUG_MODE
#define GADGET_DEBUG_CHECK_THROW(con) GADGET_CHECK_THROW(con)
#define GADGET_DEBUG_CHECK_RETURN(con, value) GADGET_CHECK_RETURN(con, value)
#define GADGET_DEBUG_CHECK_RETURN_FALSE(con) GADGET_CHECK_RETURN_FALSE(con)
#else
#define GADGET_DEBUG_CHECK_THROW(con)
#define GADGET_DEBUG_CHECK_RETURN(con, value)
#define GADGET_DEBUG_CHECK_RETURN_FALSE(con)
#endif // GADGET_DEBUG_MODE

template <typename ObjType> void matlab_printInfo(const ObjType& obj)
{
    std::ostrstream outs;
    obj.print(outs);
    outs << std::ends;
    std::string msg(outs.str());
    GADGET_MSG(msg.c_str());
}

inline void printAuthorInfo(std::stringstream& outs)
{
    using namespace std;
    outs << "---------------------------------------------------------------------" << endl;
    outs << "This software is made by: " << endl;
    outs << endl;
    outs << "\t\tHui Xue " << endl;
    outs << "Magnetic Resonance Technology Program" << endl;
    outs << "National Heart, Lung and Blood Institute" << endl;
    outs << "National Institutes of Health" << endl;
    outs << "Email: hui.xue@nih.gov" << endl;
    outs << endl;
    outs << "\t\tPeter Kellman " << endl;
    outs << "Medical Signal and Image Processing Program" << endl;
    outs << "National Heart, Lung and Blood Institute" << endl;
    outs << "National Institutes of Health" << endl;
    outs << "Email: kellmanp@nhlbi.nih.gov" << endl;
    outs << endl;
    outs << "\t\tMichael Hansen " << endl;
    outs << "Medical Signal and Image Processing Program" << endl;
    outs << "National Heart, Lung and Blood Institute" << endl;
    outs << "National Institutes of Health" << endl;
    outs << "Email: michael.hansen@nih.gov" << endl;
    outs << "---------------------------------------------------------------------" << endl;
}
