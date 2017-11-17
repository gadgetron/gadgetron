/**
\file       Matlab_info.h
\brief      Header file to store matlab mex related information
\author     Hui Xue
*/

#pragma once

inline void printAuthorInfo(std::stringstream& outs)
{
    using namespace std;
    outs << "---------------------------------------------------------------------" << endl;
    outs << "This mex software is made by: " << endl;
    outs << endl;
    outs << "\t\tHui Xue " << endl;
    outs << "Email: hui.xue@nih.gov" << endl;
    outs << "National Heart, Lung and Blood Institute" << endl;
    outs << "National Institutes of Health" << endl;
    outs << "---------------------------------------------------------------------" << endl;
}
