
#include <matrix.h>
#include <mat.h>
#ifdef MATLAB_DLL_EXPORT_SYM
#define DLL_EXPORT_SYM extern "C" __declspec(dllexport)
#endif // MATLAB_DLL_EXPORT_SYM
#include <mex.h>

// Gadgetron includes
#include <strstream>
#include "GadgetronTimer.h"
#include "MatlabUtils.h"
#include "Matlab_info.h"

#include "cmr_motion_correction.h"

#define MEXPRINTF(name) mexPrintf(#name);

#define NIn 3
#define NOut 2

static void usage()
{
    using namespace std;
    std::stringstream outs;

    outs << "==============================================================================================" << endl;
    outs << "Usage: Matlab_gt_concatenate_deform_fields \n";
    outs << "Concatenate deformation fields" << endl;
    printAuthorInfo(outs);
    outs << NIn << " Input paras:" << endl;
    outs << '\t' << "dx                                         : RO*E1*N*, in float" << endl;
    outs << '\t' << "dy                                         : RO*E1*N*, in float" << endl;
    outs << '\t' << "key_frame                                  : key frame, [0 N-1]" << endl;

    outs << "2 Output para:" << endl;
    outs << '\t' << "dx_out      : RO*E1*N array" << endl;
    outs << '\t' << "dy_out      : RO*E1*N array" << endl;
    outs << "==============================================================================================" << endl;
    outs << '\t' << "Example usage: "<< endl;

    outs << '\t' << "[dx_out, dy_out] = Matlab_gt_concatenate_deform_fields(single(dx), single(dy), key_frame); "<< endl;

    outs << "==============================================================================================" << endl;
    outs << std::ends; 

    std::string msg = outs.str();
    mexPrintf("%s\n", msg.c_str() );
}

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
    try
    {
        // ---------------------------------------------------------------
        // consistency check
        // ---------------------------------------------------------------    
        if (nrhs != NIn)
        {
            std::stringstream outs;
            outs << NIn << " input arguments are required ...";
            std::string msg = outs.str();
            mexWarnMsgTxt(msg.c_str());
            usage();
            return;
        }

        if (nlhs < NOut)
        {
            std::stringstream outs;
            outs << NOut << " output arguments are required ...";
            std::string msg = outs.str();
            mexWarnMsgTxt(msg.c_str());
            usage();
            return;
        }

        using namespace Gadgetron;

        typedef double T;

        // ---------------------------------------------------------------
        // input parameters
        // ---------------------------------------------------------------    
        // dx
        if ( !mxIsDouble(prhs[0]) )
        {
            mexWarnMsgTxt("The first input parameter should be a double array ...");
        }

        hoNDArray<T> dx;
        Gadgetron::MatlabToHoNDArray(const_cast<mxArray*>(prhs[0]), dx);

        if (!mxIsDouble(prhs[1]))
        {
            mexWarnMsgTxt("The second input parameter should be a double array ...");
        }

        hoNDArray<T> dy;
        Gadgetron::MatlabToHoNDArray(const_cast<mxArray*>(prhs[1]), dy);

        if (!dx.dimensions_equal(&dy))
        {
            mexWarnMsgTxt("!dx.dimensions_equal(dy)");
            return;
        }

        size_t N = dx.get_size(2);

        // key_frame
        unsigned int key_frame = (unsigned int)(mxGetScalar(prhs[2]));

        if (key_frame >= N)
        {
            mexWarnMsgTxt("key_frame >= N");
            return;
        }
        // ---------------------------------------------------------------
        // perform the computation
        // ---------------------------------------------------------------

        hoNDArray<double> dx_out, dy_out;
        Gadgetron::concatenate_deform_fields_2DT(dx, dy, key_frame, dx_out, dy_out);

        // ---------------------------------------------------------------
        // output parameter
        // ---------------------------------------------------------------

        mxArray* M0 = NULL;
        M0 = Gadgetron::hoNDArrayToMatlab(&dx_out);
        plhs[0] = M0;

        mxArray* M1 = NULL;
        M1 = Gadgetron::hoNDArrayToMatlab(&dy_out);
        plhs[1] = M1;
   }
    catch(...)
    {
        mexWarnMsgTxt("Exceptions happened in Matlab_gt_concatenate_deform_fields(...) ...");
        return;
    }

    return;
}
