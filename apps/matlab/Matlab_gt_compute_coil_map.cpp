
#include <matrix.h>
#include <mat.h>

#ifdef MATLAB_DLL_EXPORT_SYM
    #define DLL_EXPORT_SYM extern "C" __declspec(dllexport)
#endif // MATLAB_DLL_EXPORT_SYM
#include <mex.h>

#include <sstream>
#include "mri_core_coil_map_estimation.h"
#include "GadgetronTimer.h"
#include "MatlabUtils.h"
#include "Matlab_info.h"

#define MEXPRINTF(name) mexPrintf(#name);

static void usage()
{
    using namespace std;
    std::stringstream outs;

    outs << "==============================================================================================" << endl;
    outs << "Usage: compute_coil_map \n";
    outs << "Compute coil sensitivity maps" << endl;
    printAuthorInfo(outs);
    outs << "6 Input paras:" << endl;
    outs << '\t' << "complexIm  : RO*E1*E2*CHA*N, 2D (if E2==1) or 3D complex image array, in complex float" << endl;
    outs << '\t' << "algo       : ISMRMRD_SOUHEIL or ISMRMRD_SOUHEIL_ITER" << endl;
    outs << '\t' << "ks         : kernel size, used by both methods" << endl;
    outs << '\t' << "power      : number of times to perform power method, used by ISMRMRD_SOUHEIL" << endl;
    outs << '\t' << "iterNum    : number of maximal iteration times, used by ISMRMRD_SOUHEIL_ITER" << endl;
    outs << '\t' << "thres      : threshold of iteration, used by ISMRMRD_SOUHEIL_ITER" << endl;

    outs << "1 Output para:" << endl;
    outs << '\t' << "coilMap    : RO*E1*E2*CHA*N coil map" << endl;
    outs << "==============================================================================================" << endl;
    outs << std::ends; 

    mexPrintf("%s\n", outs.str().c_str() );
}

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
    try
    {
        // ---------------------------------------------------------------
        // consistency check
        // ---------------------------------------------------------------    
        if (nrhs != 6) 
        {
            mexWarnMsgTxt("6 input arguments are required ...");
            usage();
            return;
        }

        if (nlhs < 1 )
        {
            mexWarnMsgTxt("1 output argument is required ...");
            usage();
            return;
        }

        typedef std::complex<float> ValueType;

        Gadgetron::GadgetronTimer timer("Running coil map estimation");

        // ---------------------------------------------------------------
        // input parameters
        // ---------------------------------------------------------------    
        // target
        if ( !mxIsSingle(prhs[0]) || !mxIsComplex(prhs[0]) )
        {
            mexWarnMsgTxt("The first input parameter should be a complex single array ...");
        }

        mwSize nDim = mxGetNumberOfDimensions(prhs[0]);
        if ( nDim<4 )
        {
            mexWarnMsgTxt("1st array at least should be a 4D array");
            return;
        }

        const mwSize* dims = mxGetDimensions(prhs[0]);

        // algo
        std::string algoStr = Gadgetron::MatlabToStdString(prhs[1]);

        // ks
        unsigned long long ks = mxGetScalar(prhs[2]);

        // power
        unsigned long long power = mxGetScalar(prhs[3]);

        // iterNum
        unsigned long long iterNum = (unsigned long long)(mxGetScalar(prhs[4]));

        // iterNum
        float thres = (float)(mxGetScalar(prhs[5]));

        // ---------------------------------------------------------------
        // perform the computation
        // ---------------------------------------------------------------
        Gadgetron::hoNDArray<ValueType> complexIm;
        Gadgetron::MatlabToHoNDArray( const_cast<mxArray*>(prhs[0]), complexIm);

        Gadgetron::hoNDArray<ValueType> coilMap;

        if(algoStr == "ISMRMRD_SOUHEIL")
        {
            Gadgetron::coil_map_Inati(complexIm, coilMap, ks, ks, power);
        }
        else
        {
            Gadgetron::coil_map_Inati_Iter(complexIm, coilMap, ks, ks, iterNum, thres);
        }

        // ---------------------------------------------------------------
        // output parameter
        // ---------------------------------------------------------------
        mxArray* coilMapMx = Gadgetron::hoNDArrayToMatlab(&coilMap);
        plhs[0] = coilMapMx;
   }
    catch(...)
    {
        mexWarnMsgTxt("Exceptions happened in Matlab compute_coil_map() ...");
        return;
    }

    return;
}
