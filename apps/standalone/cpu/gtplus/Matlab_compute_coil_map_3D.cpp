
#include <matrix.h>
#include <mat.h>
#ifdef _WIN32
    #include <mexGT.h>
#else
    #include <mex.h>
#endif // _WIN32

#include "gtPlusISMRMRDReconUtil.h"
#include "gtMatlabConverter.h"
#include "gtMatlabConverterComplex.h"

#define MEXPRINTF(name) mexPrintf(#name);

static void usage()
{
    using namespace std;
    std::ostrstream outs;

    outs << "==============================================================================================" << endl;
    outs << "Usage: compute_coil_map_3D \n";
    outs << "6 Input paras:" << endl;
    outs << '\t' << "complexIm  : RO*E1*E2*CHA*N, 3D complex image array, in complex float" << endl;
    outs << '\t' << "algo       : ISMRMRD_SOUHEIL or ISMRMRD_SOUHEIL_ITER" << endl;
    outs << '\t' << "ks         : kernel size, used by both methods" << endl;
    outs << '\t' << "power      : number of times to perform power method, used by ISMRMRD_SOUHEIL" << endl;
    outs << '\t' << "iterNum    : number of maximal iteration times, used by ISMRMRD_SOUHEIL_ITER" << endl;
    outs << '\t' << "thres      : threshold of iteration, used by ISMRMRD_SOUHEIL_ITER" << endl;

    outs << "1 Output para:" << endl;
    outs << '\t' << "coilMap    : RO*E1*E2*CHA*N coil map" << endl;
    outs << "==============================================================================================" << endl;
    outs << std::ends; 

    mexPrintf("%s\n", outs.str() );
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

        Gadgetron::gtMatlabConverter<float> converter;
        Gadgetron::gtMatlabConverterComplex<ValueType> converterComplex;

        Gadgetron::gtPlus::gtPlusISMRMRDReconUtilComplex<ValueType> gtPlus_util_complex_;

        // ---------------------------------------------------------------
        // input parameters
        // ---------------------------------------------------------------    
        // target
        if ( !mxIsSingle(prhs[0]) || !mxIsComplex(prhs[0]) )
        {
            mexWarnMsgTxt("The first input parameter should be a complex single array ...");
        }

        mwSize nDim = mxGetNumberOfDimensions(prhs[0]);
        if ( nDim!=4 && nDim!=5 )
        {
            mexWarnMsgTxt("1st array is not a 4D or 5D array");
            return;
        }

        const mwSize* dims = mxGetDimensions(prhs[0]);

        // algo
        Gadgetron::gtPlus::ISMRMRDCOILMAPALGO algo = Gadgetron::gtPlus::ISMRMRD_SOUHEIL_ITER;
        std::string algoStr;
        converter.Matlab2Str(prhs[1], algoStr);
        if ( algoStr == "ISMRMRD_SOUHEIL" )
        {
            algo = Gadgetron::gtPlus::ISMRMRD_SOUHEIL;
        }

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
        converterComplex.Matlab2hoNDArray(prhs[0], complexIm);

        Gadgetron::hoNDArray<ValueType> coilMap;

        if ( !gtPlus_util_complex_.coilMap3DNIH(complexIm, coilMap, algo, ks, power, iterNum, thres) )
        {
            mexWarnMsgTxt("coilMap3DNIH(...) failed ... ");
            return;
        }

        // ---------------------------------------------------------------
        // output parameter
        // ---------------------------------------------------------------
        mxArray* coilMapMx = NULL;
        converterComplex.hoNDArray2Matlab(coilMap, coilMapMx);
        plhs[0] = coilMapMx;
   }
    catch(...)
    {
        mexWarnMsgTxt("Exceptions happened in Matlab compute_coil_map_3D() ...");
        return;
    }

    return;
}
