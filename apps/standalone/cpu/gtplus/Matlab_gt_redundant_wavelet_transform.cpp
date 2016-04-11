
#include <matrix.h>
#include <mat.h>
#ifdef _WIN32
    #include <mexGT.h>
#else
    #include <mex.h>
#endif // _WIN32

// Gadgetron includes
#include "gtMatlab.h"
#include "gtPlusISMRMRDReconUtil.h"
#include "hoNDRedundantWavelet.h"
#include "hoNDHarrWavelet.h"

#include "gtMatlabConverter.h"
#include "gtMatlabConverterComplex.h"

#define MEXPRINTF(name) mexPrintf(#name);

#define NIn 5
#define NOut 1

static void usage()
{
    using namespace std;
    std::stringstream outs;

    outs << "==============================================================================================" << endl;
    outs << "Usage: Matlab_gt_redundant_wavelet_transform \n";
    outs << "Perform gadgetron redundant wavelet forward/inverse transform for a series of 1D/2D/3D signals" << endl;
    printAuthorInfo(outs);
    outs << NIn << " Input paras:" << endl;
    outs << '\t' << "data                                       : RO*E1*E2*... , array of images or wavelet coefficients, double or complex double" << endl;
    outs << '\t' << "wavelet_name                               : wavelet name, db1, db2, db3, db4, db5" << endl;
    outs << '\t' << "transform_dim                              : 1, or 2, or 3 for 1D/2D/3D wavelet transform" << endl;
    outs << '\t' << "forward                                    : 1 for forward transform and -1 for inverse transform" << endl;
    outs << '\t' << "level                                      : number of wavelet transform levels" << endl;

    outs << "1 Output para:" << endl;
    outs << '\t' << "res                                        : wavelet coefficients or images" << endl;
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
        using namespace Gadgetron::gtPlus;

        Gadgetron::GadgetronTimer timer("Running gadgetron wavelet transform");

        // ---------------------------------------------------------------
        // input parameters
        // ---------------------------------------------------------------    
        // target
        if ( !mxIsDouble(prhs[0]) )
        {
            mexWarnMsgTxt("The first input parameter should be a complex double array ...");
            return;
        }

        bool isComplex = false;
        if ( mxIsComplex(prhs[0]) )
        {
            isComplex = true;
        }

        typedef double T;
        Gadgetron::gtMatlabConverter<T> converter;
        Gadgetron::gtMatlabConverterComplex< std::complex<T> > converterComplex;

        // wavelet_name
        std::string wavName;
        converter.Matlab2Str(prhs[1], wavName);

        // transform dimension
        size_t NDim = (size_t)mxGetScalar(prhs[2]);
        if (NDim < 1) NDim = 1;
        if (NDim > 3) NDim = 3;

        // forward
        int forward = mxGetScalar(prhs[3]);

        // level
        int level = mxGetScalar(prhs[4]);
        if (level < 1)
        {
            mexWarnMsgTxt("The level of transform should be at least 1 ...");
            return;
        }

        if (isComplex)
        {
            hoNDArray< std::complex<double> > in;
            converterComplex.Matlab2hoNDArray(prhs[0], in);

            hoNDArray< std::complex<double> > res;
            Gadgetron::hoNDHarrWavelet< std::complex<double> > harrWav;
            Gadgetron::hoNDRedundantWavelet< std::complex<double> > wav;

            if (wavName == "db1" || wavName == "haar")
            {
                harrWav.transform(in, res, NDim, level, (forward > 0));
            }
            else
            {
                wav.compute_wavelet_filter(wavName);
                wav.transform(in, res, NDim, level, (forward > 0));
            }

            mxArray* M0 = NULL;
            converterComplex.hoNDArray2Matlab(res, M0);
            plhs[0] = M0;
        }
        else
        {
            hoNDArray<T> in;
            converter.Matlab2hoNDArray(prhs[0], in);

            hoNDArray<T> res;

            Gadgetron::hoNDHarrWavelet<T> harrWav;
            Gadgetron::hoNDRedundantWavelet<T> wav;

            if (wavName == "db1" || wavName == "haar")
            {
                harrWav.transform(in, res, NDim, level, (forward > 0));
            }
            else
            {
                wav.compute_wavelet_filter(wavName);
                wav.transform(in, res, NDim, level, (forward > 0));
            }

            mxArray* M0 = NULL;
            converter.hoNDArray2Matlab(res, M0);
            plhs[0] = M0;
        }
   }
    catch(...)
    {
        mexWarnMsgTxt("Exceptions happened in Matlab_gt_redundant_wavelet_transform(...) ...");
        return;
    }

    return;
}
