
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

#include "hoImageRegContainer2DRegistration.h"

#define MEXPRINTF(name) mexPrintf(#name);

#define NIn 3
#define NOut 1

static void usage()
{
    using namespace std;
    std::stringstream outs;

    outs << "==============================================================================================" << endl;
    outs << "Usage: Matlab_gt_apply_deformation_field_reg_2D_series \n";
    outs << "Apply gadgetron 2D deformation fields for a series of 2D images" << endl;
    outs << NIn << " Input paras:" << endl;
    outs << '\t' << "data                                       : RO*E1*N, 2D array of images, in float" << endl;
    outs << '\t' << "dx                                         : deformation fields for x dimension" << endl;
    outs << '\t' << "dy                                         : deformation fields for y dimension " << endl;

    outs << "1 Output para:" << endl;
    outs << '\t' << "warpped    : RO*E1*N, warpped image series" << endl;
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

        Gadgetron::GadgetronTimer timer("Running gadgetron deformation field wrapping");

        typedef float T;

        typedef Gadgetron::hoNDImage<T, 2> Image2DType;

        // ---------------------------------------------------------------
        // input parameters
        // ---------------------------------------------------------------    
        // data
        size_t ind = 0;

        if (!mxIsSingle(prhs[ind]))
        {
            mexWarnMsgTxt("The first input parameter should be a single array ...");
        }

        hoNDArray<T> data;
        Gadgetron::MatlabToHoNDArray(const_cast<mxArray*>(prhs[ind++]), data);

        hoNDArray<T> dx, dy;

        if (!mxIsSingle(prhs[ind]))
        {
            mexWarnMsgTxt("The second input parameter should be a single array ...");
        }

        Gadgetron::MatlabToHoNDArray(const_cast<mxArray*>(prhs[ind++]), dx);

        if (!mxIsSingle(prhs[ind]))
        {
            mexWarnMsgTxt("The third input parameter should be a single array ...");
        }

        Gadgetron::MatlabToHoNDArray(const_cast<mxArray*>(prhs[ind++]), dy);

        size_t RO = data.get_size(0);
        size_t E1 = data.get_size(1);
        size_t N = data.get_size(2);

        // ---------------------------------------------------------------
        // perform the registration
        // ---------------------------------------------------------------

        std::vector<size_t> col(1, N);

        hoNDImageContainer2D<Image2DType> imageSeries;
        imageSeries.create(col);

        hoNDImageContainer2D< hoNDImage<T, 2> > deformation_field[2];
        deformation_field[0].create(col);
        deformation_field[1].create(col);

        std::vector<size_t> dimIm(2);
        dimIm[0] = RO;
        dimIm[1] = E1;

        for (size_t p = 0; p<col[0]; p++)
        {
            imageSeries(0, p).create(dimIm, data.begin() + p*RO*E1, false);
            deformation_field[0](0, p).create(dimIm, dx.begin() + p*RO*E1, false);
            deformation_field[1](0, p).create(dimIm, dy.begin() + p*RO*E1, false);
        }

        Gadgetron::hoImageRegContainer2DRegistration<T, float, 2, 2> regContainer;

        hoNDImageContainer2D< Image2DType > warppedContainer;

        regContainer.warpContainer2D(imageSeries, imageSeries, deformation_field, warppedContainer);

        hoNDArray<T> warpped;
        warpped.create(RO, E1, N);
        Gadgetron::clear(warpped);

        warppedContainer.to_NDArray(0, warpped);

        // ---------------------------------------------------------------
        // output parameter
        // ---------------------------------------------------------------

        mxArray* M0 = Gadgetron::hoNDArrayToMatlab(&warpped);
        plhs[0] = M0;
   }
    catch(...)
    {
        mexWarnMsgTxt("Exceptions happened in Matlab_gt_apply_deformation_field_reg_2D_series(...) ...");
        return;
    }

    return;
}
