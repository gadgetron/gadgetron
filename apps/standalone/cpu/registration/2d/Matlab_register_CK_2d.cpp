
#include <matrix.h>
#include <mat.h>
#include <mex.h>
#include <cmath>
#include <vector>
#include <iostream>
#include <strstream>

// Gadgetron includes
#include "hoCKOpticalFlowSolver.h"
#include "hoLinearResampleOperator_eigen.h"
#include "hoNDArray.h"
#include "hoNDArray_fileio.h"
#include "GadgetronTimer.h"
#include "parameterparser.h"

#define MEXPRINTF(name) mexPrintf(#name);

static void usage()
{
    using namespace std;
    std::ostrstream outs;

    outs << "==============================================================================================" << endl;
    outs << "Usage: register_CK_2d \n";
    outs << "5 Input paras:" << endl;
    outs << '\t' << "target     : Nfe*Npe, 2D array, target (fixed) image, in double" << endl;
    outs << '\t' << "source     : Nfe*Npe, 2D array, source (moving) image, in double" << endl;
    outs << '\t' << "alpha      : regularization parameter, alpha" << endl;
    outs << '\t' << "beta       : regularization parameter, beta" << endl;
    outs << '\t' << "level      : number of resolution levels" << endl;

    outs << "2 Output para:" << endl;
    outs << '\t' << "dx         : deformation field, along 1st dimension" << endl;
    outs << '\t' << "dy         : deformation field, along 2nd dimension" << endl;
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
        if (nrhs != 5) 
        {
            mexWarnMsgTxt("5 input arguments are required ...");
            usage();
            return;
        }

        if (nlhs < 2 )
        {
            mexWarnMsgTxt("2 output argument is required ...");
            usage();
            return;
        }

        Gadgetron::GadgetronTimer timer("Running registration");

        // ---------------------------------------------------------------
        // input parameters
        // ---------------------------------------------------------------    
        // target
        if ( !mxIsDouble(prhs[0]) )
        {
            mexWarnMsgTxt("The first input parameter should be a double array ...");
        }

        if ( !mxIsDouble(prhs[1]) )
        {
            mexWarnMsgTxt("The second input parameter should be a double array ...");
        }

        // for the image
        mwSize nDim = mxGetNumberOfDimensions(prhs[0]);
        if ( nDim!=2 )
        {
            mexWarnMsgTxt("1st array is not a 2D array");
            return;
        }

        nDim = mxGetNumberOfDimensions(prhs[1]);
        if ( nDim!=2 )
        {
            mexWarnMsgTxt("2nd array is not a 2D array");
            return;
        }

        const mwSize* dims = mxGetDimensions(prhs[0]);
        int numOfPixels = dims[0]*dims[1];

        const mwSize* dims2 = mxGetDimensions(prhs[1]);
        if ( dims[0]!=dims2[0] || dims[1]!=dims2[1] )
        {
            mexWarnMsgTxt("Input arrays have different size ... ");
            return;
        }

        double* ptrTarget = static_cast<double*>(mxGetData(prhs[0]));
        double* ptrSource = static_cast<double*>(mxGetData(prhs[1]));

        // alpha
        double alpha = mxGetScalar(prhs[2]);

        // beta
        double beta = mxGetScalar(prhs[3]);

        // level
        int level = (int)(mxGetScalar(prhs[4]));

        // ---------------------------------------------------------------
        // perform the registration
        // ---------------------------------------------------------------
        // allocate the results
        mxArray* Dx = mxCreateNumericArray(nDim, dims, mxDOUBLE_CLASS, mxREAL);
        if ( Dx == NULL )
        {
            mexWarnMsgTxt("Dx == NULL");
            return;
        }

        mxArray* Dy = mxCreateNumericArray(nDim, dims, mxDOUBLE_CLASS, mxREAL);
        if ( Dy == NULL )
        {
            mexWarnMsgTxt("Dy == NULL");
            return;
        }

        double* ptrDx = static_cast<double*>(mxGetData(Dx));
        double* ptrDy = static_cast<double*>(mxGetData(Dy));
        memset(ptrDx, 0, sizeof(double)*numOfPixels);
        memset(ptrDy, 0, sizeof(double)*numOfPixels);

        // allocate the target and source images
        typedef double _real;
        using namespace Gadgetron;

        std::vector<unsigned int> dim_array(2);
        dim_array[0] = dims[0];
        dim_array[1] = dims[1];

        boost::shared_ptr< hoNDArray<_real> > fixed_image(new hoNDArray<_real>(&dim_array));
        memcpy(fixed_image->begin(), ptrTarget, sizeof(_real)*numOfPixels);

        boost::shared_ptr< hoNDArray<_real> > moving_image(new hoNDArray<_real>(&dim_array));
        memcpy(moving_image->begin(), ptrSource, sizeof(_real)*numOfPixels);

        boost::shared_ptr< hoLinearResampleOperator_eigen<_real,2> > R( new hoLinearResampleOperator_eigen<_real,2>() );

        // Setup solver
        hoCKOpticalFlowSolver<_real,2> CK;
        CK.set_interpolator( R );
        CK.set_output_mode( hoCKOpticalFlowSolver<_real,2>::OUTPUT_VERBOSE );
        CK.set_num_multires_levels( level );
        CK.set_max_num_iterations_per_level( 500 );
        CK.set_alpha(alpha);
        CK.set_beta(beta);
        CK.set_limit(0.01f);
  
        // Run registration
        //
        boost::shared_ptr< hoNDArray<_real> > result;
        result = CK.solve( fixed_image.get(), moving_image.get() );

        if( !result.get() )
        {
            mexWarnMsgTxt("Registration solver failed. Quitting!");
            return;
        }

        memcpy(ptrDx, result->begin(), sizeof(_real)*numOfPixels);
        memcpy(ptrDy, result->begin()+numOfPixels, sizeof(_real)*numOfPixels);

        // ---------------------------------------------------------------
        // output parameter
        // ---------------------------------------------------------------
        plhs[0] = Dx;
        plhs[1] = Dy;
   }
    catch(...)
    {
        mexWarnMsgTxt("Exceptions happened in Matlab register_CK_2d() ...");
        return;
    }

    return;
}
