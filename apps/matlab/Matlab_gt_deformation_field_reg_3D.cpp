
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

// Gadgetron includes
#include "hoImageRegContainer2DRegistration.h"

#define MEXPRINTF(name) mexPrintf(#name);

#define NIn 11
#define NOut 7

static void usage()
{
    using namespace std;
    std::stringstream outs;

    outs << "==============================================================================================" << endl;
    outs << "Usage: Matlab_gt_deformation_field_reg_3D \n";
    outs << "Perform gadgetron 3D non-rigid image registration" << endl;
    printAuthorInfo(outs);
    outs << NIn << " Input paras:" << endl;
    outs << '\t' << "target                                     : Nfe*Npe*N, 3D array, target (fixed) image, in double" << endl;
    outs << '\t' << "source                                     : Nfe*Npe*N, 3D array, source (moving) image, in double" << endl;
    outs << '\t' << "dissimilarity                              : 'SSD' or 'LocalCCR' or 'MutualInformation' " << endl;
    outs << '\t' << "level                                      : number of resolution levels" << endl;
    outs << '\t' << "max_iter_num_pyramid_level                 : number of maximal iterations for every resolution level, from high resolution to low resolution, e.g. [16.0 32.0 64.0]" << endl;
    outs << '\t' << "regularization_hilbert_strength            : regularization strength for every resolution level, in pixel unit, e.g. [12.0 12.0 12.0]" << endl;
    outs << '\t' << "LocalCCR_sigmaArg                          : for LocalCCR, the locality of cost function integral, in pixel unit, e.g. [2.0 2.0 2.0]" << endl;
    outs << '\t' << "BidirectionalReg                           : 0 or 1, if 1, run the bidirectional registration" << endl;
    outs << '\t' << "DivergenceFreeReg                          : 0 or 1, if 1, apply the divergence free constraint" << endl;
    outs << '\t' << "DebugFolder                                : if not empty, intermediate results are stored in this folder" << endl;
    outs << '\t' << "verbose                                    : 0 or 1, if 1, more outputs are printed out" << endl;

    outs << NOut << " Output para:" << endl;
    outs << '\t' << "dx         : deformation field from target to source, along x dimension" << endl;
    outs << '\t' << "dy         : deformation field from target to source, along y dimension" << endl;
    outs << '\t' << "dz         : deformation field from target to source, along z dimension" << endl;
    outs << '\t' << "warped     : warped source" << endl;
    outs << '\t' << "dxInv      : deformation field from source to target, along x dimension" << endl;
    outs << '\t' << "dyInv      : deformation field from source to target, along y dimension" << endl;
    outs << '\t' << "dzInv      : deformation field from source to target, along z dimension" << endl;
    outs << "==============================================================================================" << endl;
    outs << '\t' << "Example usage: " << endl;
    outs << '\t' << "dissimilarity = 'LocalCCR'; " << endl;
    outs << '\t' << "level = 4; " << endl;
    outs << '\t' << "max_iter_num_pyramid_level = [32 64 100 100]; " << endl;
    outs << '\t' << "regularization_hilbert_strength = [12 12 12 12]; " << endl;
    outs << '\t' << "LocalCCR_sigmaArg = [2 2 2 2]; " << endl;
    outs << '\t' << "BidirectionalReg = 1; " << endl;
    outs << '\t' << "dissimilarity_thres = 1e-6; " << endl;
    outs << '\t' << "DivergenceFreeReg = 0; " << endl;
    outs << '\t' << "DebugFolder = []; " << endl;
    outs << '\t' << "verbose = 0; " << endl;

    outs << '\t' << "[dx, dy, dz, warped, dxInv, dyInv, dzInv] = Matlab_gt_deformation_field_reg_3D(single(target), single(source), dissimilarity, level, max_iter_num_pyramid_level, regularization_hilbert_strength, LocalCCR_sigmaArg, BidirectionalReg, DivergenceFreeReg, DebugFolder, verbose); " << endl;

    outs << "==============================================================================================" << endl;
    outs << std::ends;

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

        Gadgetron::GadgetronTimer timer("Running gadgetron deformation field registration");

        typedef double T;

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

        hoNDArray<T> targetArray, sourceArray;

        Gadgetron::MatlabToHoNDArray(const_cast<mxArray*>(prhs[0]), targetArray);
        Gadgetron::MatlabToHoNDArray(const_cast<mxArray*>(prhs[1]), sourceArray);

        hoMRImage<T, 3> target, source;

        target.from_NDArray(targetArray);
        source.from_NDArray(sourceArray);

        // dissimilarity
        GT_IMAGE_DISSIMILARITY dissimilarity = GT_IMAGE_DISSIMILARITY_LocalCCR;
        std::string str = MatlabToStdString(prhs[2]);
        dissimilarity = getDissimilarityType(str);

        // level
        unsigned int level = (unsigned int)(mxGetScalar(prhs[3]));

        // max_iter_num_pyramid_level
        std::vector<T> iters;
        Gadgetron::MatlabToStdVec(prhs[4], iters);

        if ( iters.size() != level )
        {
            T iter1 = iters[0];
            iters.resize(level, iter1);
        }

        // regularization_hilbert_strength
        std::vector<T> regularization_hilbert_strength;
        Gadgetron::MatlabToStdVec(prhs[5], regularization_hilbert_strength);

        if ( regularization_hilbert_strength.size() != level )
        {
            T v1 = regularization_hilbert_strength[0];
            regularization_hilbert_strength.resize(level, v1);
        }

        // LocalCCR_sigmaArg
        std::vector<T> LocalCCR_sigmaArg;
        Gadgetron::MatlabToStdVec(prhs[6], LocalCCR_sigmaArg);

        if ( LocalCCR_sigmaArg.size() != level )
        {
            T v1 = LocalCCR_sigmaArg[0];
            LocalCCR_sigmaArg.resize(level, v1);
        }

        // BidirectionalReg
        bool BidirectionalReg = true;
        if ( mxGetScalar(prhs[7]) == 0 )
        {
            BidirectionalReg = false;
        }

        // DivergenceFreeReg
        bool DivergenceFreeReg = true;
        if (mxGetScalar(prhs[8]) == 0)
        {
            DivergenceFreeReg = false;
        }

        // DebugFolder
        std::string debugFolder = MatlabToStdString(prhs[9]);

        // verbose
        bool verbose = false;
        if ( mxGetScalar(prhs[10]) > 0 )
        {
            verbose = true;
        }

        // ---------------------------------------------------------------
        // perform the registration
        // ---------------------------------------------------------------
        typedef hoMRImage<T, 3> ImageType;

        hoNDBoundaryHandlerFixedValue<ImageType> bhFixedValue;
        bhFixedValue.setArray(source);

        hoNDBoundaryHandlerPeriodic<ImageType> bhPeriodic;
        bhPeriodic.setArray(source);

        hoNDBoundaryHandlerBorderValue<ImageType> bhBorderValue;
        bhBorderValue.setArray(source);

        /// linear warp

        hoNDInterpolatorLinear<ImageType> interp;
        interp.setArray(source);
        interp.setBoundaryHandler(bhFixedValue);

        /// bspline warp
        hoNDInterpolatorBSpline<ImageType, 3> interpBSpline(5);
        interpBSpline.setArray(source);
        interpBSpline.setBoundaryHandler(bhFixedValue);

        // warpper
        hoImageRegWarper<ImageType, ImageType, T> warper;
        warper.setBackgroundValue(-1);

        ImageType warped;
        hoMRImage<double, 3> dx, dy, dz, dxInv, dyInv, dzInv;

        if ( BidirectionalReg )
        {
            Gadgetron::hoImageRegDeformationFieldBidirectionalRegister<ImageType, double> reg(level, false, -1);

            if ( !debugFolder.empty() )
            {
                reg.debugFolder_ = debugFolder;
            }

            bool use_world_coordinates = false;
            reg.setDefaultParameters(level, use_world_coordinates);

            unsigned int ii;

            for ( ii=0; ii<level; ii++ )
            {
                reg.max_iter_num_pyramid_level_[ii] = (unsigned int)iters[ii];

                reg.regularization_hilbert_strength_pyramid_level_[ii][0] = regularization_hilbert_strength[ii];
                reg.regularization_hilbert_strength_pyramid_level_[ii][1] = regularization_hilbert_strength[ii];
                reg.regularization_hilbert_strength_pyramid_level_[ii][2] = regularization_hilbert_strength[ii];

                reg.dissimilarity_LocalCCR_sigmaArg_[ii][0] = LocalCCR_sigmaArg[ii];
                reg.dissimilarity_LocalCCR_sigmaArg_[ii][1] = LocalCCR_sigmaArg[ii];
                reg.dissimilarity_LocalCCR_sigmaArg_[ii][2] = LocalCCR_sigmaArg[ii];
            }

            reg.verbose_ = verbose;

            reg.dissimilarity_type_.clear();
            reg.dissimilarity_type_.resize(level, dissimilarity);

            reg.apply_divergence_free_constraint_ = DivergenceFreeReg;

            reg.setTarget(target);
            reg.setSource(source);

            if ( verbose )
            {
                std::ostringstream outs;
                reg.print(outs);

                outs << std::ends;
                std::string msg(outs.str());
                GDEBUG_STREAM(msg.c_str());
            }

            reg.initialize();

            reg.performRegistration();

            const hoMRImage<double, 3>& dx_reg = reg.transform_->getDeformationField(0);
            const hoMRImage<double, 3>& dy_reg = reg.transform_->getDeformationField(1);
            const hoMRImage<double, 3>& dz_reg = reg.transform_->getDeformationField(2);

            const hoMRImage<double, 3>& dxInv_reg = reg.transform_inverse_->getDeformationField(0);
            const hoMRImage<double, 3>& dyInv_reg = reg.transform_inverse_->getDeformationField(1);
            const hoMRImage<double, 3>& dzInv_reg = reg.transform_inverse_->getDeformationField(2);

            dx = dx_reg;
            dy = dy_reg;
            dz = dz_reg;
            dxInv = dxInv_reg;
            dyInv = dyInv_reg;
            dzInv = dzInv_reg;

            warper.setTransformation(*reg.transform_);
            warper.setInterpolator(interpBSpline);

            warper.warp(target, source, use_world_coordinates, warped);
        }
        else
        {
            Gadgetron::hoImageRegDeformationFieldRegister<ImageType, double> reg(level, false, -1);

            if ( !debugFolder.empty() )
            {
                reg.debugFolder_ = debugFolder;
            }

            bool use_world_coordinates = false;
            reg.setDefaultParameters(level, use_world_coordinates);

            unsigned int ii;

            for ( ii=0; ii<level; ii++ )
            {
                reg.max_iter_num_pyramid_level_[ii] = (unsigned int)iters[ii];

                reg.regularization_hilbert_strength_pyramid_level_[ii][0] = regularization_hilbert_strength[ii];
                reg.regularization_hilbert_strength_pyramid_level_[ii][1] = regularization_hilbert_strength[ii];
                reg.regularization_hilbert_strength_pyramid_level_[ii][2] = regularization_hilbert_strength[ii];

                reg.dissimilarity_LocalCCR_sigmaArg_[ii][0] = LocalCCR_sigmaArg[ii];
                reg.dissimilarity_LocalCCR_sigmaArg_[ii][1] = LocalCCR_sigmaArg[ii];
                reg.dissimilarity_LocalCCR_sigmaArg_[ii][2] = LocalCCR_sigmaArg[ii];
            }

            reg.verbose_ = verbose;

            reg.dissimilarity_type_.clear();
            reg.dissimilarity_type_.resize(level, dissimilarity);

            reg.apply_divergence_free_constraint_ = DivergenceFreeReg;

            reg.setTarget(target);
            reg.setSource(source);

            if ( verbose )
            {
                std::ostringstream outs;
                reg.print(outs);
                outs << std::ends;
                std::string msg(outs.str());
                GDEBUG_STREAM(msg.c_str());
            }

            reg.initialize();

            reg.performRegistration();

            const hoMRImage<double, 3>& dx_reg = reg.transform_->getDeformationField(0);
            const hoMRImage<double, 3>& dy_reg = reg.transform_->getDeformationField(1);
            const hoMRImage<double, 3>& dz_reg = reg.transform_->getDeformationField(2);

            dx = dx_reg;
            dy = dy_reg;
            dz = dz_reg;
            dxInv = dx; Gadgetron::clear(dxInv);
            dyInv = dy; Gadgetron::clear(dyInv);
            dzInv = dz; Gadgetron::clear(dzInv);

            warper.setTransformation(*reg.transform_);
            warper.setInterpolator(interpBSpline);

            warper.warp(target, source, use_world_coordinates, warped);
        }

        // ---------------------------------------------------------------
        // output parameter
        // ---------------------------------------------------------------

        hoNDArray<double> deformField;

        mxArray* M0 = NULL;
        dx.to_NDArray(deformField);
        M0 = Gadgetron::hoNDArrayToMatlab(&deformField);
        plhs[0] = M0;

        mxArray* M1 = NULL;
        dy.to_NDArray(deformField);
        M1 = Gadgetron::hoNDArrayToMatlab(&deformField);
        plhs[1] = M1;

        mxArray* M2 = NULL;
        dz.to_NDArray(deformField);
        M2 = Gadgetron::hoNDArrayToMatlab(&deformField);
        plhs[2] = M2;

        // ---------------------------------------------------------------

        hoNDArray<T> out;

        mxArray* M3 = NULL;
        warped.to_NDArray(out);
        M3 = Gadgetron::hoNDArrayToMatlab(&out);
        plhs[3] = M3;

        // ---------------------------------------------------------------

        mxArray* M4 = NULL;
        dxInv.to_NDArray(deformField);
        M4 = Gadgetron::hoNDArrayToMatlab(&deformField);
        plhs[4] = M4;

        mxArray* M5 = NULL;
        dyInv.to_NDArray(deformField);
        M5 = Gadgetron::hoNDArrayToMatlab(&deformField);
        plhs[5] = M5;

        mxArray* M6 = NULL;
        dzInv.to_NDArray(deformField);
        M6 = Gadgetron::hoNDArrayToMatlab(&deformField);
        plhs[6] = M6;
   }
    catch(...)
    {
        mexWarnMsgTxt("Exceptions happened in Matlab_gt_deformation_field_reg_3D(...) ...");
        return;
    }

    return;
}
