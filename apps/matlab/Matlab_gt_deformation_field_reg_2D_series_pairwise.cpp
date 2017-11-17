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

// Gadgetron includes
#include "hoImageRegContainer2DRegistration.h"

#define MEXPRINTF(name) mexPrintf(#name);

#define NIn 15
#define NOut 2

static void usage()
{
    using namespace std;
    std::stringstream outs;

    outs << "==============================================================================================" << endl;
    outs << "Usage: Matlab_gt_deformation_field_reg_2D_series_pairwise \n";
    outs << "Perform gadgetron 2D non-rigid image registration for two series of 2D images" << endl;
    outs << NIn << " Input paras:" << endl;
    outs << '\t' << "target                                     : RO*E1*N, 2D target array of images, in float" << endl;
    outs << '\t' << "source                                     : RO*E1*N, 2D source array of images, in float" << endl;
    outs << '\t' << "dissimilarity                              : 'SSD' or 'LocalCCR' or 'MutualInformation' " << endl;
    outs << '\t' << "level                                      : number of resolution levels" << endl;
    outs << '\t' << "max_iter_num_pyramid_level                 : number of maximal iterations for every resolution level, from high resolution to low resolution, e.g. [16.0 32.0 64.0]" << endl;
    outs << '\t' << "regularization_hilbert_strength            : regularization strength for every resolution level, in pixel unit, e.g. [12.0 12.0 12.0]" << endl;
    outs << '\t' << "LocalCCR_sigmaArg                          : for LocalCCR, the locality of cost function integral, in pixel unit, e.g. [2.0 2.0 2.0]" << endl;
    outs << '\t' << "BidirectionalReg                           : 0 or 1, if 1, run the bidirectional registration" << endl;
    outs << '\t' << "dissimilarity_thres                        : threshold for dissimilarity measures, e.g. 1e-5" << endl;
    outs << '\t' << "div_num                                    : number of search steps, e.g. 3" << endl;
    outs << '\t' << "inverse_deform_enforce_iter                : number of bidirectional iterations, e.g. 10" << endl;
    outs << '\t' << "inverse_deform_enforce_weight              : weight for bidirectional iteration, e.g. 0.5" << endl;
    outs << '\t' << "DivergenceFreeReg                          : 0 or 1, if 1, apply the divergence free constraint" << endl;
    outs << '\t' << "DebugFolder                                : if not empty, intermediate results are stored in this folder" << endl;
    outs << '\t' << "verbose                                    : 0 or 1, if 1, more outputs are printed out" << endl;

    outs << "2 - 5 Output para:" << endl;
    outs << '\t' << "dx         : RO*E1*N, deformation field from target to source, along x dimension" << endl;
    outs << '\t' << "dy         : RO*E1*N, deformation field from target to source, along y dimension" << endl;
    outs << '\t' << "warped     : RO*E1*N, warped images" << endl;
    outs << '\t' << "dxInv      : RO*E1*N, deformation field from source to target, along x dimension" << endl;
    outs << '\t' << "dyInv      : RO*E1*N, deformation field from source to target, along y dimension" << endl;
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

        Gadgetron::GadgetronTimer timer("Running gadgetron deformation field registration");

        typedef float T;

        // ---------------------------------------------------------------
        // input parameters
        // ---------------------------------------------------------------    
        // target
        if ( !mxIsSingle(prhs[0]) )
        {
            mexWarnMsgTxt("The first input parameter should be a single array ...");
        }

        Gadgetron::hoNDArray<T> targetArray;
        Gadgetron::MatlabToHoNDArray(const_cast<mxArray*>(prhs[0]), targetArray);

        Gadgetron::hoNDImage<T, 3> target;
        target.from_NDArray(targetArray);

        // key frame
        if ( !mxIsSingle(prhs[1]) )
        {
            mexWarnMsgTxt("The second input parameter should be a single array ...");
        }

        Gadgetron::hoNDArray<T> sourceArray;
        Gadgetron::MatlabToHoNDArray(const_cast<mxArray*>(prhs[1]), sourceArray);

        Gadgetron::hoNDImage<T, 3> source;
        source.from_NDArray(sourceArray);

        // dissimilarity
        Gadgetron::GT_IMAGE_DISSIMILARITY dissimilarity = Gadgetron::GT_IMAGE_DISSIMILARITY_LocalCCR;
        std::string str2;
        str2 = Gadgetron::MatlabToStdString(prhs[2]);
        dissimilarity = Gadgetron::getDissimilarityType(str2);

        // level
        unsigned int level = (unsigned int)(mxGetScalar(prhs[3]));

        // max_iter_num_pyramid_level
        std::vector<T> itersFloat;
        Gadgetron::MatlabToStdVec(prhs[4], itersFloat);

        std::vector<unsigned int> iters(itersFloat.size());
        std::copy(itersFloat.begin(), itersFloat.end(), iters.begin());

        if ( iters.size() != level )
        {
            unsigned int iter1 = iters[0];
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

        // dissimilarity_thres
        T dissimilarity_thres = (T)(mxGetScalar(prhs[8]));

        // div_num
        unsigned int div_num = (unsigned int)(mxGetScalar(prhs[9]));

        // inverse_deform_enforce_iter
        unsigned int inverse_deform_enforce_iter = (unsigned int)(mxGetScalar(prhs[10]));

        // inverse_deform_enforce_weight
        T inverse_deform_enforce_weight = (T)(mxGetScalar(prhs[11]));

        // DivergenceFreeReg
        bool DivergenceFreeReg = true;
        if (mxGetScalar(prhs[12]) == 0)
        {
            DivergenceFreeReg = false;
        }

        // DebugFolder
        std::string debugFolder;
        debugFolder = Gadgetron::MatlabToStdString(prhs[13]);

        // verbose
        bool verbose = false;
        if ( mxGetScalar(prhs[14]) > 0 )
        {
            verbose = true;
        }

        // ---------------------------------------------------------------
        // perform the registration
        // ---------------------------------------------------------------
        typedef float T;

        typedef Gadgetron::hoNDImage<T, 2> Image2DType;

        std::vector<size_t> dim(3);
        dim[0] = target.get_size(0);
        dim[1] = target.get_size(1);
        dim[2] = target.get_size(2);

        Gadgetron::hoNDImageContainer2D<Image2DType> targetSeries;
        targetSeries.create(target.begin(), dim);

        Gadgetron::hoNDImageContainer2D<Image2DType> sourceSeries;
        sourceSeries.create(source.begin(), dim);

        Gadgetron::hoImageRegContainer2DRegistration<Image2DType, Image2DType, float> regContainer;

        regContainer.setDefaultParameters(level, false);

        regContainer.container_reg_mode_ = Gadgetron::GT_IMAGE_REG_CONTAINER_PAIR_WISE;

        regContainer.container_reg_transformation_ = Gadgetron::GT_IMAGE_REG_TRANSFORMATION_DEFORMATION_FIELD;
        if ( BidirectionalReg )
        {
            regContainer.container_reg_transformation_ = Gadgetron::GT_IMAGE_REG_TRANSFORMATION_DEFORMATION_FIELD_BIDIRECTIONAL;
        }

        regContainer.bg_value_ = -1;
        regContainer.use_world_coordinates_ = false;
        regContainer.resolution_pyramid_levels_ = level;
        regContainer.max_iter_num_pyramid_level_ = iters;
        regContainer.dissimilarity_type_ = dissimilarity;

        regContainer.dissimilarity_thres_pyramid_level_.clear();
        regContainer.dissimilarity_thres_pyramid_level_.resize(level, dissimilarity_thres);

        regContainer.div_num_pyramid_level_.clear();
        regContainer.div_num_pyramid_level_.resize(level, div_num);

        regContainer.regularization_hilbert_strength_world_coordinate_ = false;

        // regContainer.apply_divergence_free_constraint_ = DivergenceFreeReg;

        unsigned int ii;
        for ( ii=0; ii<level; ii++ )
        {
            regContainer.regularization_hilbert_strength_pyramid_level_[ii][0] = regularization_hilbert_strength[ii];
            regContainer.regularization_hilbert_strength_pyramid_level_[ii][1] = regularization_hilbert_strength[ii];

            regContainer.dissimilarity_LocalCCR_sigmaArg_[ii][0] = LocalCCR_sigmaArg[ii];
            regContainer.dissimilarity_LocalCCR_sigmaArg_[ii][1] = LocalCCR_sigmaArg[ii];
        }

        regContainer.inverse_deform_enforce_iter_pyramid_level_.clear();
        regContainer.inverse_deform_enforce_iter_pyramid_level_.resize(level, inverse_deform_enforce_iter);

        regContainer.inverse_deform_enforce_weight_pyramid_level_.clear();
        regContainer.inverse_deform_enforce_weight_pyramid_level_.resize(level, inverse_deform_enforce_weight);

        regContainer.debugFolder_ = debugFolder;
        regContainer.verbose_ = verbose;

        /*if ( verbose )
        {
            matlab_printInfo(regContainer);
        }*/

        bool warped = false;
        if ( nlhs > 3 )
        {
            warped = true;
        }

        regContainer.registerOverContainer2DPairWise(targetSeries, sourceSeries, warped);

        // ---------------------------------------------------------------
        // output parameter
        // ---------------------------------------------------------------

        Gadgetron::hoNDArray<float> deformField;

        mxArray* M0 = NULL;
        regContainer.deformation_field_[0].to_NDArray(0, deformField);
        M0 = Gadgetron::hoNDArrayToMatlab(&deformField);
        plhs[0] = M0;

        mxArray* M1 = NULL;
        regContainer.deformation_field_[1].to_NDArray(0, deformField);
        M1 = Gadgetron::hoNDArrayToMatlab(&deformField);
        plhs[1] = M1;

        if ( nlhs > 2 )
        {
            Gadgetron::hoNDArray<T> out;

            mxArray* M2 = NULL;
            regContainer.warped_container_.to_NDArray(0, out);
            M2 = Gadgetron::hoNDArrayToMatlab(&out);
            plhs[2] = M2;
        }

        if ( nlhs > 3 )
        {
            mxArray* M3 = NULL;
            regContainer.deformation_field_inverse_[0].to_NDArray(0, deformField);
            M3 = Gadgetron::hoNDArrayToMatlab(&deformField);
            plhs[3] = M3;
        }

        if ( nlhs > 4 )
        {
            mxArray* M4 = NULL;
            regContainer.deformation_field_inverse_[1].to_NDArray(0, deformField);
            M4 = Gadgetron::hoNDArrayToMatlab(&deformField);
            plhs[4] = M4;
        }
   }
    catch(...)
    {
        mexWarnMsgTxt("Exceptions happened in Matlab_gt_deformation_field_reg_2D_series_pairwise(...) ...");
        return;
    }

    return;
}
