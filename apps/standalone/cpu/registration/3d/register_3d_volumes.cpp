/*
  An example of how to register two 3d volumes using Cornelius-Kanade optical flow
*/

// Gadgetron includes
#include "hoCKOpticalFlowSolver.h"
#include "hoImageRegContainer2DRegistration.h"
#include "hoLinearResampleOperator.h"
#include "hoNDArray.h"
#include "hoNDArray_fileio.h"
#include "parameterparser.h"

// Std includes
#include <iostream>

using namespace Gadgetron;
using namespace std;

// Define desired precision
typedef float _real;

int main(int argc, char** argv) {
    //
    // Parse command line
    //

    ParameterParser parms;
    parms.add_parameter('f', COMMAND_LINE_STRING, 1, "Fixed image file name (.real)", true);
    parms.add_parameter('m', COMMAND_LINE_STRING, 1, "Moving image file name (.real)", true);
    parms.add_parameter('r', COMMAND_LINE_STRING, 1, "Result file name", true, "displacement_field.real");
    parms.add_parameter('i', COMMAND_LINE_STRING, 1, "Result image name", true, "displacement_image.real");
    parms.add_parameter('a', COMMAND_LINE_STRING, 1, "Sigma file name", true);
    parms.add_parameter('b', COMMAND_LINE_STRING, 1, "Regularization weight (beta) file", true);
    parms.add_parameter('l', COMMAND_LINE_INT, 1, "Number of multiresolution levels", true, "3");
    parms.add_parameter('t', COMMAND_LINE_STRING, 1, "Number of iterations file", true);
    parms.add_parameter('d', COMMAND_LINE_STRING, 1, "Dissimilarity: 'SSD' or 'LocalCCR' or 'MutualInformation'", true,
                        "LocalCCR");

    parms.parse_parameter_list(argc, argv);
    if (parms.all_required_parameters_set()) {
        cout << " Running registration with the following parameters: " << endl;
        parms.print_parameter_list();
    } else {
        cout << " Some required parameters are missing: " << endl;
        parms.print_parameter_list();
        parms.print_usage();
        return 1;
    }

    // Load sample data from disk
    //

    boost::shared_ptr<hoNDArray<double>> sigmas =
        read_nd_array<double>((char*)parms.get_parameter('a')->get_string_value());

    boost::shared_ptr<hoNDArray<double>> regularizers =
        read_nd_array<double>((char*)parms.get_parameter('b')->get_string_value());

    boost::shared_ptr<hoNDArray<double>> iterations =
        read_nd_array<double>((char*)parms.get_parameter('t')->get_string_value());

    typedef double T;

    hoMRImage<T, 3> target, source;

    boost::shared_ptr<hoNDArray<double>> fixed_image =
        read_nd_array<double>((char*)parms.get_parameter('f')->get_string_value());

    boost::shared_ptr<hoNDArray<double>> moving_image =
        read_nd_array<double>((char*)parms.get_parameter('m')->get_string_value());

    if (!fixed_image.get() || !moving_image.get()) {
        cout << endl << "One of the input images is not found. Quitting!\n" << endl;
        return 1;
    }

    target.from_NDArray(*fixed_image);
    source.from_NDArray(*moving_image);

    // dissimilarity
    GT_IMAGE_DISSIMILARITY dissimilarity = GT_IMAGE_DISSIMILARITY_LocalCCR;
    std::string str = parms.get_parameter('d')->get_string_value();
    dissimilarity = getDissimilarityType(str);

    size_t num_fixed_dims = fixed_image->get_number_of_dimensions();
    size_t num_moving_dims = moving_image->get_number_of_dimensions();

    if (!(num_fixed_dims == 3 || num_fixed_dims == 4)) {
        cout << endl << "The fixed image is not three- or four-dimensional. Quitting!\n" << endl;
        return 1;
    }

    if (!(num_moving_dims == 3 || num_moving_dims == 4)) {
        cout << endl << "The moving image is not three- or four-dimensional. Quitting!\n" << endl;
        return 1;
    }

    unsigned int multires_levels = parms.get_parameter('l')->get_int_value();
    // level
    unsigned int level = multires_levels;

    std::vector<T> iters(iterations.get()->get_data_ptr(),iterations.get()->get_data_ptr()+level);
    std::vector<T> regularization_hilbert_strength(regularizers.get()->get_data_ptr(),regularizers.get()->get_data_ptr()+level);
    std::vector<T> LocalCCR_sigmaArg(sigmas.get()->get_data_ptr(),sigmas.get()->get_data_ptr()+level);


    if (iters.size() != level) {
        T iter1 = iters[0];
        iters.resize(level, iter1);
    }
    if (iters.size() != level) {
        T iter1 = iters[0];
        iters.resize(level, iter1);
    }

    // regularization_hilbert_strength

    if (regularization_hilbert_strength.size() != level) {
        T v1 = regularization_hilbert_strength[0];
        regularization_hilbert_strength.resize(level, v1);
    }

    // LocalCCR_sigmaArg

    if (LocalCCR_sigmaArg.size() != level) {
        T v1 = LocalCCR_sigmaArg[0];
        LocalCCR_sigmaArg.resize(level, v1);
    }

    // BidirectionalReg
    bool BidirectionalReg = true;

    // DivergenceFreeReg
    bool DivergenceFreeReg = false;

    // verbose
    bool verbose = false;

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
    hoMRImage<T, 3> dx, dy, dz, dxInv, dyInv, dzInv;
    // Use trilinear interpolation for resampling
    //
    Gadgetron::hoImageRegDeformationFieldBidirectionalRegister<ImageType, T> reg(level, false, -1);

    bool use_world_coordinates = false;
    reg.setDefaultParameters(level, use_world_coordinates);

    unsigned int ii;

    for (ii = 0; ii < level; ii++) {
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

    if (verbose) {
        std::ostringstream outs;
        reg.print(outs);

        outs << std::ends;
        std::string msg(outs.str());
        GDEBUG_STREAM(msg.c_str());
    }
    reg.initialize();

    reg.performRegistration();

    const hoMRImage<T, 3>& dx_reg = reg.transform_->getDeformationField(0);
    const hoMRImage<T, 3>& dy_reg = reg.transform_->getDeformationField(1);
    const hoMRImage<T, 3>& dz_reg = reg.transform_->getDeformationField(2);

    const hoMRImage<T, 3>& dxInv_reg = reg.transform_inverse_->getDeformationField(0);
    const hoMRImage<T, 3>& dyInv_reg = reg.transform_inverse_->getDeformationField(1);
    const hoMRImage<T, 3>& dzInv_reg = reg.transform_inverse_->getDeformationField(2);

    dx = dx_reg;
    dy = dy_reg;
    dz = dz_reg;
    dxInv = dxInv_reg;
    dyInv = dyInv_reg;
    dzInv = dzInv_reg;

    warper.setTransformation(*reg.transform_);
    warper.setInterpolator(interpBSpline);

    warper.warp(target, source, use_world_coordinates, warped);

    // ---------------------------------------------------------------
    // output parameter
    // ---------------------------------------------------------------
    auto size_deformation = *fixed_image->get_dimensions();
    size_deformation.push_back(target.get_number_of_dimensions());
    hoNDArray<T> deformation;
    hoNDArray<T> inv_deformation;
    deformation.create(size_deformation);
    inv_deformation.create(size_deformation);
    
    using namespace Gadgetron::Indexing;
    
    hoNDArray<T> deformField;
    auto result_str = std::string(parms.get_parameter('r')->get_string_value());

    dx_reg.to_NDArray(deformField);
    deformation(slice, slice, slice, 0) = deformField;
    dy_reg.to_NDArray(deformField);
    deformation(slice, slice, slice, 1) = deformField;
    dz_reg.to_NDArray(deformField);
    deformation(slice, slice, slice, 2) = deformField;

    dxInv_reg.to_NDArray(deformField);
    inv_deformation(slice, slice, slice, 0) = deformField;
    dyInv_reg.to_NDArray(deformField);
    inv_deformation(slice, slice, slice, 1) = deformField;
    dzInv_reg.to_NDArray(deformField);
    inv_deformation(slice, slice, slice, 2) = deformField;

    write_nd_array<double>(&deformation, (result_str + std::string("_deformation.double")).c_str());

    write_nd_array<double>(&inv_deformation, (result_str + std::string("_Invdeformation.double")).c_str());

    // ---------------------------------------------------------------

    hoNDArray<T> out;

    warped.to_NDArray(out);
    write_nd_array<double>(&out, (result_str + std::string("_warped_image.double")).c_str());

    // ---------------------------------------------------------------


    return 0;
}
