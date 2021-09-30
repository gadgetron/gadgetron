/*
  An example of how to register two 3d volumes using Cornelius-Kanade optical flow
*/

// Gadgetron includes
#include "cuNDArray.h"
#include "cuNonCartesianMOCOOperator.h"
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
    parms.add_parameter('m', COMMAND_LINE_STRING, 1, "Moving image file name (.complex)", true);
    parms.add_parameter('v', COMMAND_LINE_STRING, 1, "Vfield (.real)", true);
    parms.add_parameter('r', COMMAND_LINE_STRING, 1, "Result file name", true, "path/prefix");
    parms.add_parameter('g', COMMAND_LINE_INT, 1, "use GPU", false, "1");

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
    auto useGPU = parms.get_parameter('g')->get_int_value();
    hoNDArray<float_complext> output;
    boost::shared_ptr<hoNDArray<std::complex<float>>> moving_image =
        read_nd_array<std::complex<float>>((char*)parms.get_parameter('m')->get_string_value());

    boost::shared_ptr<hoNDArray<float>> vfield =
        read_nd_array<float>((char*)parms.get_parameter('v')->get_string_value());

    if (useGPU) {

        auto mimage = cuNDArray<float_complext>(hoNDArray<float_complext>(*moving_image));
        auto cuvfield = cuNDArray<float>(*vfield);
        boost::shared_ptr<cuNonCartesianMOCOOperator<float, 3>> E3_;
        E3_ = boost::shared_ptr<cuNonCartesianMOCOOperator<float, 3>>(
            new cuNonCartesianMOCOOperator<float, 3>(ConvolutionType::ATOMIC));

        E3_->deform_image(&mimage, cuvfield);
        output = hoNDArray<float_complext>(
        std::move(*boost::reinterpret_pointer_cast<hoNDArray<std::complex<float>>>(mimage.to_host())));
    } else {
        typedef hoMRImage<double, 3> ImageType;
        ImageType target, source, deftemp;

        using namespace Gadgetron::Indexing;
        Gadgetron::hoImageRegDeformationField<double, 3> deformation_IM;

        // auto tdef = hoNDArray<double>(hoNDArray<float>((*vfield)(slice,slice,slice,0)));
        // deftemp.from_NDArray(std::move(tdef));
        deftemp.from_NDArray(hoNDArray<double>(hoNDArray<float>((*vfield)(slice, slice, slice, 0))));
        deformation_IM.setDeformationField(deftemp, 0);
        deftemp.from_NDArray(hoNDArray<double>(hoNDArray<float>((*vfield)(slice, slice, slice, 1))));
        deformation_IM.setDeformationField(deftemp, 1);
        deftemp.from_NDArray(hoNDArray<double>(hoNDArray<float>((*vfield)(slice, slice, slice, 2))));
        deformation_IM.setDeformationField(deftemp, 2);
        cout << "made it here" << endl;

        hoNDBoundaryHandlerFixedValue<ImageType> bhFixedValue;
        bhFixedValue.setArray(source);

        hoNDBoundaryHandlerPeriodic<ImageType> bhPeriodic;
        bhPeriodic.setArray(source);

        hoNDBoundaryHandlerBorderValue<ImageType> bhBorderValue;
        bhBorderValue.setArray(source);
        cout << "made it 2" << endl;

        /// bspline warp
        hoNDInterpolatorBSpline<ImageType, 3> interpBSpline(5);
        interpBSpline.setBoundaryHandler(bhFixedValue);
        cout << "made it 3" << endl;

        ImageType warped;
        hoImageRegWarper<ImageType, ImageType, double> warper;
        warper.setBackgroundValue(-1);
        warper.setTransformation(deformation_IM);
        warper.setInterpolator(interpBSpline);

        cout << "made it 4" << endl;
        hoNDArray<double> realOut(moving_image->get_dimensions());
        hoNDArray<double> imagOut(moving_image->get_dimensions());

        auto temp = hoNDArray<double>(real(*moving_image));
        source.from_NDArray(temp);
        interpBSpline.setArray(source);

        warper.warp(source, source, false, warped);
        warped.to_NDArray(realOut);

        temp = hoNDArray<double>(imag(*moving_image));
        source.from_NDArray(temp);
        interpBSpline.setArray(source);

        warper.warp(source, source, false, warped);
        warped.to_NDArray(imagOut);

        output = hoNDArray<float_complext>(*real_imag_to_complex<double_complext>(&realOut, &imagOut));
    }

    auto result_str = std::string(parms.get_parameter('r')->get_string_value());

    write_nd_array<float_complext>(&output, (result_str + std::string("_warped_image.complex")).c_str());

    // ---------------------------------------------------------------

    return 0;
}
