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
#include "cuNonCartesianMOCOOperator.h"
#include "cuNDArray.h"
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
     boost::shared_ptr<hoNDArray<std::complex<float>>> moving_image =
        read_nd_array<std::complex<float>>((char*)parms.get_parameter('m')->get_string_value());

    boost::shared_ptr<hoNDArray<float>> vfield =
        read_nd_array<float>((char*)parms.get_parameter('v')->get_string_value());

    auto mimage = cuNDArray<float_complext>(hoNDArray<float_complext>(*moving_image));
    auto cuvfield = cuNDArray<float>(*vfield);

    boost::shared_ptr<cuNonCartesianMOCOOperator<float, 3>> E3_;
    E3_ = boost::shared_ptr<cuNonCartesianMOCOOperator<float, 3>>(new cuNonCartesianMOCOOperator<float, 3>(ConvolutionType::ATOMIC));

    E3_->deform_image(&mimage,cuvfield);

    hoNDArray<float_complext> out;
    auto result_str = std::string(parms.get_parameter('r')->get_string_value());

    write_nd_array<float_complext>(&out, (result_str + std::string("_warped_image.complex")).c_str());

    // ---------------------------------------------------------------


    return 0;
}
