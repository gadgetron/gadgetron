/*
  An example of how to estimate DCF
*/
// Gadgetron includes

//#include <gadgetron/hoNDArray_utils.h>

#include "cuNFFT.h"
#include "cuNDFFT.h"
#include "cuNDArray_math.h"
#include "cuNDArray.h"
#include "cuNDArray_math.h"
#include "cuNDArray_operators.h"
#include "NFFTOperator.h"

// Std includes
#include <iostream>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include "GadgetronTimer.h"
#include <boost/program_options.hpp>
#include <boost/asio.hpp>
#include "hoArmadillo.h"
#include <complex>

using namespace Gadgetron;
using namespace Gadgetron::Core;
using namespace std;

// Define desired precision
typedef float _real;
namespace po = boost::program_options;
uint64d2 image_dims_os_;
std::vector<size_t> image_dims_;
float kernel_width_, oversampling_factor_;
boost::shared_ptr<cuNFFT_plan<float, 2>> nfft_plan_;
boost::shared_ptr<cuNDArray<float_complext>> cuData;

int main(int argc, char **argv)
{
    size_t RO, INT, Kz, CHA; 
    RO = 2500;
    INT = 1100;
    CHA = 8;

    kernel_width_ = 3;
    oversampling_factor_ = 2.1;
    size_t xsize_ = 256;
    size_t ysize_ = 256;
    
    unsigned int warp_size = cudaDeviceManager::Instance()->warp_size();

    std::vector<size_t> data_dims = {RO, INT};
    cuNDArray<float_complext> fake_data(data_dims); 
    
    hoNDArray<floatd2> fake_traj_ho(data_dims);
    cuNDArray<float>   fake_dcw(data_dims);
    
    vector_td<float, 2> init_val;
    init_val[0]= 0.1f;
    init_val[1]= 0.1f;
    
    
    fake_traj_ho.fill(init_val);

    fill(&fake_dcw,1.0f);
    fill(&fake_data,complext(1.0f,0.0f));

    auto fake_traj = cuNDArray<floatd2>(fake_traj_ho);


    image_dims_.push_back(xsize_);
    image_dims_.push_back(ysize_);

    image_dims_os_ = uint64d2(((static_cast<size_t>(std::ceil(image_dims_[0] * oversampling_factor_)) + warp_size - 1) / warp_size) * warp_size,
                              ((static_cast<size_t>(std::ceil(image_dims_[1] * oversampling_factor_)) + warp_size - 1) / warp_size) * warp_size); // No oversampling is needed in the z-direction for SOS

    nfft_plan_ = NFFT<cuNDArray, float, 2>::make_plan(from_std_vector<size_t, 2>(image_dims_), image_dims_os_, kernel_width_, ConvolutionType::STANDARD);

    std::vector<size_t> flat_dims = {fake_traj.get_number_of_elements()};
    cuNDArray<floatd2> flat_traj(flat_dims, fake_traj.get_data_ptr());
    nfft_plan_->preprocess(flat_traj, NFFT_prep_mode::NC2C);

     std::vector<size_t> recon_dims = {image_dims_[0], image_dims_[1]};
     auto temp = boost::make_shared<cuNDArray<float_complext>>(recon_dims);

//     //cuNDArray<float_complext> tempdata_gpu1({RO, E1E2, CHA}, 1); // Tricks to save memory and allow calculations

//     using namespace Gadgetron::Indexing;
    
 
 {
     GadgetronTimer timer("Reconstruct");
     nfft_plan_->compute(&fake_data, *temp, &fake_dcw, NFFT_comp_mode::BACKWARDS_NC2C);
 }

    std::exit(0);
    return 0;
}
