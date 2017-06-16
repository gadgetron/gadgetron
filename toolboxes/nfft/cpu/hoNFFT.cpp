#include "hoNFFT.h"
#include "hoNDFFT.h"
#include "hoNDArray_elemwise.h"
#include "hoNDArray_utils.h"
#include "vector_td_utilities.h"
#include "vector_td_io.h"

#include <stdio.h>
#include <math.h>
#include <cmath>
#include <sstream>
#include <stdexcept>

using std::vector;
using namespace Gadgetron;

// #include "KaiserBessel_kernel.cpp"
// #include "NFFT_C2NC_conv_kernel.cpp"
// #include "NFFT_NC2C_conv_kernel.cpp"
// #include "NFFT_preprocess_kernel.cpp"

template<class REAL, unsigned int D> 
Gadgetron::hoNFFT_plan<REAL, D>::hoNFFT_plan(){
	barebones();
}

template<class REAL, unsigned int D>
Gadgetron::hoNFFT_plan<REAL, D>::hoNFFT_plan(
	typename uint64d<D>::Type matrix_size,
	typename uint64d<D>::Type matrix_size_os,
	REAL W
){
	std::cout << "\nREACHED HERE" << std::endl;
	std::cout << std::endl << (matrix_size) << std::endl << std::endl;
	barebones();
	setup(matrix_size, matrix_size_os, W);
}

template<class REAL, unsigned int D>
void Gadgetron::hoNFFT_plan<REAL, D>::setup(
	typename uint64d<D>::Type matrix_size,
	typename uint64d<D>::Type matrix_size_os,
	REAL W
){
	// todo
}

template<class REAL, unsigned int D>
void Gadgetron::hoNFFT_plan<REAL, D>::preprocess(
	hoNDArray<typename reald<REAL, D>::Type> *trajectory,
	NFFT_prep_mode mode
){
	// todo
}

template<class REAL, unsigned int D>
void Gadgetron::hoNFFT_plan<REAL, D>::compute(
	hoNDArray<complext<REAL>> *in,
	hoNDArray<complext<REAL>> *out,
	hoNDArray<REAL> *dcw,
	NFFT_comp_mode mode
){
	// todo
}

template<class REAL, unsigned int D>
void Gadgetron::hoNFFT_plan<REAL, D>::mult_MH_M(
	hoNDArray<complext<REAL>> *in,
	hoNDArray<complext<REAL>> *out,
	hoNDArray<REAL> *dcw,
	std::vector<size_t> halfway_dims
){
	// todo	
}

template<class REAL, unsigned int D>
void Gadgetron::hoNFFT_plan<REAL, D>::convolve(
	hoNDArray<complext<REAL>> *in,
	hoNDArray<complext<REAL>> *out,
	hoNDArray<REAL>* dcw,
	NFFT_conv_mode mode,
	bool accumulate
){
	// todo
}

template<class REAL, unsigned int D>
void Gadgetron::hoNFFT_plan<REAL, D>::fft(
	hoNDArray<complext<REAL>> *data,
	NFFT_fft_mode mode,
	bool do_scale
){
	// todo
}

template<class REAL, unsigned int D>
void Gadgetron::hoNFFT_plan<REAL, D>::deapodize(
	hoNDArray<complext<REAL>> *image,
	bool fourier_domain
){
	// todo	
}

template<class REAL, unsigned int D>
void Gadgetron::hoNFFT_plan<REAL, D>::check_consistency(
	hoNDArray<complext<REAL>> *samples,
	hoNDArray<complext<REAL>> *image,
	hoNDArray<REAL> *dcw
){
	// todo
}

template<class REAL, unsigned int D>
void Gadgetron::hoNFFT_plan<REAL, D>::barebones(){
	std::cout << "\n\n\n\n barebones \n\n\n\n" << std::endl;
}

template<class REAL, unsigned int D>
void Gadgetron::hoNFFT_plan<REAL, D>::compute_beta(){
	// todo
}

template<class REAL, unsigned int D>
boost::shared_ptr<hoNDArray<complext<REAL>>> Gadgetron::hoNFFT_plan<REAL, D>::compute_deapodization_filter(bool FFTed){
	// todo 
}

template<class REAL, unsigned int D>
void Gadgetron::hoNFFT_plan<REAL, D>::compute_NFFT_C2NC(
	hoNDArray<complext<REAL>> *in,
	hoNDArray<complext<REAL>> *out
){
	deapodize(in);
	fft(in, NFFT_FORWARDS);
	convolve(in, out, 0x0, NFFT_CONV_C2NC);
}

template<class REAL, unsigned int D>
void Gadgetron::hoNFFT_plan<REAL, D>::compute_NFFTH_NC2C(
	hoNDArray<complext<REAL>> *in,
	hoNDArray<complext<REAL>> *out
){
	convolve(in, out, 0x0, NFFT_CONV_NC2C);
	fft(out, NFFT_BACKWARDS);
	deapodize(out);
}

template<class REAL, unsigned int D>
void Gadgetron::hoNFFT_plan<REAL, D>::compute_NFFTH_C2NC(
	hoNDArray<complext<REAL>> *in,
	hoNDArray<complext<REAL>> *out
){
	deapodize(in, true);
	fft(in, NFFT_BACKWARDS);
	convolve(in, out, 0x0, NFFT_CONV_C2NC);
}

template<class REAL, unsigned int D>
void Gadgetron::hoNFFT_plan<REAL, D>::compute_NFFT_NC2C(
	hoNDArray<complext<REAL>> *in,
	hoNDArray<complext<REAL>> *out
){
	convolve(in, out, 0x0, NFFT_CONV_NC2C);
	fft(out, NFFT_FORWARDS);
	deapodize(out, true);
}

template<class REAL, unsigned int D>
void Gadgetron::hoNFFT_plan<REAL, D>::convolve_NFFT_C2NC(
	hoNDArray<complext<REAL>> *in,
	hoNDArray<complext<REAL>> *out,
	bool accumulate
){
	// todo	
}

template<class REAL, unsigned int D>
void Gadgetron::hoNFFT_plan<REAL, D>::convolve_NFFT_NC2C(
	hoNDArray<complext<REAL>> *in,
	hoNDArray<complext<REAL>> *out,
	bool accumulate
){
	// todo
}

template<class REAL, unsigned int D>
void Gadgetron::hoNFFT_plan<REAL, D>::image_wrap(
	hoNDArray<complext<REAL>> *in,
	hoNDArray<complext<REAL>> *out,
	bool accumulate
){
	// todo
}

template EXPORTCPUNFFT class Gadgetron::hoNFFT_plan<float, 1>;
template EXPORTCPUNFFT class Gadgetron::hoNFFT_plan<double, 1>;
template EXPORTCPUNFFT class Gadgetron::hoNFFT_plan<float, 2>;
template EXPORTCPUNFFT class Gadgetron::hoNFFT_plan<double, 2>;
template EXPORTCPUNFFT class Gadgetron::hoNFFT_plan<float, 3>;
template EXPORTCPUNFFT class Gadgetron::hoNFFT_plan<double, 3>;
template EXPORTCPUNFFT class Gadgetron::hoNFFT_plan<float, 4>;
template EXPORTCPUNFFT class Gadgetron::hoNFFT_plan<double, 4>;
