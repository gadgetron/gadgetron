#include "hoNFFT.h"
#include "hoNDFFT.h"
#include "hoNDArray_elemwise.h"
#include "hoNDArray_fileio.h"
#include "hoNDArray_reductions.h"
#include "hoNDArray_utils.h"
#include "vector_td_utilities.h"
#include "vector_td_operators.h"
#include "vector_td_io.h"

#include <math.h>
#include <algorithm>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <boost/make_shared.hpp>

using namespace std;
using namespace Gadgetron;

void printHighlighted(string s){
	std::cout << "=============================" << std::endl;
	std::cout << s << std::endl;
	std::cout << "=============================" << std::endl;
}

void printDivider(){
	cout << "===========================" << endl;
}

template<typename T>
void printTail(T v){
	for(auto it = v.end()-10; it != v.end(); it++)
		cout << *it << ", ";
	cout << endl;
}

template<class Real>
Real bessi0(Real x){
	Real denominator;
	Real numerator;
	Real z;
	if (x == 0.0) {
  	return 1.0;
  } else {
  	z = x * x;
    numerator = (z* (z* (z* (z* (z* (z* (z* (z* (z* (z* (z* (z* (z* 
    	(z* 0.210580722890567e-22  + 0.380715242345326e-19 ) +
     	0.479440257548300e-16) + 0.435125971262668e-13 ) +
      0.300931127112960e-10) + 0.160224679395361e-7  ) +
      0.654858370096785e-5)  + 0.202591084143397e-2  ) +
      0.463076284721000e0)   + 0.754337328948189e2   ) +
      0.830792541809429e4)   + 0.571661130563785e6   ) +
      0.216415572361227e8)   + 0.356644482244025e9   ) +
     	0.144048298227235e10);

    denominator = (z*(z*(z-0.307646912682801e4)+
      0.347626332405882e7)-0.144048298227235e10);
	}

  return -numerator/denominator;
}

template<class Real, unsigned int D>
hoNDArray<complext<Real>> Gadgetron::hoNFFT_plan<Real, D>::processTesting(
	hoNDArray<typename reald<Real, D>::Type> traj,
	hoNDArray<complext<Real>> data,
	hoNDArray<Real> weights,
	size_t n,
	Real osf,
	Real wg
){
	hoNDArray<complext<Real>> output(osf*n, osf*n);
	output.fill(0);
	hoNDArray<complext<Real>> m(n*osf, n*osf);
	data *= weights;
	
	for(size_t i = 0; i < traj.get_number_of_elements(); i++){
		complext<Real> dw = data[i]; 
		for(int lx = -kwidth; lx < kwidth+1; lx++){
			for(int ly = -kwidth; ly < kwidth+1; ly++){
				Real nxt = std::round(nx[i]+lx);
				Real nyt = std::round(ny[i]+ly);
				
				Real kkx = std::min(
					std::round(kosf*std::abs(nx[i]-nxt)),
					std::floor(kosf*kwidth)
				);
				Real kky = std::min(
					std::round(kosf*std::abs(ny[i]-nyt)),
					std::floor(kosf*kwidth)
				);
				Real kwx = p[kkx]; Real kwy = p[kky];

				nxt = std::max(nxt, Real(0)); nxt = std::min(nxt, osf*n-1);
				nyt = std::max(nyt, Real(0)); nyt = std::min(nyt, osf*n-1);

				output[(size_t)(nxt+nyt*osf*n)] += dw*kwx*kwy;
			}
		}
	}

	for(size_t i = 0; i < n*osf; i++){
		output[i] = 0;
		output[n*osf*i] = 0;
		output[n*osf*(n*osf-1)+i] = 0;
		output[n*osf*i+(n*osf-1)] = 0;
	}
	
	hoNDFFT<Real>::instance()->ifft(&output);
	output /= da;
	multiplyConj(output, output, output);
	
	return output;
}

template<class T>
void printArray(T *arr){
	for(auto it = arr->begin(); it < arr->end(); it += (arr->get_number_of_elements()/100)){
		cout << *it << ", ";
	}
	cout << endl;
}

template<class Real, unsigned int D>
hoNDArray<complext<Real>> Gadgetron::hoNFFT_plan<Real, D>::processAll(
	hoNDArray<typename reald<Real, D>::Type> traj,
	hoNDArray<complext<Real>> data,
	hoNDArray<Real> weights, 
	size_t n,
	Real osf,
	Real wg
){
	prepare(traj, data, weights, n, osf, wg);
	hoNDArray<complext<Real>> output(matrixSizeOs[0], matrixSizeOs[1]);
	size_t trajSize = traj.get_number_of_elements();
	size_t dataSize = data.get_number_of_elements();
	size_t nCoils = dataSize/trajSize;
	for(int i = 0; i < nCoils; i++){
		hoNDArray<complext<Real>> tmpData(trajSize);
		std::copy(data.begin()+(size_t)(i*trajSize), data.begin()+(size_t)((i+1)*trajSize), tmpData.begin());
		auto result = processTesting(traj, tmpData, weights, n, osf, wg);
		add(&output, &result, &output);
	}
	sqrt_inplace(&output);
	auto res = boost::make_shared<hoNDArray<complext<Real>>>(output);
	write_nd_array<complext<Real>>(res.get(), "/home/shresthab2/workspace/tests/output.cplx");
	return output;
}
	


template<class Real, unsigned int D>
Gadgetron::hoNFFT_plan<Real, D>::hoNFFT_plan(){}

template<class Real, unsigned int D>
Gadgetron::hoNFFT_plan<Real, D>::hoNFFT_plan(
	typename uint64d<D>::Type matrixSize,
	typename uint64d<D>::Type matrixSizeOs,
	Real kernelWindowSize
){
	this->matrixSize = matrixSize;
	this->matrixSizeOs = matrixSizeOs;
	this->kernelWindowSize = kernelWindowSize;
}

template<class Real, unsigned int D>
Gadgetron::hoNFFT_plan<Real, D>::~hoNFFT_plan(){}

template<class Real, unsigned int D>
void Gadgetron::hoNFFT_plan<Real, D>::preprocess(
	hoNDArray<typename reald<Real,D>::Type> *trajectory,
	NFFT_prep_mode mode
){
}

template<class Real, unsigned int D>
void Gadgetron::hoNFFT_plan<Real, D>::compute(
	hoNDArray<complext<Real>> *in,
	hoNDArray<complext<Real>> *out,
	hoNDArray<Real> *dcw,
	NFFT_computeMode mode
){
	switch(mode){
	case NFFT_BACKWARDS_NC2C:{
		*in *= *dcw;
		compute_NFFTH_NC2C(in, out);
		break;
	}
	case NFFT_FORWARDS_NC2C:
		// todo 
		break;
	case NFFT_BACKWARDS_C2NC:
		// todo
		break;
	case NFFT_FORWARDS_C2NC:
		// todo 
		break;
	}
}


// Convolution for different dimensiosn

template<class Real, unsigned int D>
void Gadgetron::hoNFFT_plan<Real, D>::convolveHelper(
	hoNDArray<complext<Real>> *input,
	hoNDArray<complext<Real>> *output,
	hoNDArray<Real> *dcw,
	NFFT_convolveMode mode
){
	switch(D){
	case 2:{

		break;
	}
	case 3:
		// todo
		break;
	}
}

template<class Real, unsigned int D>
void Gadgetron::hoNFFT_plan<Real, D>::convolve(
	hoNDArray<complext<Real>> *input,
	hoNDArray<complext<Real>> *output,
	hoNDArray<Real> *dcw,
	NFFT_convolveMode mode
){
	switch(mode){
	case NFFT_CONV_NC2C:{
		convolveHelper(input, output, dcw, mode);
		break;
	}
	case NFFT_CONV_C2NC:
		// todo
		break;
	}
}

template<class Real, unsigned int D>
void Gadgetron::hoNFFT_plan<Real, D>::fft(
	hoNDArray<complext<Real>> *data,
	NFFT_fftMode mode,
	bool doScale
){
	typename uint64d<D>::Type _dimsToTransform  = counting_vec<size_t, D>();
	vector<size_t> dimsToTransform = to_std_vector(_dimsToTransform);
	if(mode == NFFT_FORWARDS){
		hoNDFFT<Real>::instance()->fft(data);
	}else{
		hoNDFFT<Real>::instance()->ifft(data);
	}
}

template<class Real, unsigned int D>
void Gadgetron::hoNFFT_plan<Real, D>::deapodize(
	hoNDArray<complext<Real>> *input,
	bool fourierDomain
){
}


template<class Real, unsigned int D>
void Gadgetron::hoNFFT_plan<Real, D>::compute_NFFTH_NC2C(
	hoNDArray<complext<Real>> *input,
	hoNDArray<complext<Real>> *output
){
}

template<class Real, unsigned int D>
void Gadgetron::hoNFFT_plan<Real, D>::prepare(
	hoNDArray<typename reald<Real, D>::Type> traj,
	hoNDArray<complext<Real>> data,
	hoNDArray<Real> weights,
	size_t n,
	Real osf,
	Real wg
){
	kw = wg/osf;
	kosf = std::floor(0.91/(osf*1e-3));
	kwidth = osf*kw/2;

	Real tmp = kw*(osf-0.5);
	beta = M_PI*std::sqrt(tmp*tmp-0.8);
	
	// Compute kernel
	p.create(kosf*kwidth+1);
	for(size_t i = 0; i < kosf*kwidth+1; i++){
		Real om = Real(i)/(kosf*kwidth);
		p[i] = bessi0(beta*std::sqrt(1-om*om));
	}
	Real pConst = p[0];
	for(auto it = p.begin(); it != p.end(); it++)
		*it /= pConst;
	p[kosf*kwidth] = 0;

	// Convert kspace to matrix indices
	nx.create(traj.get_number_of_elements());
	ny.create(traj.get_number_of_elements());
	for(size_t i = 0; i < traj.get_number_of_elements(); i++){
		nx[i] = (n*osf/2)+osf*n*traj[i][0];
		ny[i] = (n*osf/2)+osf*n*traj[i][1];
	}
	
	// Compute deapodiztion filter
	hoNDArray<Real> dax(osf*n);
	for(int i = 0; i < osf*n; i++){
		Real x = (i-osf*n/2)/n;
		Real tmp = M_PI*M_PI*kw*kw*x*x-beta*beta;
		auto sqa = std::sqrt(complex<Real>(tmp, 0));
		dax[i] = (std::sin(sqa)/sqa).real();
	}
	auto daxConst = dax[osf*n/2-1];
	for(auto it = dax.begin(); it != dax.end(); it++)
		*it /= daxConst;
	da.create(osf*n, osf*n);
	for(size_t i = 0; i < osf*n; i++)
		for(size_t j = 0; j < osf*n; j++)
			da[i+j*n*osf] = dax[i]*dax[j];
}
// template instansiation

template class EXPORTCPUNFFT Gadgetron::hoNFFT_plan<float, 2>;

