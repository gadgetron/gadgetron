#include "hoNFFTOperator.h"
#include "hoNDArray_elemwise.h"
#include "hoNDArray_reductions.h"
#include "hoNDArray_utils.h"

using namespace std;

namespace Gadgetron{
	template<class Real, unsigned int D>
	void hoNFFTOperator<Real, D>::mult_M(
		hoNDArray<complext<Real>> *in,
		hoNDArray<complext<Real>> *out,
		bool accumulate
	){
		cout << "mult M" << endl;
		plan.compute(*in, *out, dcw, hoNFFT_plan<Real, D>::NFFT_FORWARDS_C2NC);
		cout << "  Norm (in): " << Gadgetron::nrm2(in) << endl;
		cout << "  Norm (out): " << Gadgetron::nrm2(out) << endl;
	}

	template<class Real, unsigned int D>
	void hoNFFTOperator<Real, D>::mult_MH(
		hoNDArray<complext<Real>> *in,
		hoNDArray<complext<Real>> *out,
		bool accumulate
	){
		cout << "mult MH" << endl;
		plan.compute(*in, *out, dcw, hoNFFT_plan<Real, D>::NFFT_BACKWARDS_NC2C);
		cout << "  Norm (in): " << Gadgetron::nrm2(in) << endl;
		cout << "  Norm (out): " << Gadgetron::nrm2(out) << endl;
	}

	template<class Real, unsigned int D>
	void hoNFFTOperator<Real, D>::mult_MH_M(
		hoNDArray<complext<Real>> *in,
		hoNDArray<complext<Real>> *out,
		bool accumulate
	){
		cout << "mult MHM" << endl;
		plan.mult_MH_M(*in, *out);
	}

	template<class Real, unsigned int D>
	void hoNFFTOperator<Real, D>::setup(
		typename uint64d<D>::Type n,
		Real osf,
		Real wg
	){
		plan = hoNFFT_plan<Real, D>(n, osf, wg);
	}

	template<class Real, unsigned int D>
	void hoNFFTOperator<Real, D>::preprocess(
		hoNDArray<typename reald<Real, D>::Type> k
	){
		plan.preprocess(k);
	}
	
	template class EXPORTCPUNFFT hoNFFTOperator<float, 2>;
	template class EXPORTCPUNFFT hoNFFTOperator<double, 2>;
}

