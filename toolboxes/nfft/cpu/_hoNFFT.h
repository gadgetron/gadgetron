#pragma once 

#include "hoNDArray.h"
#include "vector_td.h"
#include "complext.h"
#include "cpunfft_export.h"

#include <boost/shared_ptr.hpp>

template<class Real, unsigned int D>
struct NC2CConvolution;

namespace Gadgetron{
	template<class Real, unsigned int D>
	class EXPORTCPUNFFT hoNFFT_plan{
	public:
		hoNFFT_plan();
		hoNFFT_plan(
			typename uint64d<D>::Type matrixSize,
			typename uint64d<D>::Type matrixSizeOs,
			Real kernelWindowSize
		);
		~hoNFFT_plan();
		hoNDArray<complext<Real>> processAll(
			hoNDArray<typename reald<Real, D>::Type> traj,
			hoNDArray<complext<Real>> data,
			hoNDArray<Real> weights,
			size_t matrixSize,
			Real osf,
			Real kernelWidth
		);
		hoNDArray<complext<Real>> processTesting(
			hoNDArray<typename reald<Real, D>::Type> traj,
			hoNDArray<complext<Real>> data,
			hoNDArray<Real> weights,
			size_t matrixSize,
			Real osf,
			Real kernelWidth
		);

		enum NFFT_prep_mode{
			NFFT_PREP_C2NC,
			NFFT_PREP_NC2C,
			NFFT_PREP_ALL
		};
		void preprocess(
			hoNDArray<typename reald<Real, D>::Type> *trajectory, 
			NFFT_prep_mode mode
		);
		
		enum NFFT_computeMode{
			NFFT_FORWARDS_C2NC,
			NFFT_FORWARDS_NC2C,
			NFFT_BACKWARDS_C2NC,
			NFFT_BACKWARDS_NC2C
		};
		void compute(	
			hoNDArray<complext<Real>> *input,
			hoNDArray<complext<Real>> *output,
			hoNDArray<Real> *densityCompensation,
			NFFT_computeMode mode
		);	

		enum NFFT_convolveMode{
			NFFT_CONV_C2NC,
			NFFT_CONV_NC2C
		};
		void convolve(
			hoNDArray<complext<Real>> *input,
			hoNDArray<complext<Real>> *output,
			hoNDArray<Real> *densityCompensation,
			NFFT_convolveMode mode
		);
		void convolveHelper(
			hoNDArray<complext<Real>> *input,
			hoNDArray<complext<Real>> *output,
			hoNDArray<Real> *dcw,
			NFFT_convolveMode mode
		);

		enum NFFT_fftMode{
			NFFT_FORWARDS,
			NFFT_BACKWARDS
		};
		void fft(
			hoNDArray<complext<Real>> *data,
			NFFT_fftMode mode,
			bool doScale = true
		);

		void deapodize(
			hoNDArray<complext<Real>> *input,
			bool fourierDomain = false
		);
		
		boost::shared_ptr<hoNDArray<complext<Real>>>
		computeDeapodizationFilter(bool fourierDomain);
		
	private:
		void compute_NFFTH_NC2C(
			hoNDArray<complext<Real>> *input,
			hoNDArray<complext<Real>> *output
		);

		void convolve_NFFT_NC2C(
			hoNDArray<complext<Real>> *input,
			hoNDArray<complext<Real>> *output,
			bool accumulate
		);
		
		void prepare(
			hoNDArray<typename reald<Real, D>::Type> traj,
			hoNDArray<complext<Real>> data,
			hoNDArray<Real> weights,
			size_t matrixSize,
			Real osf,
			Real kernelWidth
		);

		typename uint64d<D>::Type matrixSize;
		typename uint64d<D>::Type matrixSizeOs;
		Real kernelWindowSize;
		
		Real kw;
		Real kosf;
		Real kwidth;
		Real halfKw;
		Real beta;

		hoNDArray<Real> p;
		hoNDArray<Real> da;
		hoNDArray<Real> nx;
		hoNDArray<Real> ny;
	};
}
