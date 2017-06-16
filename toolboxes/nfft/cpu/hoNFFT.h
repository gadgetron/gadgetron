#pragma once 

#include "hoNDArray.h"
#include "cpunfft_export.h"
#include "vector_td.h"
#include "complext.h"

#include <boost/shared_ptr.hpp>

namespace Gadgetron{
	template<class REAL, unsigned int D>
	class EXPORTCPUNFFT hoNFFT_plan{
	public:
		hoNFFT_plan();
		hoNFFT_plan(
			typename uint64d<D>::Type matrix_size,
			typename uint64d<D>::Type matrix_size_os,
			REAL W
		);
		virtual ~hoNFFT_plan();

		void setup(
			typename uint64d<D>::Type matrix_size,
			typename uint64d<D>::Type matrix_size_os,
			REAL W
		);

		enum NFFT_prep_mode{NFFT_PREP_C2NC, NFFT_PREP_NC2C, NFFT_PREP_ALL};
		void preprocess(
			hoNDArray<typename reald<REAL, D>::Type> *trajectory,
			NFFT_prep_mode mode
		);

		enum NFFT_comp_mode{NFFT_FORWARDS_C2NC, NFFT_FORWARDS_NC2C, NFFT_BACKWARDS_C2NC, NFFT_BACKWARDS_NC2C};
		void compute(
			hoNDArray<complext<REAL>> *in,
			hoNDArray<complext<REAL>> *out,
			hoNDArray<REAL> *dcw,
			NFFT_comp_mode mode
		);

		void mult_MH_M(
			hoNDArray<complext<REAL>> *in,
			hoNDArray<complext<REAL>> *out,
			hoNDArray<REAL> *dcw,
			std::vector<size_t> halfway_dims
		);

	public:
		enum NFFT_conv_mode{NFFT_CONV_C2NC, NFFT_CONV_NC2C};
		void convolve(
			hoNDArray<complext<REAL>> *in,
			hoNDArray<complext<REAL>> *out,
			hoNDArray<REAL> *dcw,
			NFFT_conv_mode mode,
			bool accumulate = false
		);

		enum NFFT_fft_mode{NFFT_FORWARDS, NFFT_BACKWARDS};
		void fft(
			hoNDArray<complext<REAL>> *data,
			NFFT_fft_mode mode,
			bool do_scale = true 
		);

		void deapodize(
			hoNDArray<complext<REAL>> *image,
			bool fourier_domain = false
		);

	public:
		inline typename uint64d<D>::Type get_matrix_size(){return matrix_size;}
		inline typename uint64d<D>::Type get_matris_size_os(){return matrix_size_os;}
		inline REAL get_W(){return W;}
		inline bool is_setup(){return initialized;}

	private:
		enum NFFT_components{_NFFT_CONV_C2NC = 1, _NFFT_CONV_NC2C = 2, _NFFT_FFT = 4, _NFFT_DEAPODIZATION = 8};
		void check_consistency(
			hoNDArray<complext<REAL>> *samples,
			hoNDArray<complext<REAL>> *image,
			hoNDArray<REAL> *dcw
		);

		void barebones();
		void compute_beta();
		boost::shared_ptr<hoNDArray<complext<REAL>>> compute_deapodization_filter(
			bool FFTed = false
		);

		void compute_NFFT_C2NC(
			hoNDArray<complext<REAL>> *in,
			hoNDArray<complext<REAL>> *out
		);

		void compute_NFFT_NC2C(
			hoNDArray<complext<REAL>> *in,
			hoNDArray<complext<REAL>> *out
		);

		void compute_NFFTH_NC2C(
			hoNDArray<complext<REAL>> *in,
			hoNDArray<complext<REAL>> *out
		);

		void compute_NFFTH_C2NC(
			hoNDArray<complext<REAL>> *in,
			hoNDArray<complext<REAL>> *out
		);
		
		void convolve_NFFT_C2NC(
			hoNDArray<complext<REAL>> *in,
			hoNDArray<complext<REAL>> *out,
			bool accumulate
		);

		void convolve_NFFT_NC2C(
			hoNDArray<complext<REAL>> *in,
			hoNDArray<complext<REAL>> *out,
			bool accumulate
		);

		void image_wrap(
			hoNDArray<complext<REAL>> *in,
			hoNDArray<complext<REAL>> *out,
			bool accumulate
		);

	private:
		typename uint64d<D>::Type matrix_size;
		typename uint64d<D>::Type matrix_size_os;
		typename uint64d<D>::Type matrix_size_wrap;

		typename reald<REAL, D>::Type alpha;
		typename reald<REAL, D>::Type beta;

		REAL W;
		unsigned int number_of_samples;
		unsigned int number_of_frames;
		
		boost::shared_ptr<hoNDArray<complext<REAL>>> deapodization_filter;
		boost::shared_ptr<hoNDArray<complext<REAL>>> deapodization_filterFFT;

		bool processed_C2NC, processed_NC2C;
		bool initialized;
	};
}
