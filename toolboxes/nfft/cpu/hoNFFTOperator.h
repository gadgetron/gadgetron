#pragma once 

#include "linearOperator.h"
#include "hoNFFT.h"
#include "cpunfft_export.h"

namespace Gadgetron{
	template<class Real, unsigned int D>
	class EXPORTCPUNFFT hoNFFTOperator :
	public virtual linearOperator<hoNDArray<complext<Real>>>{
		public:
			hoNFFTOperator(): linearOperator<hoNDArray<complext<Real>>>(){}

			virtual ~hoNFFTOperator(){}
			
			virtual void mult_M(
				hoNDArray<complext<Real>> *in,
				hoNDArray<complext<Real>> *out,
				bool accumulate = false
			);

			virtual void mult_MH(
				hoNDArray<complext<Real>> *in,
				hoNDArray<complext<Real>> *out,
				bool accumulate = false
			);

			virtual void mult_MH_M(
				hoNDArray<complext<Real>> *in,
				hoNDArray<complext<Real>> *out,
				bool accumulate = false
			);

			virtual void setup(
				typename uint64d<D>::Type n,
				Real osf,
				Real wg
			);

			virtual void preprocess(
				hoNDArray<typename reald<Real, D>::Type> traj
			);

		protected:
			hoNFFT_plan<Real, D> plan;
			hoNDArray<Real> dcw;
	};
}
