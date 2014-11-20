#include "CSIOperator.h"

using namespace Gadgetron;


template<class T> static __global__ void dft(__restrict__ complext<T>* kspace, const __restrict__ complext<T>* tspace, T* frequencies, unsigned int spiral_length, unsigned int echoes, unsigned int nfreqs,T dte, T dtt){
	const int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x;
	if (idx < spiral_length*nfreqs ){
		complext<T> result = 0;
		T frequency = frequencies[idx/spiral_length];
		T time_offset = dtt*idx%spiral_length;
		unsigned int kpoint = idx%spiral_length;
		for (unsigned int i =0; i < echoes; i++){
			result += exp(complext<T>(0,-frequency*(dte*i+time_offset)))*tspace[kpoint+i*spiral_length];
		}
		tspace[idx] = result;
	}
}

template<class T> static __global__ void dftH(const __restrict__ complext<T>* kspace, __restrict__ complext<T>* tspace, T* frequencies, unsigned int spiral_length, unsigned int echoes, unsigned int nfreqs,T dte, T dtt){
	const int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x;
	if (idx < spiral_length*echoes ){
		complext<T> result = 0;
		unsigned int kpoint = idx%spiral_length;
		T timeshift = dte*idx/spiral_length+dtt*kpoint;

		for (unsigned int i =0; i < nfreqs; i++){
			result += exp(complext<T>(0,frequencies[i])*(timeshift))*kspace[kpoint+i*spiral_length];

		}
		tspace[idx] = result;
	}
}


