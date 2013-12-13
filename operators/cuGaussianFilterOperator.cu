#include "cuGaussianFilterOperator.h"
#include "setup_grid.h"
#include "gpuoperators_export.h"
#include "vector_td_utilities.h"
#include "math_constants.h"
using namespace Gadgetron;


template<class T, unsigned int D> static __global__ void gauss_kernel(T* in, T* out, const vector_td<unsigned int,D> dims,T sigma,int dim){

	const int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x;

	if (idx < prod(dims)){
		vector_td<unsigned int,D> co = idx_to_co<D>(idx, dims);
		vector_td<unsigned int,D> co2 = co;
		int steps = ceil(sigma*6);
		int length = dims[dim];
		T res = T(0);

		T norm = 0;
		for (int i = -steps; i < steps; i++){
			int index = co[dim]+i;
			T weight = exp(-0.5*T(i*i)/(2*sigma*sigma));
			norm += weight;
			if ((index > 0) && (index <length)){
				co2[dim]=index;

				res += in[co_to_idx<D>(co2,dims)]*weight;

			}
		}

		out[idx] = res/norm;
	}

}








template<class T, unsigned int D> void cuGaussianFilterOperator<T,D>::mult_M(cuNDArray<T>* in,cuNDArray<T>* out, bool accumulate){

	dim3 gridDim;
	dim3 blockDim;
	setup_grid(in->get_number_of_elements(),&blockDim,&gridDim);
	vector_td<unsigned int,D> dims = from_std_vector<unsigned int,D>(*(in->get_dimensions()));
	cuNDArray<T> in2(*in);
	cuNDArray<T> out2(*out);

	std::vector<unsigned int> batch_dim = to_std_vector(dims);
	size_t elements = prod(dims);
	for (int batch =0; batch < in->get_number_of_elements()/elements; batch++){
		cuNDArray<T> in_view2(&batch_dim,in2.get_data_ptr()+elements*batch);
		cuNDArray<T> out_view2(&batch_dim,out2.get_data_ptr()+elements*batch);
		cuNDArray<T> out_view(&batch_dim,out->get_data_ptr()+elements*batch);

		cuNDArray<T>* inptr = &in_view2;
		cuNDArray<T>* outptr = &out_view2;
		for (int d = 0; d < D; d++){
			gauss_kernel<T,D><<<gridDim,blockDim>>>(inptr->get_data_ptr(),outptr->get_data_ptr(),dims,_sigma,d);
			std::swap(inptr,outptr);
		}
		if (accumulate) out_view += *inptr;
		else out_view = *inptr;

	}



}


template<class T, unsigned int D> void cuGaussianFilterOperator<T,D>::mult_MH(cuNDArray<T>* in,cuNDArray<T>* out, bool accumulate){

	this->mult_M(in,out,accumulate);

}


template EXPORTGPUOPERATORS class cuGaussianFilterOperator<float,1>;
template EXPORTGPUOPERATORS class cuGaussianFilterOperator<float,2>;
template EXPORTGPUOPERATORS class cuGaussianFilterOperator<float,3>;
template EXPORTGPUOPERATORS class cuGaussianFilterOperator<float,4>;

template EXPORTGPUOPERATORS class cuGaussianFilterOperator<double,1>;
template EXPORTGPUOPERATORS class cuGaussianFilterOperator<double,2>;
template EXPORTGPUOPERATORS class cuGaussianFilterOperator<double,3>;
template EXPORTGPUOPERATORS class cuGaussianFilterOperator<double,4>;

