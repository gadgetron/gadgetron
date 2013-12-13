#include "cuDemonsSolver.h"
#include "setup_grid.h"
#include "gpureg_export.h"
using namespace Gadgetron;

template< class T, unsigned int D> static inline  __device__ void partialDerivs(const T* in, const vector_td<unsigned int,D>& dims, vector_td<unsigned int,D>& co, T * out)
{

	T xi = in[co_to_idx<D>((co+dims)%dims,dims)];
	for (int i = 0; i < D; i++){
		co[i]+=1;
		T dt = in[co_to_idx<D>((co+dims)%dims,dims)];
		out[i] = dt-xi;
		co[i]-=1;
	}
}
/***
 *
 * @param fixed The fixed image
 * @param moving The Moving image
 * @param tot_elemens Total number of elements in fixed (and moving)
 * @param dims Dimensions of the subspace into which the convolution needs to be done
 * @param out Output vector field. Must have same dimensions as fixed and moving + an additional D dimension
 * @param alpha Regularization weight
 * @param beta Small constant added to prevent division by 0.
 */

template<class T, unsigned int D> static __global__ void demons_kernel(T* fixed, T* moving,  size_t tot_elements,const vector_td<unsigned int,D> dims, T* out,T alpha, T beta){

	const int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x;


	if (idx <  tot_elements){

		size_t elements= prod(dims);

		unsigned int batch = idx/elements;
		T * fixed_batch = fixed+elements*batch;
		T * moving_batch = moving+elements*batch;


		T dmov[D];
		T dfix[D];

		vector_td<unsigned int,D> co = idx_to_co<D>(idx, dims);

		partialDerivs(fixed_batch,dims,co,dfix);
		partialDerivs(moving_batch,dims,co,dmov);
		T It = moving_batch[idx]-fixed_batch[idx];

		T gradNorm1 = 0;
		T gradNorm2 = 0;
		for (int i = 0; i < D; i++){
			gradNorm1 += dmov[i]*dmov[i];
			gradNorm2 += dfix[i]*dfix[i];
		}



		for(int i = 0; i < D; i++){
			out[idx+i*tot_elements] = It*(dmov[i]/(gradNorm1+alpha*alpha*It*It+beta)+dfix[i]/(gradNorm2+alpha*alpha*It*It+beta));

		}
	}

}



template<class T, unsigned int D>  boost::shared_ptr<cuNDArray<T> > cuDemonsSolver<T,D>::demonicStep(cuNDArray<T>* fixed,cuNDArray<T>* moving){



	std::vector<unsigned int> dims = *fixed->get_dimensions();
	dims.push_back(D);

	vector_td<unsigned int,D> idims = from_std_vector<unsigned int,D>(dims);
	boost::shared_ptr<cuNDArray<T> > out(new cuNDArray<T>(&dims));
	clear(out.get());

	dim3 gridDim;
	dim3 blockDim;
	setup_grid(fixed->get_number_of_elements(),&blockDim,&gridDim);

	demons_kernel<T,D><<< gridDim,blockDim>>>(fixed->get_data_ptr(),moving->get_data_ptr(),fixed->get_number_of_elements(),idims,out->get_data_ptr(),alpha,beta);


	return out;
}


template class EXPORTGPUREG cuDemonsSolver<float, 1>;
template class EXPORTGPUREG cuDemonsSolver<float, 2>;
template class EXPORTGPUREG cuDemonsSolver<float, 3>;
template class EXPORTGPUREG cuDemonsSolver<float, 4>;

template class EXPORTGPUREG cuDemonsSolver<double, 1>;
template class EXPORTGPUREG cuDemonsSolver<double, 2>;
template class EXPORTGPUREG cuDemonsSolver<double, 3>;
template class EXPORTGPUREG cuDemonsSolver<double, 4>;
