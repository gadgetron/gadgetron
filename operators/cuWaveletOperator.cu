#include "cuWaveletOperator.h"
#include "cuNDArray_operators.h"
#include "cuNDArray_elemwise.h"
#include "cuNDArray_utils.h"
#include "vector_td_utilities.h"
#include "vector_td_io.h"
#include "complext.h"
#include <iostream>
#include "check_CUDA.h"
#include "cudaDeviceManager.h"
#include <stdio.h>
#include <iostream>
using namespace Gadgetron;


// Template Power function
template<unsigned int i, unsigned int j>
struct Pow
{
  enum { Value = i*Pow<i,j-1>::Value};
};

template <unsigned int i>
struct Pow<i,1>
{
  enum { Value = i};
};


template<class T, unsigned int D, unsigned int N>  struct haahr{

	static inline __device__ void loadData(T* elements, T* data, const vector_td<int,D>& dims){
		int offset = 1;
		for (int i = 0; i < N; i++)
			offset *= dims[i];

		haahr<T,D,N-1>::loadData(elements,data,dims);
		haahr<T,D,N-1>::loadData(elements+Pow<2,N>::Value,data+offset,dims);
	}

	static inline __device__ void predict(T* elements){
		haahr<T,D,N-1>::predict(elements);
		haahr<T,D,N-1>::predict(elements+Pow<2,N>::Value);
		for (int i = 0; i < Pow<2,N>::Value; i++){
			T d = (elements[i]+elements[Pow<2,N>::Value+i])*T(0.5);
			T c = (elements[i]-elements[Pow<2,N>::Value+i])*T(2);
			elements[i] = d;
			elements[Pow<2,N>::Value+i] =c;
		}
	}



	static inline __device__ void ipredict(T* elements){
		for (int i = 0; i < Pow<2,N>::Value; i++){
			T d = (elements[i]+elements[Pow<2,N>::Value+i]/T(4));
			T c = (elements[i]-elements[Pow<2,N>::Value+i]/T(4));
			elements[i] = d;
			elements[Pow<2,N>::Value+i] =c;
		}
		haahr<T,D,N-1>::ipredict(elements);
		haahr<T,D,N-1>::ipredict(elements+Pow<2,N>::Value);
	}

	static inline __device__ void saveData(T* elements, T* data, const vector_td<int,D>& dims){
			int offset = 1;
			for (int i = 0; i < N; i++)
				offset *= dims[i];

			haahr<T,D,N-1>::saveData(elements,data,dims);
			haahr<T,D,N-1>::saveData(elements+Pow<2,N>::Value,data+offset,dims);
		}


};

template<class T, unsigned int D> struct haahr<T,D,0>{

	static inline __device__ void loadData(T* elements, T* data, const vector_td<int,D>& dims){
		elements[0] = data[0];
		elements[1] = data[1];
	}

	static inline __device__ void predict(T* elements){
		T d = (elements[0]+elements[1])*T(0.5);
		T c = (elements[0]-elements[1])*T(2);
		elements[0] = d;
		elements[1] = c;

	}

	static inline __device__ void ipredict(T* elements){
			T d = (elements[0]+elements[1]/T(4));
			T c = (elements[0]-elements[1]/T(4));
			elements[0] = d;
			elements[1] = c;

		}
	static inline __device__ void saveData(T* elements, T* data, const vector_td<int,D>& dims){
		data[0] = elements[0];
		data[1] = elements[1];
	}
};


template<class T, unsigned int D> __global__ void waveletKernel(T* in, T* out, const vector_td<int,D> dims){

	T elements[Pow<2,D>::Value];
	int newsize = prod(dims)/Pow<2,D>::Value;
	const vector_td<int,D> dims2 = dims/2;

	const int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x;
	if( idx < newsize ){
		vector_td<int,D> co = idx_to_co<D>(idx,dims2);
		co *= 2;
		haahr<T,D,D-1>::loadData(elements,in+co_to_idx<D>(co,dims),dims);
		haahr<T,D,D-1>::predict(elements);

		for (int i = 0; i < Pow<2,D>::Value; i++){
			out[i*newsize+idx] = elements[i];
		}

	}
}


template<class T, unsigned int D> __global__ void inv_waveletKernel(T* in, T* out, const vector_td<int,D> dims){

	T elements[Pow<2,D>::Value];

	const vector_td<int,D> dims2 = dims/2;
	int oldsize = prod(dims2);
	const int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x;
	if( idx < oldsize ){

		for (int i = 0; i < Pow<2,D>::Value; i++){
			elements[i] = in[i*oldsize+idx] ;
		}
		haahr<T,D,D-1>::ipredict(elements);
		vector_td<int,D> co = idx_to_co<D>(idx,dims2);
		co *= 2;

		haahr<T,D,D-1>::saveData(elements,out+co_to_idx<D>(co,dims),dims);



	}
}

bool isPowerOfTwo (unsigned int x)
{
  return ((x != 0) && ((x & (~x + 1)) == x));
}

template<typename T> T next_power2(T value)
{
    --value;
    for(size_t i = 1; i < sizeof(T) * CHAR_BIT; i*=2)
        value |= value >> i;
    return value+1;
}


template<class T, unsigned int D> void cuWaveletOperator<T,D>::set_domain_dimensions(std::vector<unsigned int>* dims){

	linearOperator<cuNDArray<T> >::set_domain_dimensions(dims);
	/*std::vector<unsigned int> newdims;
	for (int i = 0; i < dims->size(); i++){
		if (isPowerOfTwo(dims->at(i)))
			newdims.push_back(dims->at(i));
		else
			newdims.push_back(next_power2(dims->at(i)));
	}*/

	//BAD BAD BAD IDEAS are going here. This is a dirty dirty hack
	unsigned int m=0;
	for (int i =0; i < dims->size(); i++)
		if (dims->at(i) > m) m = dims->at(i);

	if (!isPowerOfTwo(m))
		m = next_power2(m);

	std::vector<unsigned int> newdims;

	for (int i =0; i < dims->size(); i++)
		newdims.push_back(m);

	linearOperator<cuNDArray<T> >::set_codomain_dimensions(&newdims);

}

template<class T, unsigned int D> void cuWaveletOperator<T,D>::mult_M(cuNDArray<T>* in, cuNDArray<T>* out, bool accumulate ){

	if (! in->dimensions_equal(this->get_domain_dimensions().get()))
		throw std::runtime_error("cuWaveletOperator::mult_M: size of input array does not match operator domain size.");
	if (! out->dimensions_equal(this->get_codomain_dimensions().get()))
			throw std::runtime_error("cuWaveletOperator::mult_M: size of output array does not match operator domain size.");


		boost::shared_ptr<cuNDArray<T> > tmp_in(new cuNDArray<T>(this->get_codomain_dimensions()));

		if (in->dimensions_equal(tmp_in.get()))
			*tmp_in = *in;
		else
			pad<T,D>(in,tmp_in.get());


		cuNDArray<T>* tmp_out;
		if (accumulate)
			tmp_out =new cuNDArray<T>(tmp_in->get_dimensions());
		else
			tmp_out = out;

	  typename intd<D>::Type dims = to_intd( from_std_vector<unsigned int,D>(*(tmp_in->get_dimensions())));

	  typename intd<D>::Type dims2;
	  	  for (int i = 0; i < D; i++) dims2[i] = 2;

	  while(prod(dims) >= Pow<2,D>::Value){
			int threadsPerBlock =std::min(prod(dims)/Pow<2,D>::Value,cudaDeviceManager::Instance()->max_blockdim());
			dim3 dimBlock( threadsPerBlock);
			int totalBlocksPerGrid = std::max(1,prod(dims)/Pow<2,D>::Value/cudaDeviceManager::Instance()->max_blockdim());
			dim3 dimGrid(totalBlocksPerGrid);


			waveletKernel<T,D><<<dimGrid,dimBlock>>>(tmp_in->get_data_ptr(),tmp_out->get_data_ptr(),dims);
			dims /= 2;
			*tmp_in = *tmp_out;
	  }

	  if (dims == dims2){

	  	std::vector<unsigned int> sdim = to_std_vector(to_uintd(dims));
	  	cuNDArray<T> smallArray(&sdim,tmp_in->get_data_ptr());
	  	smallArray.squeeze();
	  	cuNDArray<T> smallTmp(smallArray);
	  	linearOperator<cuNDArray<T> >* smallWave;

	  	switch(smallArray.get_number_of_dimensions()){
	  		case 1:
	  			smallWave = new cuWaveletOperator<T,1>;
	  			break;
	  		case 2:
	  			smallWave = new cuWaveletOperator<T,2>;
					break;
	  		case 3:
					smallWave = new cuWaveletOperator<T,3>;
					break;
	  		default:
	  			throw std::logic_error("cuWaveletOperator::mult_M: Illegal number of input dimensions given");

	  	}

	  	smallWave->mult_M(&smallTmp,&smallArray,false);



	  }
	  if (accumulate){
	  	*out += *tmp_in;
	  	delete tmp_out;
	  }

}


template<class T, unsigned int D> void cuWaveletOperator<T,D>::mult_MH(cuNDArray<T>* in, cuNDArray<T>* out, bool accumulate ){


	if (! out->dimensions_equal(this->get_domain_dimensions().get()))
			throw std::runtime_error("cuWaveletOperator::mult_MH: size of output array does not match operator domain size.");
		if (! in->dimensions_equal(this->get_codomain_dimensions().get()))
				throw std::runtime_error("cuWaveletOperator::mult_MH: size of input array does not match operator domain size.");
		cuNDArray<T>* tmp_out = new cuNDArray<T>(*in);

		cuNDArray<T>* tmp_in = new cuNDArray<T>(*in);


	  const typename intd<D>::Type dims = to_intd( from_std_vector<unsigned int,D>(*(in->get_dimensions())));

	  typename intd<D>::Type dims2;
	  for (int i = 0; i < D; i++) dims2[i] = 2;

	  while(dims2 <= dims){

			int threadsPerBlock =std::min(prod(dims2)/Pow<2,D>::Value,cudaDeviceManager::Instance()->max_blockdim());
			dim3 dimBlock( threadsPerBlock);
			int totalBlocksPerGrid = std::max(1,prod(dims2)/Pow<2,D>::Value/cudaDeviceManager::Instance()->max_blockdim());
			dim3 dimGrid(totalBlocksPerGrid);
			inv_waveletKernel<T,D><<<dimGrid,dimBlock>>>(tmp_in->get_data_ptr(),tmp_out->get_data_ptr(),dims2);
			dims2 *= 2;
			*tmp_in = *tmp_out;
	  }

	  if (!in->dimensions_equal(&this->domain_dims_)){
	  	delete tmp_in;
	  	tmp_in = new cuNDArray<T>(&this->domain_dims_);
	  	vector_td<unsigned int,D> offset;
	  	for (int i = 0; i < D; i++ ) offset[i] = (this->codomain_dims_[i]-this->domain_dims_[i])/2;
	  	crop<T,D>(offset,tmp_out,tmp_in);
	  }

	  if (accumulate){
	  	*out += *tmp_in;
	  } else {
	  	*out = *tmp_in;
	  }
	  delete tmp_in;
	  delete tmp_out;

}

template class  cuWaveletOperator<float,1>;
template class  cuWaveletOperator<float,2>;
template class  cuWaveletOperator<float,3>;
template class  cuWaveletOperator<float,4>;

template class  cuWaveletOperator<double,1>;
template class  cuWaveletOperator<double,2>;
template class  cuWaveletOperator<double,3>;
template class  cuWaveletOperator<double,4>;

template class  cuWaveletOperator<float_complext,1>;
template class  cuWaveletOperator<float_complext,2>;
template class  cuWaveletOperator<float_complext,3>;
template class  cuWaveletOperator<float_complext,4>;

template class  cuWaveletOperator<double_complext,1>;
template class  cuWaveletOperator<double_complext,2>;
template class  cuWaveletOperator<double_complext,3>;
template class  cuWaveletOperator<double_complext,4>;


