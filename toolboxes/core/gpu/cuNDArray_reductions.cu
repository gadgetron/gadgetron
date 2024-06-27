#include "cuNDArray_reductions.h"
#include "setup_grid.h"
#include <thrust/extrema.h>

namespace Gadgetron {

  template<class T> static void 
  find_stride( cuNDArray<T> *in, size_t dim, size_t *stride, std::vector<size_t> *dims )
  {
    *stride = 1;
    for( unsigned int i=0; i<in->get_number_of_dimensions(); i++ ){
      if( i != dim )
        dims->push_back(in->get_size(i));
      if( i < dim )
        *stride *= in->get_size(i);
    }
  }
  
  // Sum
  //
  template<class T> 
  __global__ void sum_kernel( T *in, T *out, 
                              unsigned int stride, unsigned int number_of_batches, unsigned int number_of_elements )
  {
    const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;

    if( idx < number_of_elements ){

      unsigned int in_idx = (idx/stride)*stride*number_of_batches+(idx%stride);

      T val = in[in_idx];

      for( unsigned int i=1; i<number_of_batches; i++ ) 
        val += in[i*stride+in_idx];

      out[idx] = val; 
    }
  }

  // Sum
  //
  template<class T>  boost::shared_ptr< cuNDArray<T> > sum( cuNDArray<T> *in, unsigned int dim )
  {
    // Some validity checks
    if( !(in->get_number_of_dimensions()>1) ){
      throw std::runtime_error("sum: underdimensioned.");;
    }

    if( dim > in->get_number_of_dimensions()-1 ){
      throw std::runtime_error( "sum: dimension out of range.");;
    }

    unsigned int number_of_batches = in->get_size(dim);
    unsigned int number_of_elements = in->get_number_of_elements()/number_of_batches;

    // Setup block/grid dimensions
    dim3 blockDim; dim3 gridDim;
    setup_grid( number_of_elements, &blockDim, &gridDim );

    // Find element stride
    size_t stride; std::vector<size_t> dims;
    find_stride<T>( in, dim, &stride, &dims );

    // Invoke kernel
    boost::shared_ptr< cuNDArray<T> > out(new cuNDArray<T>());
    out->create(dims);

    sum_kernel<T><<< gridDim, blockDim >>>( in->get_data_ptr(), out->get_data_ptr(), stride, number_of_batches, number_of_elements );

    CHECK_FOR_CUDA_ERROR();
    return out;
  }

  template<class T> T mean(cuNDArray<T>* in)
  {
    return thrust::reduce(in->begin(),in->end(),T(0),thrust::plus<T>())/T(in->get_number_of_elements());
  }

  template<class T> T min(cuNDArray<T>* in)
	{
  	return *thrust::min_element(in->begin(),in->end());
	}

  template<class T> T max(cuNDArray<T>* in)
	{
		return *thrust::max_element(in->begin(),in->end());
	}

  template boost::shared_ptr< cuNDArray<float> > sum<float>( cuNDArray<float>*, unsigned int);
  template boost::shared_ptr< cuNDArray<double> > sum<double>( cuNDArray<double>*, unsigned int);
  template boost::shared_ptr< cuNDArray<float_complext> > sum<float_complext>( cuNDArray<float_complext>*, unsigned int);
  template boost::shared_ptr< cuNDArray<double_complext> > sum<double_complext>( cuNDArray<double_complext>*, unsigned int);  

  template float mean<float>(cuNDArray<float>*);
  template float_complext mean<float_complext>(cuNDArray<float_complext>*);
  template double mean<double>(cuNDArray<double>*);
  template double_complext mean<double_complext>(cuNDArray<double_complext>*);

  template float min<float>(cuNDArray<float>*);
  template float max<float>(cuNDArray<float>*);
  template double min<double>(cuNDArray<double>*);
	template double max<double>(cuNDArray<double>*);
}
