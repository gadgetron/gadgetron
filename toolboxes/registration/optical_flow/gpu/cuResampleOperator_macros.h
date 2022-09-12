#pragma once

#include "cuNDArray_operators.h"
#include "cuNDArray_elemwise.h"
#include "cuNDArray_utils.h"
#include "complext.h"

/* 
   This macro definition is a workaround 
   for missing pure virtual device function support in Cuda.
   
   We provide this macro to avoid explicitly duplicating 
   the code below in every "cuResampleOperator-inherited" class.
*/

#define DECLARE_CU_RESAMPLE_OPERATOR_SUPPORT(COMPONENT)                 \
                                                                        \
  template<class T, unsigned int D> __global__ void                     \
  mult_M_kernel_batch( T *in, T *out,                                   \
                       typename realType<T>::Type *displacements,       \
                       typename uintd<D>::Type matrix_size, unsigned int num_batches ) \
  {                                                                     \
    typedef typename realType<T>::Type REAL;                            \
    const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x; \
    const unsigned int num_elements = prod(matrix_size);                \
                                                                        \
    if( idx < num_elements*num_batches ){                               \
                                                                        \
      const unsigned int batch_no = idx/num_elements;                   \
      const unsigned int idx_in_batch = idx-batch_no*num_elements;      \
      const typename uintd<D>::Type co = idx_to_co( idx_in_batch, matrix_size ); \
                                                                        \
      typename reald<REAL,D>::Type co_disp = vector_td<REAL,D>(co); \
      for( unsigned int dim=0; dim<D; dim++ )                           \
        co_disp.vec[dim] +=  displacements[dim*num_elements+idx_in_batch]; \
                                                                        \
      out[idx] = interpolate<T,D>( batch_no, co_disp, matrix_size, in ); \
    }                                                                   \
  }                                                                     \
                                                                        \
  template<class T, unsigned int D> __global__ void                     \
  mult_M_kernel_extended( T *in, T *out,                                \
                          typename realType<T>::Type *displacements,    \
                          typename uintd<D>::Type matrix_size,          \
                          unsigned int num_elements_in,                 \
                          unsigned int extended_size )                  \
  {                                                                     \
    typedef typename realType<T>::Type REAL;                            \
    const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x; \
    const unsigned int num_elements_mat = prod(matrix_size);            \
    const unsigned int num_elements_ext = prod(matrix_size)*extended_size; \
                                                                        \
    if( idx < num_elements_ext ){                                       \
                                                                        \
      const unsigned int batch_no = idx/num_elements_mat;               \
      const unsigned int idx_in_batch = idx-batch_no*num_elements_mat;	\
                                                                        \
      const typename uintd<D>::Type co = idx_to_co( idx_in_batch, matrix_size ); \
                                                                        \
      typename reald<REAL,D>::Type co_disp = vector_td<REAL,D>(co); \
      for( unsigned int dim=0; dim<D; dim++ )                           \
        co_disp.vec[dim] +=  displacements[dim*num_elements_ext+batch_no*num_elements_mat+idx_in_batch]; \
                                                                        \
      out[idx] = interpolate<T,D>( (idx >= num_elements_in) ? 0 : batch_no, \
                                   co_disp, matrix_size, in );          \
    }                                                                   \
  }                                                                     \
                                                                        \
  template<class T, unsigned int D> void                                \
  cu##COMPONENT<T,D>::mult_M( cuNDArray<T> *in, cuNDArray<T> *out, bool accumulate ) \
  {                                                                     \
    if( !in || !out ){                                                  \
      throw cuda_error("cuResampleOperator::mult_M(): illegal input/output array."); \
    }                                                                   \
                                                                        \
    if( !this->offsets_.get() ){                                        \
      throw cuda_error("cuResampleOperator::mult_M(): displacement field not set."); \
    }                                                                   \
                                                                        \
    cuNDArray<T> tmp;                                                   \
    if( accumulate ){                                                   \
      tmp = *out;                                                       \
    }                                                                   \
                                                                        \
    unsigned int num_disp_vectors = this->get_number_of_displacement_vectors(); \
    int surplus = this->offsets_->get_number_of_dimensions()-D;         \
                                                                        \
    if( !( surplus == 1 || surplus == 2) || this->offsets_->get_size(D-1+surplus) < D ){ \
      throw cuda_error("cuResampleOperator::mult_M(): unexpected dimensions of displacement field."); \
    }                                                                   \
                                                                        \
    if( surplus == 1 ){                                                 \
      if( in->get_number_of_elements() != out->get_number_of_elements() ){ \
        throw cuda_error("cuResampleOperator::mult_M(): in/out array dimensions mismatch (1)."); \
      }                                                                 \
      if( (in->get_number_of_elements() % num_disp_vectors ) != 0 ){    \
        throw cuda_error("cuResampleOperator::mult_M(): in/out array dimensions mismatch displacement field."); \
      }                                                                 \
    }                                                                   \
                                                                        \
    if( surplus == 2 ){                                                 \
      if( (out->get_number_of_elements() % in->get_number_of_elements()) != 0 ){ \
        throw cuda_error("cuResampleOperator::mult_M(): in/out array dimensions mismatch (2)."); \
      }                                                                 \
      if( out->get_number_of_dimensions() != (D+1) || out->get_number_of_elements() != num_disp_vectors ){ \
        throw cuda_error("cuResampleOperator::mult_M(): output array dimensions mismatch displacement field."); \
      }                                                                 \
    }                                                                   \
                                                                        \
    typename uint64d<D>::Type matrix_size = from_std_vector<size_t,D>(*in->get_dimensions().get()); \
    unsigned int num_elements_mat = prod(matrix_size);                  \
    unsigned int num_batches = (surplus == 2) ? 1 : in->get_number_of_elements() / num_elements_mat; \
    unsigned int extended_dim = (surplus == 1) ? 1 : out->get_size(D);	\
                                                                        \
    dim3 blockDim, gridDim;                                             \
                                                                        \
    if( surplus == 1 ){                                                 \
      setup_grid( num_elements_mat, &blockDim, &gridDim, num_batches ); \
    }                                                                   \
    else{                                                               \
      setup_grid( num_elements_mat*extended_dim, &blockDim, &gridDim ); \
    }                                                                   \
                                                                        \
    if( surplus == 1 ) {                                                \
      mult_M_kernel_batch<T,D><<< gridDim, blockDim >>>                 \
        ( in->get_data_ptr(), out->get_data_ptr(),                      \
          this->offsets_->get_data_ptr(), vector_td<unsigned int,D>(matrix_size), num_batches ); \
    }                                                                   \
    else{                                                               \
      mult_M_kernel_extended<T,D><<< gridDim, blockDim >>>              \
        ( in->get_data_ptr(), out->get_data_ptr(), this->offsets_->get_data_ptr(), \
          vector_td<unsigned int,D>(matrix_size), in->get_number_of_elements(), extended_dim ); \
    }                                                                   \
                                                                        \
    CHECK_FOR_CUDA_ERROR();                                             \
                                                                        \
    if( accumulate ){                                                   \
      *out += tmp;                                                      \
    }                                                                   \
  }                                                                     \
                                                                        \
  template<class T, unsigned int D> __global__ void                     \
  mult_MH_kernel( T *in, T *out, typename realType<T>::Type *weights,   \
                  unsigned int *indices, unsigned int *lower_bounds, unsigned int *upper_bounds, \
                  unsigned int num_elements, unsigned int num_batches ) \
  {                                                                     \
    typedef typename realType<T>::Type REAL;                            \
    const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x; \
                                                                        \
    if( idx < num_elements*num_batches ){                               \
                                                                        \
      const unsigned int batch_no = idx/num_elements;                   \
      const unsigned int idx_in_batch = idx-batch_no*num_elements;      \
                                                                        \
      const unsigned int lower_bound = lower_bounds[idx_in_batch];      \
      const unsigned int upper_bound = upper_bounds[idx_in_batch];      \
                                                                        \
      T val = T(0);                                                     \
                                                                        \
      if( lower_bound > upper_bound ||                                  \
          lower_bound >= (_get_num_neighbors<D>()*num_elements) ||      \
          upper_bound >= (_get_num_neighbors<D>()*num_elements) ){      \
                                                                        \
        out[idx] = T(0);                                                \
        return;                                                         \
      }                                                                 \
                                                                        \
      for( unsigned int i=lower_bound; i<upper_bound; i++ ){            \
        unsigned int in_idx = indices[i];                               \
        if( in_idx >= num_elements ){                                   \
          val = T(0);                                                   \
          continue;                                                     \
        }                                                               \
        REAL weight = weights[i];                                       \
        val += (in[in_idx+batch_no*num_elements]*weight);               \
      }                                                                 \
      out[idx] = val;                                                   \
    }                                                                   \
  }                                                                     \
                                                                        \
  template<class T, unsigned int D> void                                \
  cu##COMPONENT<T,D>::mult_MH( cuNDArray<T> *in, cuNDArray<T> *out, bool accumulate ) \
  {                                                                     \
    if( !in || !out ){                                                  \
      throw cuda_error("cuResampleOperator::mult_MH(): illegal input/output array."); \
    }                                                                   \
                                                                        \
    if( !this->preprocessed_ ){                                         \
      throw cuda_error("cuResampleOperator::mult_MH(): preprocessing has not been performed."); \
    }                                                                   \
                                                                        \
    cuNDArray<T> tmp;                                                   \
    if( accumulate ){                                                   \
      tmp = *out;                                                       \
    }                                                                   \
                                                                        \
    unsigned int num_disp_vectors = this->get_number_of_displacement_vectors(); \
    int surplus = this->offsets_->get_number_of_dimensions()-D;         \
                                                                        \
    if( surplus == 1 ){                                                 \
      if( in->get_number_of_elements() != out->get_number_of_elements() ){ \
        throw cuda_error("cuResampleOperator::mult_MH(): in/out array dimensions mismatch (1)."); \
      }                                                                 \
      if( (in->get_number_of_elements() % num_disp_vectors ) != 0 ){    \
        throw cuda_error("cuResampleOperator::mult_MH(): in/out array dimensions mismatch displacement field (1)."); \
      }                                                                 \
    }                                                                   \
                                                                        \
    if( surplus == 2 ){                                                 \
      if( (in->get_number_of_elements() % out->get_number_of_elements()) != 0 ){ \
        throw cuda_error("cuResampleOperator::mult_MH(): in/out array dimensions mismatch (2)."); \
      }                                                                 \
      if( in->get_number_of_dimensions() != (D+1) || in->get_number_of_elements() != num_disp_vectors ){ \
        throw cuda_error("cuResampleOperator::mult_MH(): output array dimensions mismatch displacement field."); \
      }                                                                 \
    }                                                                   \
                                                                        \
    cuNDArray<T> *tmp_out = out; bool mod_out = false;                  \
    if( surplus == 2 && (in->get_number_of_elements()/out->get_number_of_elements()) > 1 ){ \
      mod_out = true;                                                   \
      tmp_out = new cuNDArray<T>(in->get_dimensions().get());           \
    }                                                                   \
                                                                        \
    typename uint64d<D>::Type matrix_size = from_std_vector<size_t,D>( *this->offsets_->get_dimensions().get() ); \
    unsigned int num_batches = (surplus == 2) ? 1 : in->get_number_of_elements() / prod(matrix_size); \
    unsigned int extended_dim = (surplus == 1) ? 1 : in->get_size(D);   \
    unsigned int num_elements = prod(matrix_size)*extended_dim;         \
                                                                        \
    dim3 blockDim, gridDim;                                             \
                                                                        \
    setup_grid( num_elements, &blockDim, &gridDim, num_batches );       \
    mult_MH_kernel<T,D><<< gridDim, blockDim >>>                        \
      ( in->get_data_ptr(), tmp_out->get_data_ptr(),                    \
        raw_pointer_cast(&this->weights_[0]), raw_pointer_cast(&this->indices_[0]), \
        raw_pointer_cast(&this->lower_bounds_[0]), raw_pointer_cast(&this->upper_bounds_[0]), \
        num_elements, num_batches );                                    \
                                                                        \
    if( mod_out ){                                                      \
      *out = *sum<T>( tmp_out, D );                                     \
      delete tmp_out;                                                   \
    }                                                                   \
                                                                        \
    CHECK_FOR_CUDA_ERROR();                                             \
                                                                        \
    if( accumulate ){                                                   \
      *out += tmp;                                                      \
    }                                                                   \
  }
