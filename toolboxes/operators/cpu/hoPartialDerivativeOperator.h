/** \file hoPartialDerivativeOperator.h
\brief Partial derivative regularization operator, CPU based.
*/

#pragma once

#include "partialDerivativeOperator.h"
#include "hoNDArray_math.h"
#include "vector_td_utilities.h"

#ifdef USE_OMP
#include <omp.h>
#endif

namespace Gadgetron{

    /** \class hoPartialDerivativeOperator
    \brief CPU implementation of device dependent portions of the partialDerivative operator.
    */
    template <class T, unsigned int D> class hoPartialDerivativeOperator
        : public partialDerivativeOperator<D, hoNDArray<T> >
    {
    public:

        hoPartialDerivativeOperator() : 
          partialDerivativeOperator< D, hoNDArray<T> >(0) {}

          hoPartialDerivativeOperator( size_t dimension ) : 
          partialDerivativeOperator<D, hoNDArray<T> >( dimension ) {}

          virtual ~hoPartialDerivativeOperator() {}

          virtual void compute_partial_derivative( typename int64d<D>::Type stride, hoNDArray<T> *in,
              hoNDArray<T> *out, bool accumulate )
          {
              if( !in || !out || in->get_number_of_elements() != out->get_number_of_elements() ){
                  throw std::runtime_error( "hoPartialDerivativeOperator::compute_partial_derivative : array dimensions mismatch.");
              }

              if( in->get_number_of_dimensions() != D || out->get_number_of_dimensions() != D ){
                  throw std::runtime_error("hoPartialDerivativeOperator::compute_partial_derivative : dimensionality mismatch");
              }

              typename int64d<D>::Type dims = vector_td<long long,D>( from_std_vector<size_t,D>( *(in->get_dimensions().get()) ));

#ifdef USE_OMP
#pragma omp parallel for
#endif
      for( long long idx=0; idx<in->get_number_of_elements(); idx++ ) {

                  T valN, valC;

                  typename int64d<D>::Type co = idx_to_co(idx, dims);
                  typename int64d<D>::Type coN = (co+dims+stride)%dims;

                  valN = in->get_data_ptr()[co_to_idx(coN, dims)];
                  valC = in->get_data_ptr()[co_to_idx(co, dims)];

                  T val = valN-valC;

                  if( accumulate )
                      out->get_data_ptr()[idx] += val;
                  else
                      out->get_data_ptr()[idx] = val;
              }
          }

          virtual void compute_second_order_partial_derivative( typename int64d<D>::Type forwards_stride,
              typename int64d<D>::Type adjoint_stride, 
              hoNDArray<T> *in, hoNDArray<T> *out, bool accumulate )
          {
              if( !in || !out || in->get_number_of_elements() != out->get_number_of_elements() ){
                  throw std::runtime_error( "hoPartialDerivativeOperator::compute_second_order_partial_derivative : array dimensions mismatch.");
              }

              if( in->get_number_of_dimensions() != D || out->get_number_of_dimensions() != D ){
                  throw std::runtime_error( "hoPartialDerivativeOperator::compute_second_order_partial_derivative : dimensionality mismatch");
              }

              typename int64d<D>::Type dims = vector_td<long long,D>( from_std_vector<size_t,D>( *(in->get_dimensions().get()) ));

#ifdef USE_OMP
#pragma omp parallel for
#endif
      for( long long idx=0; idx<in->get_number_of_elements(); idx++ ) {

                  T valN1, valN2, valC;

                  typename int64d<D>::Type co = idx_to_co(idx, dims);
                  typename int64d<D>::Type coN1 = (co+dims+forwards_stride)%dims;
                  typename int64d<D>::Type coN2 = (co+dims+adjoint_stride)%dims;

                  valN1 = in->get_data_ptr()[co_to_idx(coN1, dims)];
                  valN2 = in->get_data_ptr()[co_to_idx(coN2, dims)];
                  valC = in->get_data_ptr()[co_to_idx(co, dims)];

                  T val = valC+valC-valN1-valN2;

                  if( accumulate )
                      out->get_data_ptr()[idx] += val;
                  else
                      out->get_data_ptr()[idx] = val;
              }
          }

    };
}
