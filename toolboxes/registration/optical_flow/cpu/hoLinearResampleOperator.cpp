#include "hoLinearResampleOperator.h"
#include "vector_td_utilities.h"
#include "vector_td_operators.h"
#include "hoArmadillo.h"

#include <stdio.h>
#include <cmath>

namespace Gadgetron{

  template <class T, unsigned int D>
  void hoLinearResampleOperator<T,D>::mult_M( hoNDArray<T> *in, hoNDArray<T> *out, bool accumulate )
  {
    if( !this->preprocessed_ ){
      throw std::runtime_error("hoLinearResampleOperator::mult_M(): displacements not set." );
    }
  
    if( !in || !in->get_data_ptr() || !out || !out->get_data_ptr() ){
      throw std::runtime_error("hoLinearResampleOperator::mult_M(): illegal input/output array." );
    }
  
    arma::Row<typename stdType<T>::Type > in_vec = as_arma_row(*in);
    arma::Row<typename stdType<T>::Type > out_vec = as_arma_row(*out);
    out_vec = in_vec*R_T_;
  }

  template <class T, unsigned int D>
  void hoLinearResampleOperator<T,D>::mult_MH( hoNDArray<T> *in, hoNDArray<T> *out, bool accumulate )
  {
    if( !this->preprocessed_ ){
      throw std::runtime_error("hoLinearResampleOperator::mult_M(): displacements not set." );
    }
  
    if( !in || !in->get_data_ptr() || !out || !out->get_data_ptr() ){
      throw std::runtime_error("hoLinearResampleOperator::mult_M(): illegal input/output array." );
    }

    arma::Col<typename stdType<T>::Type > in_vec = as_arma_col(*in);
    arma::Col<typename stdType<T>::Type > out_vec = as_arma_col(*out);
    out_vec = R_T_ * in_vec;
  }
  
  template <class T, unsigned int D>
  void hoLinearResampleOperator<T,D>::reset()
  {
    R_T_.reset();
    resampleOperator< hoNDArray<typename realType<T>::Type>, hoNDArray<T> >::reset();
  }
  
  template <class T, unsigned int D> void
  hoLinearResampleOperator<T,D>::set_displacement_field( boost::shared_ptr< hoNDArray<typename realType<T>::Type> > displacements )
  {
    typedef typename realType<T>::Type REAL;
    this->preprocessed_ = false;

    if( displacements.get() == 0x0 ){
      throw std::runtime_error("hoLinearResampleOperator::set_displacement_field : displacements ptr is 0x0." );
    }  
  
    const int surplus = displacements->get_number_of_dimensions()-D;
  
    if( !( surplus == 1 || surplus == 2 ) ){
      throw std::runtime_error("hoLinearResampleOperator::set_displacement_field : unexpected array dimensionality." );
    }  
  
    // Determine the number of registrations performed
    const unsigned int extended_dim = (surplus == 1) ? 1 : displacements->get_size(D); 
    const unsigned int field_dim = (surplus == 1) ? displacements->get_size(D) : displacements->get_size(D+1);

    if( !(field_dim == D || field_dim == D+1 )){
      throw std::runtime_error("hoLinearResampleOperator::set_displacement_field : illegal tailing array dim" );
    }
  
    const typename uint64d<D>::Type matrix_size = from_std_vector<size_t,D>( displacements->get_dimensions());
    const size_t num_elements_mat = prod(matrix_size);
    const size_t num_elements_ext = prod(matrix_size)*extended_dim;
    
    const unsigned int num_neighbors = this->get_num_neighbors();
    arma::umat locations(2,num_elements_ext*num_neighbors);
    arma::Col<typename realType<T>::Type > values(num_elements_ext*num_neighbors);
    size_t location_index = 0;

    for( size_t idx=0; idx<num_elements_ext; idx++ ){
    
      const size_t batch_no = idx/num_elements_mat;
      const size_t idx_in_batch = idx-batch_no*num_elements_mat;
    
      const typename uint64d<D>::Type co = idx_to_co( idx_in_batch, matrix_size );
      typename reald<REAL,D>::Type co_disp = vector_td<REAL,D>(co);
      for( unsigned int dim=0; dim<D; dim++ ){
        REAL tmp = displacements->get_data_ptr()[dim*num_elements_ext+batch_no*num_elements_mat+idx_in_batch];
        co_disp.vec[dim] += tmp;
      } 
    
      // Determine the number of neighbors
      //
    
      const typename uint64d<D>::Type twos(2);
    
      // Weights are non-zero only if all neighbors exist
      //
    
      if( this->is_border_pixel(co_disp, matrix_size) )
        continue;
    
      // Iterate over all neighbors
      //
    
      size_t mat_j = idx;
      size_t mat_i;
    
      for( unsigned int i=0; i<num_neighbors; i++ ){
      
        // Determine image coordinate of current neighbor
        //
        
        const typename uint64d<D>::Type stride = idx_to_co<size_t,D>( i, twos );

        if( weak_greater_equal( stride, matrix_size ) ) continue; // For dimensions of size 1
        
        typename reald<REAL,D>::Type co_stride;
      
        for( unsigned int dim=0; dim<D; dim++ ){
          if( stride.vec[dim] == 0 ){
            co_stride.vec[dim] = std::floor(co_disp.vec[dim]);
          }
          else{
            co_stride.vec[dim] = std::ceil(co_disp.vec[dim]);
            if( co_stride.vec[dim] == co_disp.vec[dim] )
              co_stride.vec[dim] += REAL(1.0);
          }
        }

        // Validate that the coordinate is within the expected range
        //

        typename uint64d<D>::Type ones(1);
        typename uint64d<D>::Type co_stride_uint64d = vector_td<size_t,D>(co_stride);

        if( weak_greater( co_stride_uint64d, matrix_size-ones ) ){

          for( unsigned int dim=0; dim<D; dim++ ){
            if( co_stride[dim] < REAL(0) )
              co_stride_uint64d[dim] = 0;
            if( co_stride[dim] > (REAL(matrix_size[dim])-REAL(1)) )
              co_stride_uint64d[dim] = matrix_size[dim]-1;
          }
        }
	
        mat_i = co_to_idx(co_stride_uint64d, matrix_size)+batch_no*num_elements_mat;
      
        // Determine weight
        //
      
        REAL weight = REAL(1);
      
        for( unsigned int dim=0; dim<D; dim++ ){	  
          if( stride.vec[dim] == 0 ){
            weight *= (REAL(1.0)-(co_disp.vec[dim]-co_stride.vec[dim])); }
          else{
            weight *= (REAL(1.0)-(co_stride.vec[dim]-co_disp.vec[dim])); }
        }
      
        locations(0,location_index) = mat_i;
        locations(1,location_index) = mat_j;
        values(location_index) = weight;
        location_index++;
      }
    }
    locations.resize(2,location_index);
    values.resize(location_index);
    R_T_ = arma::SpMat<REAL>( locations, values, num_elements_mat*extended_dim, num_elements_ext, false );
    this->preprocessed_ = true;
  }

  template <class T, unsigned int D> bool
  hoLinearResampleOperator<T,D>::is_border_pixel( typename reald<typename realType<T>::Type,D>::Type co, typename uint64d<D>::Type dims )
  {
    typedef typename realType<T>::Type REAL;

    for( unsigned int dim=0; dim<D; dim++ ){
      if( dims[dim] > 1 && ( co[dim] < REAL(0) || co[dim] >= (REAL(dims[dim])-REAL(1)) ) )
        return true;
    }
    return false;
  }

  template <class T, unsigned int D> unsigned int
  hoLinearResampleOperator<T,D>::get_num_neighbors()
  {
    return 1 << D;
  }
  
  template class EXPORTCPUREG hoLinearResampleOperator<float,1>;
  template class EXPORTCPUREG hoLinearResampleOperator<float,2>;
  template class EXPORTCPUREG hoLinearResampleOperator<float,3>;
  template class EXPORTCPUREG hoLinearResampleOperator<float,4>;

  template class EXPORTCPUREG hoLinearResampleOperator<double,1>;
  template class EXPORTCPUREG hoLinearResampleOperator<double,2>;
  template class EXPORTCPUREG hoLinearResampleOperator<double,3>;
  template class EXPORTCPUREG hoLinearResampleOperator<double,4>;
}
