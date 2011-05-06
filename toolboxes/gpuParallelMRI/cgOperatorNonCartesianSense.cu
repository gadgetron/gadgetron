#include "cgOperatorNonCartesianSense.h"
#include "vector_td_utilities.h"
#include "ndarray_vector_td_utilities.h"
#include <cublas_v2.h>

template<class REAL, unsigned int D> int 
cgOperatorNonCartesianSense<REAL,D>::mult_M( cuNDArray<_complext>* in, cuNDArray<_complext>* out, bool accumulate )
{
  if (!(in->dimensions_equal(this->dimensions_)) || !(out->dimensions_equal(this->dimensions_out_)) ) {
    std::cerr << "cgOperatorNonCartesianSense::mult_M dimensions mismatch" << std::endl;
    return -1;
  }
  
  cuNDArray<_complext> tmp;
  std::vector<unsigned int> full_dimensions = this->dimensions_;
  full_dimensions.push_back(this->ncoils_);

  if( !tmp.create(full_dimensions) ) {
    std::cerr << "cgOperatorNonCartesianSense::mult_M unable to allocate temp array" << std::endl;
    return -2;    
  }
  
  if( mult_csm( in,&tmp ) < 0 ) {
    std::cerr << "cgOperatorNonCartesianSense::mult_M : Unable to multiply with coil sensitivities" << std::endl;
    return -3;
  }

  std::vector<unsigned int> out_dims = out->get_dimensions();
  std::vector<unsigned int> tmp_dims = tmp.get_dimensions();

  if (this->ncoils_ == 1) {
    out->squeeze();
    tmp.squeeze();
  }

  //Do the NFFT
  if( !plan_.compute( out, &tmp, weights_, NFFT_plan<REAL,D>::NFFT_FORWARDS )) {
    std::cerr << "cgOperatorNonCartesianSense::mult_M : failed during NFFT" << std::endl;
    return -4;
  }

  if (this->ncoils_ == 1) {
    out->reshape(out_dims);
  }

  return 0;
}

template<class REAL, unsigned int D> int 
cgOperatorNonCartesianSense<REAL,D>::mult_MH( cuNDArray<_complext>* in, cuNDArray<_complext>* out, bool accumulate )
{
  if( !(out->dimensions_equal(this->dimensions_)) || !(in->dimensions_equal(this->dimensions_out_)) ) {
    std::cerr << "cgOperatorNonCartesianSense::mult_MH dimensions mismatch" << std::endl;
    return -1;
  }

  std::vector<unsigned int> tmp_dimensions = this->dimensions_;
  tmp_dimensions.push_back(this->ncoils_);

  cuNDArray<_complext> tmp;
  if( !tmp.create(tmp_dimensions) ) {
    std::cerr << "cgOperatorNonCartesianSense::mult_MH: Unable to create temp storage" << std::endl;
    return -2;
  }

  //  if( !weights_ ) {
  //  std::cout << "cgOperatorNonCartesianSense::mult_MH : gridding weights are zero, aborting" << std::endl;
  //  return -1;
  // }
  
  std::vector<unsigned int> tmp_dims_in = in->get_dimensions();
  if (this->ncoils_ == 1) {
    in->squeeze();
    tmp.squeeze();
  }

  // Do the NFFT
  if( !plan_.compute( in, &tmp, weights_, NFFT_plan<REAL,D>::NFFT_BACKWARDS )) {
    std::cerr << "cgOperatorNonCartesianSense::mult_MH : failed during NFFT" << std::endl;
    return -3;
  }
  
  if( this->ncoils_ == 1 ) {
    in->reshape(tmp_dims_in);
  }

  if( !accumulate )
    cuNDA_clear<_complext>( out );
  
  if( mult_csm_conj_sum( &tmp, out ) < 0 ) {
    std::cerr << "cgOperatorNonCartesianSense::mult_MH: Unable to multiply with conjugate of sensitivity maps and sum" << std::endl;
    return -4; 
  }
  
  return 0;
}

template<class REAL, unsigned int D> int 
cgOperatorNonCartesianSense<REAL,D>::setup( _uintd matrix_size, _uintd matrix_size_os, REAL W )
{  
  _uintd fixed_dims;
  to_vector_td<unsigned int,D>( fixed_dims, 0 );
  
  if( !plan_.setup( matrix_size, matrix_size_os, fixed_dims, W )) {
    std::cerr << "cgOperatorNonCartesianSense: failed to setup plan" << std::endl;
    return -1;
  }
  
  return 0;
}

template<class REAL, unsigned int D> int 
cgOperatorNonCartesianSense<REAL,D>::set_trajectory( cuNDArray<_reald>* trajectory ) 
{
  // TEMPORARY
  if( trajectory && trajectory->get_number_of_dimensions() > 1 ){
    std::cout << "\nTEMPORARILY the trajectory array must be one-dimensional" << std::endl;
    return -1;
  }

  if (trajectory) {
    trajectory_ = trajectory;
    this->nsamples_ = trajectory->get_number_of_elements();
    this->dimensions_out_.clear();
    this->dimensions_out_.push_back(this->nsamples_);
    this->dimensions_out_.push_back(this->ncoils_);
    
    if( !plan_.preprocess( trajectory_, NFFT_plan<REAL,D>::NFFT_PREP_ALL )) {
      std::cerr << "cgOperatorNonCartesianSense: failed to run preprocess" << std::endl;
      return -1;
    }
  } 
  else {
    return -2;
  }

  return 0;
}

template<class REAL, unsigned int D> int 
cgOperatorNonCartesianSense<REAL,D>::set_weights( cuNDArray<REAL>* w ) 
{
  if (!trajectory_) {
    std::cerr << "cgOperatorNonCartesianSense::set_weights : Error setting weights, trajectory not set" << std::endl;
    return -1;
  }
  if (w->get_number_of_elements() != trajectory_->get_number_of_elements()) {
    std::cerr << "cgOperatorNonCartesianSense::set_weights : weights dimensionality do not match trajectory" << std::endl;
    return -2;
  }

  weights_ = w;

  return 0;
}

//
// Instantiations
//

template class cgOperatorNonCartesianSense<float,2>;
