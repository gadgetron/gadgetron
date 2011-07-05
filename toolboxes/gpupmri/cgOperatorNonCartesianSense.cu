#include "cgOperatorNonCartesianSense.h"
#include "vector_td_utilities.h"
#include "ndarray_vector_td_utilities.h"

static unsigned int prod( std::vector<unsigned int> &vec )
{
  unsigned int result = 1;
  for( unsigned int i=0; i<vec.size(); i++ ){
    result *= vec[i];
  }
  return result;
}

template<class REAL, unsigned int D> int 
cgOperatorNonCartesianSense<REAL,D>::mult_M( cuNDArray<_complext>* in, cuNDArray<_complext>* out, bool accumulate )
{
  if( (in->get_number_of_elements() != prod(this->dimensionsI_)) || (out->get_number_of_elements() != prod(this->dimensionsK_)) ) {
    std::cerr << "cgOperatorNonCartesianSense::mult_M dimensions mismatch" << std::endl;
    return -1;
  }

  cuNDArray<_complext> tmp;
  std::vector<unsigned int> full_dimensions = this->dimensionsI_;
  full_dimensions.push_back(this->ncoils_);

  if( !tmp.create(&full_dimensions) ) {
    std::cerr << "cgOperatorNonCartesianSense::mult_M unable to allocate temp array" << std::endl;
    return -2;    
  }
  
  if( mult_csm( in, &tmp ) < 0 ) {
    std::cerr << "cgOperatorNonCartesianSense::mult_M : Unable to multiply with coil sensitivities" << std::endl;
    return -3;
  }

  // Forwards NFFT
  if( !plan_.compute( out, &tmp, dcw_.get(), NFFT_plan<REAL,D>::NFFT_FORWARDS )) {
    std::cerr << "cgOperatorNonCartesianSense::mult_M : failed during NFFT" << std::endl;
    return -4;
  }

  return 0;
}

template<class REAL, unsigned int D> int 
cgOperatorNonCartesianSense<REAL,D>::mult_MH( cuNDArray<_complext>* in, cuNDArray<_complext>* out, bool accumulate )
{
  if( (out->get_number_of_elements() != prod(this->dimensionsI_)) || (in->get_number_of_elements() != prod(this->dimensionsK_)) ) {
    std::cerr << "cgOperatorNonCartesianSense::mult_MH dimensions mismatch" << std::endl;
    return -1;
  }

  std::vector<unsigned int> tmp_dimensions = this->dimensionsI_;
  tmp_dimensions.push_back(this->ncoils_);

  cuNDArray<_complext> tmp;
  if( !tmp.create(&tmp_dimensions) ) {
    std::cerr << "cgOperatorNonCartesianSense::mult_MH: Unable to create temp storage" << std::endl;
    return -2;
  }

  // Do the NFFT
  if( !plan_.compute( in, &tmp, dcw_.get(), NFFT_plan<REAL,D>::NFFT_BACKWARDS )) {
    std::cerr << "cgOperatorNonCartesianSense::mult_MH : failed during NFFT" << std::endl;
    return -3;
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
  if( !plan_.setup( matrix_size, matrix_size_os, W )) {
    std::cerr << "cgOperatorNonCartesianSense: failed to setup plan" << std::endl;
    return -1;
  }
  
  return 0;
}

template<class REAL, unsigned int D> int 
cgOperatorNonCartesianSense<REAL,D>::preprocess( cuNDArray<_reald> *trajectory ) 
{
  if( this->csm_->get_number_of_elements() == 0 ) {
    std::cerr << "cgOperatorNonCartesianSense::set_trajectory : Error setting trajectory, csm not set" << std::endl;
    return -1;
  }
  
  if( trajectory ){

    unsigned int num_frames = trajectory->get_number_of_elements()/trajectory->get_size(0);
    this->dimensionsK_.clear();
    this->dimensionsK_ = *trajectory->get_dimensions();
    this->dimensionsK_.push_back(this->ncoils_);
    
    this->dimensionsI_.clear();
    this->dimensionsI_ = *this->csm_->get_dimensions();
    this->dimensionsI_.pop_back();
    this->dimensionsI_.push_back(num_frames);
    
    if( !plan_.preprocess( trajectory, NFFT_plan<REAL,D>::NFFT_PREP_ALL )) {
      std::cerr << "cgOperatorNonCartesianSense: failed to run preprocess" << std::endl;
      return -2;
    }
  }
  else {
    std::cerr << "cgOperatorNonCartesianSense: cannot set trajectory to 0x0." << std::endl;
    return -3;
  }
  
  return 0;
}

template<class REAL, unsigned int D> int 
cgOperatorNonCartesianSense<REAL,D>::set_dcw( boost::shared_ptr< cuNDArray<REAL> > dcw ) 
{
  dcw_ = dcw;
  
  return 0;
}

//
// Instantiations
//

template class EXPORTGPUPMRI cgOperatorNonCartesianSense<float,2>;
