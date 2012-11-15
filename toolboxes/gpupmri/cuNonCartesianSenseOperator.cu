#include "cuNonCartesianSenseOperator.h"
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

template<class REAL, unsigned int D, bool ATOMICS> void
cuNonCartesianSenseOperator<REAL,D,ATOMICS>::mult_M( cuNDArray<_complext>* in, cuNDArray<_complext>* out, bool accumulate )
{
  if( (in->get_number_of_elements() != prod(*this->get_domain_dimensions())) || 
      (out->get_number_of_elements() != prod(*this->get_codomain_dimensions())) ) {
    throw std::runtime_error( "cuNonCartesianSenseOperator::mult_M: dimensions mismatch");
  }

  if( !ready_ ) {
  	throw std::runtime_error( "cuNonCartesianSenseOperator::mult_M: plan has not been set up");

  }


  std::vector<unsigned int> full_dimensions = *this->get_domain_dimensions();
  full_dimensions.push_back(this->ncoils_);
  cuNDArray<_complext> tmp(&full_dimensions, this->device_);

  mult_csm( in, &tmp );

  if( accumulate ){
    throw std::runtime_error( "cuNonCartesianSenseOperator::mult_M: accumulation mode not (yet) supported");
  }
  
  // Forwards NFFT
  plan_->compute( &tmp, out, dcw_.get(), NFFT_plan<REAL,D,ATOMICS>::NFFT_FORWARDS_C2NC );


}

template<class REAL, unsigned int D, bool ATOMICS> void
cuNonCartesianSenseOperator<REAL,D,ATOMICS>::mult_MH( cuNDArray<_complext>* in, cuNDArray<_complext>* out, bool accumulate )
{
  if( (out->get_number_of_elements() != prod(*this->get_domain_dimensions())) || 
      (in->get_number_of_elements() != prod(*this->get_codomain_dimensions())) ) {
  	throw std::runtime_error("cuNonCartesianSenseOperator::mult_MH: dimensions mismatch");
  }

  if( !ready_ ) {
  	throw std::runtime_error("cuNonCartesianSenseOperator::mult_MH: plan has not been set up" );
  }

  std::vector<unsigned int> tmp_dimensions = *this->get_domain_dimensions();
  tmp_dimensions.push_back(this->ncoils_);

  cuNDArray<_complext> tmp(&tmp_dimensions, this->device_);
  
  // Do the NFFT
  plan_->compute( in, &tmp, dcw_.get(), NFFT_plan<REAL,D,ATOMICS>::NFFT_BACKWARDS_NC2C );

  if( !accumulate ){
    this->_set_device();
    out->clear();

    this->_restore_device();
    
  }
  
  mult_csm_conj_sum( &tmp, out );
  
}

template<class REAL, unsigned int D, bool ATOMICS> void
cuNonCartesianSenseOperator<REAL,D,ATOMICS>::setup( _uintd matrix_size, _uintd matrix_size_os, REAL W )
{  
  plan_->setup( matrix_size, matrix_size_os, W );
  
  ready_ = true;
}

template<class REAL, unsigned int D, bool ATOMICS> void
cuNonCartesianSenseOperator<REAL,D,ATOMICS>::preprocess( cuNDArray<_reald> *trajectory ) 
{
  if( !this->csm_.get() || this->csm_->get_number_of_elements() == 0 ) {
    throw std::runtime_error( "cuNonCartesianSenseOperator::set_trajectory : Error setting trajectory, csm not set");
  }
  
  if( !ready_ ) {
  	throw std::runtime_error( "cuNonCartesianSenseOperator::preprocess: plan has not been set up");

  }

  if( trajectory ){

    unsigned int num_frames = trajectory->get_number_of_elements()/trajectory->get_size(0);
    std::vector<unsigned int> tmp_dims;
    tmp_dims = *trajectory->get_dimensions();
    tmp_dims.push_back(this->ncoils_);
    this->set_codomain_dimensions(&tmp_dims);
    
    tmp_dims.clear();
    tmp_dims = *this->csm_->get_dimensions();
    tmp_dims.pop_back();
    tmp_dims.push_back(num_frames);
    this->set_domain_dimensions(&tmp_dims);
    
    plan_->preprocess( trajectory, NFFT_plan<REAL,D,ATOMICS>::NFFT_PREP_ALL );
  }
  else {
    throw std::runtime_error( "cuNonCartesianSenseOperator: cannot set trajectory to 0x0.");

  }
  
}

template<class REAL, unsigned int D, bool ATOMICS> void
cuNonCartesianSenseOperator<REAL,D,ATOMICS>::set_dcw( boost::shared_ptr< cuNDArray<REAL> > dcw ) 
{
  if( dcw->get_device() != this->device_ ){
    this->_set_device();
    dcw_ = boost::shared_ptr< cuNDArray<REAL> >(new cuNDArray<REAL>(*dcw.get()));
    this->_restore_device();
  }
  else    
    dcw_ = dcw;
  
}

//
// Instantiations
//

template class EXPORTGPUPMRI cuNonCartesianSenseOperator<float,2,true>;
template class EXPORTGPUPMRI cuNonCartesianSenseOperator<float,2,false>;

template class EXPORTGPUPMRI cuNonCartesianSenseOperator<float,3,true>;
template class EXPORTGPUPMRI cuNonCartesianSenseOperator<float,3,false>;

template class EXPORTGPUPMRI cuNonCartesianSenseOperator<float,4,true>;
template class EXPORTGPUPMRI cuNonCartesianSenseOperator<float,4,false>;


template class EXPORTGPUPMRI cuNonCartesianSenseOperator<double,2,false>;
template class EXPORTGPUPMRI cuNonCartesianSenseOperator<double,3,false>;
template class EXPORTGPUPMRI cuNonCartesianSenseOperator<double,4,false>;
