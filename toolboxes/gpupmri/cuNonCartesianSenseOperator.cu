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

template<class REAL, unsigned int D, bool ATOMICS> int 
cuNonCartesianSenseOperator<REAL,D,ATOMICS>::mult_M( cuNDArray<_complext>* in, cuNDArray<_complext>* out, bool accumulate )
{
  if( (in->get_number_of_elements() != prod(this->dimensionsI_)) || (out->get_number_of_elements() != prod(this->dimensionsK_)) ) {
    std::cerr << "cuNonCartesianSenseOperator::mult_M: dimensions mismatch" << std::endl;
    return -1;
  }

  if( !ready_ ) {
    std::cerr << "cuNonCartesianSenseOperator::mult_M: plan has not been set up" << std::endl;
    return -1;
  }

  cuNDArray<_complext> tmp;
  std::vector<unsigned int> full_dimensions = this->dimensionsI_;
  full_dimensions.push_back(this->ncoils_);

  if( !tmp.create( &full_dimensions, this->device_ ) ) {
    std::cerr << "cuNonCartesianSenseOperator::mult_M: unable to allocate temp array" << std::endl;
    return -2;    
  }
  
  if( mult_csm( in, &tmp ) < 0 ) {
    std::cerr << "cuNonCartesianSenseOperator::mult_M: unable to multiply with coil sensitivities" << std::endl;
    return -3;
  }

  if( accumulate ){
    std::cerr << "cuNonCartesianSenseOperator::mult_M: accumulation mode not (yet) supported" << std::endl;
    return -4; 
  }
  
  // Forwards NFFT
  if( !plan_->compute( &tmp, out, dcw_.get(), NFFT_plan<REAL,D,ATOMICS>::NFFT_FORWARDS_C2NC )) {
    std::cerr << "cuNonCartesianSenseOperator::mult_M : failed during NFFT" << std::endl;
    return -4;
  }

  return 0;
}

template<class REAL, unsigned int D, bool ATOMICS> int 
cuNonCartesianSenseOperator<REAL,D,ATOMICS>::mult_MH( cuNDArray<_complext>* in, cuNDArray<_complext>* out, bool accumulate )
{

  if( (out->get_number_of_elements() != prod(this->dimensionsI_)) || 
      (in->get_number_of_elements() != prod(this->dimensionsK_)) ) {
    std::cerr << "cuNonCartesianSenseOperator::mult_MH: dimensions mismatch" << std::endl;
    return -1;
  }

  if( !ready_ ) {
    std::cerr << "cuNonCartesianSenseOperator::mult_MH: plan has not been set up" << std::endl;
    return -1;
  }

  std::vector<unsigned int> tmp_dimensions = this->dimensionsI_;
  tmp_dimensions.push_back(this->ncoils_);

  cuNDArray<_complext> tmp;
  if( !tmp.create( &tmp_dimensions, this->device_ ) ) {
    std::cerr << "cuNonCartesianSenseOperator::mult_MH: Unable to create temp storage" << std::endl;
    return -2;
  }
  
  // Do the NFFT
  if( !plan_->compute( in, &tmp, dcw_.get(), NFFT_plan<REAL,D,ATOMICS>::NFFT_BACKWARDS_NC2C )) {
    std::cerr << "cuNonCartesianSenseOperator::mult_MH: NFFT failed" << std::endl;
    return -3;
  }

  if( !accumulate ){
    int ret1 = this->_set_device();
    if( ret1 == 0 ) cuNDA_clear<_complext>( out );
    int ret2 = this->_restore_device();
    
    if( ret1<0 || ret2<0 ){
      std::cerr << "cuNonCartesianSenseOperator::mult_MH: device error" << std::endl;
      return -4; 
    }    
  }
  
  if( mult_csm_conj_sum( &tmp, out ) < 0 ) {
    std::cerr << "cuNonCartesianSenseOperator::mult_MH: unable to multiply with conjugate of sensitivity maps and sum" << std::endl;
    return -5; 
  }
  
  return 0;
}

template<class REAL, unsigned int D, bool ATOMICS> int 
cuNonCartesianSenseOperator<REAL,D,ATOMICS>::setup( _uintd matrix_size, _uintd matrix_size_os, REAL W )
{  
  if( !plan_->setup( matrix_size, matrix_size_os, W )) {
    std::cerr << "cuNonCartesianSenseOperator: failed to setup plan" << std::endl;
    return -1;
  }
  
  ready_ = true;
  return 0;
}

template<class REAL, unsigned int D, bool ATOMICS> int 
cuNonCartesianSenseOperator<REAL,D,ATOMICS>::preprocess( cuNDArray<_reald> *trajectory ) 
{
  if( !this->csm_.get() || this->csm_->get_number_of_elements() == 0 ) {
    std::cerr << "cuNonCartesianSenseOperator::set_trajectory : Error setting trajectory, csm not set" << std::endl;
    return -1;
  }
  
  if( !ready_ ) {
    std::cerr << "cuNonCartesianSenseOperator::preprocess: plan has not been set up" << std::endl;
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
    
    if( !plan_->preprocess( trajectory, NFFT_plan<REAL,D,ATOMICS>::NFFT_PREP_ALL )) {
      std::cerr << "cuNonCartesianSenseOperator: failed to run preprocess" << std::endl;
      return -2;
    }
  }
  else {
    std::cerr << "cuNonCartesianSenseOperator: cannot set trajectory to 0x0." << std::endl;
    return -3;
  }
  
  return 0;
}

template<class REAL, unsigned int D, bool ATOMICS> int 
cuNonCartesianSenseOperator<REAL,D,ATOMICS>::set_dcw( boost::shared_ptr< cuNDArray<REAL> > dcw ) 
{
  int ret1 = 0, ret2 = 0;
  
  if( dcw->get_device() != this->device_ ){
    ret1 = this->_set_device();
    if( ret1 == 0 ) dcw_ = boost::shared_ptr< cuNDArray<REAL> >(new cuNDArray<REAL>(*dcw.get()));
    ret2 = this->_restore_device();
  }
  else    
    dcw_ = dcw;
  
  if( ret1 == 0 && ret2 == 0 )
    return 0;
  else{
    std::cout << std::endl << "cuNonCartesianSenseOperator::set_dcw failed" << std::endl;
    return -1;
  }
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

template class EXPORTGPUPMRI cuNonCartesianSenseOperator<double,2,true>;
template class EXPORTGPUPMRI cuNonCartesianSenseOperator<double,2,false>;

template class EXPORTGPUPMRI cuNonCartesianSenseOperator<double,3,true>;
template class EXPORTGPUPMRI cuNonCartesianSenseOperator<double,3,false>;

template class EXPORTGPUPMRI cuNonCartesianSenseOperator<double,4,true>;
template class EXPORTGPUPMRI cuNonCartesianSenseOperator<double,4,false>;
