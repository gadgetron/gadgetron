#include "cuNFFTOperator.h"
#include "ndarray_vector_td_utilities.h"

template<class REAL, unsigned int D> int 
cuNFFTOperator<REAL,D>::mult_M( cuNDArray<complext<REAL> > *in, cuNDArray<complext<REAL> > *out, bool accumulate )
{
  if( !in || !out ){
    std::cerr << "Error: cuNFFTOperator::mult_M : 0x0 input/output not accepted" << std::endl;
    return -1;
  }

  if( !ready_ ) {
    std::cerr << "Error: cuNFFTOperator::mult_M : plan has not been set up" << std::endl;
    return -1;
  }

  if( dimensionsK_.size() == 0 ){
    std::cerr << "Error: cuNFFTOperator::mult_M : preprocessing has not be performed" << std::endl;
    return -1;
  }

  cuNDArray<complext<REAL> > *tmp_out;

  if( accumulate ){
    tmp_out = new cuNDArray<complext<REAL> >();    
    if( !tmp_out->create(out->get_dimensions().get()) ){
      std::cerr << "Error: cuNFFTOperator::mult_M : memory allocation failed" << std::endl;
      return -1;
    }
  }
  else{
    tmp_out = out;
  }
  
  if( !plan_->compute( in, tmp_out, dcw_.get(), NFFT_plan<REAL,D>::NFFT_FORWARDS_C2NC )) {
    std::cerr << "Error: cuNFFTOperator::mult_M : NFFT failed" << std::endl;
    return -1;
  }

  if( accumulate ){
    if( !cuNDA_axpy< complext<REAL> >( REAL(1), tmp_out, out ) ){
      std::cerr << "Error: cuNFFTOperator::mult_M : NFFT accumulation failed" << std::endl;
      return -1;
    }
    delete tmp_out;
  }
  
  return 0;
}

template<class REAL, unsigned int D> int 
cuNFFTOperator<REAL,D>::mult_MH( cuNDArray<complext<REAL> > *in, cuNDArray<complext<REAL> > *out, bool accumulate )
{
  if( !in || !out ){
    std::cerr << "Error: cuNFFTOperator::mult_MH : 0x0 input/output not accepted" << std::endl;
    return -1;
  }

  if( !ready_ ) {
    std::cerr << "Error: cuNFFTOperator::mult_MH : plan has not been set up" << std::endl;
    return -1;
  }

  if( dimensionsK_.size() == 0 ){
    std::cerr << "Error: cuNFFTOperator::mult_MH : preprocessing has not be performed" << std::endl;
    return -1;
  }

  cuNDArray<complext<REAL> > *tmp_out;

  if( accumulate ){
    tmp_out = new cuNDArray<complext<REAL> >();
    if( !tmp_out->create(out->get_dimensions().get()) ){
      std::cerr << "Error: cuNFFTOperator::mult_MH : memory allocation failed" << std::endl;
      return -1;
    }
  }
  else{
    tmp_out = out;
  }

  if( !plan_->compute( in, tmp_out, dcw_.get(), NFFT_plan<REAL,D>::NFFT_BACKWARDS_NC2C )) {
    std::cerr << "Error: cuNFFTOperator::mult_MH : NFFT failed" << std::endl;
    return -1;
  }

  if( accumulate ){
    if( !cuNDA_axpy< complext<REAL> >( REAL(1), tmp_out, out ) ){
      std::cerr << "Error: cuNFFTOperator::mult_MH : NFFT accumulation failed" << std::endl;
      return -1;
    }
    delete tmp_out;
  }
  
  return 0;
}

template<class REAL, unsigned int D> int 
cuNFFTOperator<REAL,D>::mult_MH_M( cuNDArray<complext<REAL> > *in, cuNDArray<complext<REAL> > *out, bool accumulate )
{
  if( !in || !out ){
    std::cerr << "Error: cuNFFTOperator::mult_MH_M : 0x0 input/output not accepted" << std::endl;
    return -1;
  }
  
  if( !ready_ ) {
    std::cerr << "Error: cuNFFTOperator::mult_MH_M : plan has not been set up" << std::endl;
    return -1;
  }
  
  if( dimensionsK_.size() == 0 ){
    std::cerr << "Error: cuNFFTOperator::mult_MH_M : preprocessing has not be performed" << std::endl;
    return -1;
  }

  cuNDArray<complext<REAL> > *tmp_out;

  if( accumulate ){
    tmp_out = new cuNDArray<complext<REAL> >();
    if( !tmp_out->create(out->get_dimensions().get()) ){
      std::cerr << "Error: cuNFFTOperator::mult_MH_M : memory allocation failed (2)" << std::endl;
      return -1;
    }
  }
  else{
    tmp_out = out;
  }

  if( !plan_->mult_MH_M( in, tmp_out, dcw_.get(), dimensionsK_ )) {
    std::cerr << "Error: cuNFFTOperator::mult_MH_M: NFFT failed" << std::endl;
    return -1;
  }
  
  if( accumulate ){
    if( !cuNDA_axpy< complext<REAL> >( REAL(1), tmp_out, out ) ){
      std::cerr << "Error: cuNFFTOperator::mult_MH_M : NFFT accumulation failed" << std::endl;
      return -1;
    }
    delete tmp_out;
  }
  
  return 0;
}

template<class REAL, unsigned int D> int 
cuNFFTOperator<REAL,D>::setup( typename uintd<D>::Type matrix_size, typename uintd<D>::Type matrix_size_os, REAL W )
{  
  if( !plan_->setup( matrix_size, matrix_size_os, W )) {
    std::cerr << "Error: cuNFFTOperator : failed to setup plan" << std::endl;
    return -1;
  }
  
  ready_ = true;
  return 0;
}

template<class REAL, unsigned int D> int 
cuNFFTOperator<REAL,D>::preprocess( cuNDArray<typename reald<REAL,D>::Type> *trajectory ) 
{
  if( !ready_ ) {
    std::cerr << "Error: cuNFFTOperator::preprocess : plan has not been set up" << std::endl;
    return -1;
  }

  if( trajectory ){

    dimensionsK_.clear();
    dimensionsK_ = *trajectory->get_dimensions();
    
    if( !plan_->preprocess( trajectory, NFFT_plan<REAL,D>::NFFT_PREP_ALL )) {
      std::cerr << "Error: cuNFFTOperator : preprocessing failed" << std::endl;
      return -1;
    }
  }
  else {
    std::cerr << "Error: cuNFFTOperator : cannot set trajectory to 0x0." << std::endl;
    return -1;
  }
  
  return 0;
}


//
// Instantiations
//

template class EXPORTGPUNFFT cuNFFTOperator<float,1>;
template class EXPORTGPUNFFT cuNFFTOperator<float,2>;
template class EXPORTGPUNFFT cuNFFTOperator<float,3>;
template class EXPORTGPUNFFT cuNFFTOperator<float,4>;

template class EXPORTGPUNFFT cuNFFTOperator<double,1>;
template class EXPORTGPUNFFT cuNFFTOperator<double,2>;
template class EXPORTGPUNFFT cuNFFTOperator<double,3>;
template class EXPORTGPUNFFT cuNFFTOperator<double,4>;
