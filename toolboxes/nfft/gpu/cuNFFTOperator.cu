#include "cuNFFTOperator.h"
#include "cuNDArray_operators.h"
#include "cuNDArray_elemwise.h"

namespace Gadgetron{

  template<class REAL, unsigned int D> void
  cuNFFTOperator<REAL,D>::mult_M( cuNDArray<complext<REAL> > *in, cuNDArray<complext<REAL> > *out, bool accumulate )
  {
    if( !in || !out ){
      BOOST_THROW_EXCEPTION(runtime_error("Error: cuNFFTOperator::mult_M : 0x0 input/output not accepted"));

    }

    if( !ready_ ) {
      BOOST_THROW_EXCEPTION(runtime_error("Error: cuNFFTOperator::mult_M : plan has not been set up"));

    }

    if( dimensionsK_.size() == 0 ){
      BOOST_THROW_EXCEPTION(runtime_error("Error: cuNFFTOperator::mult_M : preprocessing has not be performed"));

    }

    cuNDArray<complext<REAL> > *tmp_out;

    if( accumulate ){
      tmp_out = new cuNDArray<complext<REAL> >(out->get_dimensions());
    }
    else{
      tmp_out = out;
    }
  
    plan_->compute( in, tmp_out, dcw_.get(), cuNFFT_plan<REAL,D>::NFFT_FORWARDS_C2NC );

    if( accumulate ){
      *out += *tmp_out;
      delete tmp_out;
    }
  

  }

  template<class REAL, unsigned int D> void
  cuNFFTOperator<REAL,D>::mult_MH( cuNDArray<complext<REAL> > *in, cuNDArray<complext<REAL> > *out, bool accumulate )
  {
    if( !in || !out ){
      BOOST_THROW_EXCEPTION(runtime_error("Error: cuNFFTOperator::mult_MH : 0x0 input/output not accepted"));

    }

    if( !ready_ ) {
      BOOST_THROW_EXCEPTION(runtime_error("Error: cuNFFTOperator::mult_MH : plan has not been set up"));

    }

    if( dimensionsK_.size() == 0 ){
      BOOST_THROW_EXCEPTION(runtime_error("Error: cuNFFTOperator::mult_MH : preprocessing has not be performed"));

    }

    cuNDArray<complext<REAL> > *tmp_out;

    if( accumulate ){
      tmp_out = new cuNDArray<complext<REAL> >(out->get_dimensions());
    }
    else{
      tmp_out = out;
    }

    plan_->compute( in, tmp_out, dcw_.get(), cuNFFT_plan<REAL,D>::NFFT_BACKWARDS_NC2C );
    if( accumulate ){
      *out += *tmp_out;
      delete tmp_out;
    }
  

  }

  template<class REAL, unsigned int D> void
  cuNFFTOperator<REAL,D>::mult_MH_M( cuNDArray<complext<REAL> > *in, cuNDArray<complext<REAL> > *out, bool accumulate )
  {
    if( !in || !out ){
      BOOST_THROW_EXCEPTION(runtime_error("Error: cuNFFTOperator::mult_MH_M : 0x0 input/output not accepted"));
    }
    
    if( !ready_ ) {
      BOOST_THROW_EXCEPTION(runtime_error("Error: cuNFFTOperator::mult_MH_M : plan has not been set up"));
    }
    
    if( dimensionsK_.size() == 0 ){
      BOOST_THROW_EXCEPTION(runtime_error("Error: cuNFFTOperator::mult_MH_M : preprocessing has not be performed"));
    }
    
    cuNDArray<complext<REAL> > *tmp_out;
    
    if( accumulate ){
      tmp_out = new cuNDArray<complext<REAL> >(out->get_dimensions());
    }
    else{
      tmp_out = out;
    }
    
    plan_->mult_MH_M( in, tmp_out, dcw_.get(), dimensionsK_ );
    
    if( accumulate ){
      *out += *tmp_out;
      delete tmp_out;
    } 
  }
  
  template<class REAL, unsigned int D> void
  cuNFFTOperator<REAL,D>::setup( typename uintd<D>::Type matrix_size, typename uintd<D>::Type matrix_size_os, REAL W )
  {  
    plan_->setup( matrix_size, matrix_size_os, W );  
    ready_ = true;
  }

  template<class REAL, unsigned int D> void
  cuNFFTOperator<REAL,D>::preprocess( cuNDArray<typename reald<REAL,D>::Type> *trajectory ) 
  {
    if( !ready_ ) {
      BOOST_THROW_EXCEPTION(runtime_error("Error: cuNFFTOperator::preprocess : plan has not been set up"));
    }

    if( trajectory ){
      dimensionsK_.clear();
      dimensionsK_ = *trajectory->get_dimensions();    
      plan_->preprocess( trajectory, cuNFFT_plan<REAL,D>::NFFT_PREP_ALL );
    }
    else {
      BOOST_THROW_EXCEPTION(runtime_error("Error: cuNFFTOperator : cannot set trajectory to 0x0."));
    }
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
}
