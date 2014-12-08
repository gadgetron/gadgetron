#include "cuNFFTOperator.h"
#include "cuNDArray_operators.h"
#include "cuNDArray_elemwise.h"
#include "cuNDArray_blas.h"

namespace Gadgetron{

  template<class REAL, unsigned int D> void
  cuNFFTOperator<REAL,D>::mult_M( cuNDArray<complext<REAL> > *in, cuNDArray<complext<REAL> > *out, bool accumulate )
  {
    if( !in || !out ){
      throw std::runtime_error("cuNFFTOperator::mult_M : 0x0 input/output not accepted");
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
      throw std::runtime_error("cuNFFTOperator::mult_MH : 0x0 input/output not accepted");
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
      throw std::runtime_error("cuNFFTOperator::mult_MH_M : 0x0 input/output not accepted");
    }
    
    boost::shared_ptr< std::vector<size_t> > codomain_dims = this->get_codomain_dimensions();
    if( codomain_dims.get() == 0x0 || codomain_dims->size() == 0 ){
      throw std::runtime_error("cuNFFTOperator::mult_MH_M : operator codomain dimensions not set");
    }

    cuNDArray<complext<REAL> > *tmp_out;
    
    if( accumulate ){
      tmp_out = new cuNDArray<complext<REAL> >(out->get_dimensions());
    }
    else{
      tmp_out = out;
    }
    
    plan_->mult_MH_M( in, tmp_out, dcw_.get(), *codomain_dims );
    
    if( accumulate ){
      *out += *tmp_out;
      delete tmp_out;
    } 
  }
  
  template<class REAL, unsigned int D> void
  cuNFFTOperator<REAL,D>::setup( typename uint64d<D>::Type matrix_size, typename uint64d<D>::Type matrix_size_os, REAL W )
  {  
    plan_->setup( matrix_size, matrix_size_os, W );  
  }

  template<class REAL, unsigned int D> void
  cuNFFTOperator<REAL,D>::preprocess( cuNDArray<typename reald<REAL,D>::Type> *trajectory ) 
  {
    if( trajectory == 0x0 ){
      throw std::runtime_error("cuNFFTOperator::preprocess : 0x0 trajectory provided.");
    }
    
    plan_->preprocess( trajectory, cuNFFT_plan<REAL,D>::NFFT_PREP_ALL );
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
