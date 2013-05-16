#include "cuNonCartesianSenseOperator.h"
#include "vector_td_utilities.h"

using namespace Gadgetron;

static unsigned int prodv( std::vector<unsigned int> &vec )
{
  unsigned int result = 1;
  for( unsigned int i=0; i<vec.size(); i++ ){
    result *= vec[i];
  }
  return result;
}

template<class REAL, unsigned int D, bool ATOMICS> void
cuNonCartesianSenseOperator<REAL,D,ATOMICS>::mult_M( cuNDArray< complext<REAL> >* in, cuNDArray< complext<REAL> >* out, bool accumulate )
{
  if( !in || !out ){
    BOOST_THROW_EXCEPTION(runtime_error("cuNonCartesianSenseOperator::mult_M : 0x0 input/output not accepted"));
  }
  
  if( (in->get_number_of_elements() != prodv(*this->get_domain_dimensions())) ||
      (out->get_number_of_elements() != prodv(*this->get_codomain_dimensions())) ) {
    BOOST_THROW_EXCEPTION(runtime_error("cuNonCartesianSenseOperator::mult_M: dimensions mismatch"));
  }

  std::vector<unsigned int> full_dimensions = *this->get_domain_dimensions();
  full_dimensions.push_back(this->ncoils_);
  cuNDArray< complext<REAL> > tmp(&full_dimensions);  
  mult_csm( in, &tmp );
  
  // Forwards NFFT

  if( accumulate ){
    cuNDArray< complext<REAL> > tmp_out(out->get_dimensions());
    plan_->compute( &tmp, &tmp_out, dcw_.get(), cuNFFT_plan<REAL,D,ATOMICS>::NFFT_FORWARDS_C2NC );
    *out += tmp_out;
  }
  else
    plan_->compute( &tmp, out, dcw_.get(), cuNFFT_plan<REAL,D,ATOMICS>::NFFT_FORWARDS_C2NC );
}

template<class REAL, unsigned int D, bool ATOMICS> void
cuNonCartesianSenseOperator<REAL,D,ATOMICS>::mult_MH( cuNDArray< complext<REAL> >* in, cuNDArray< complext<REAL> >* out, bool accumulate )
{
  if( !in || !out ){
    BOOST_THROW_EXCEPTION(runtime_error("cuNonCartesianSenseOperator::mult_MH : 0x0 input/output not accepted"));
  }
  
  if( (out->get_number_of_elements() != prodv(*this->get_domain_dimensions())) ||
      (in->get_number_of_elements() != prodv(*this->get_codomain_dimensions())) ) {
    throw std::runtime_error("cuNonCartesianSenseOperator::mult_MH: dimensions mismatch");
  }

  std::vector<unsigned int> tmp_dimensions = *this->get_domain_dimensions();
  tmp_dimensions.push_back(this->ncoils_);
  cuNDArray< complext<REAL> > tmp(&tmp_dimensions);
  
  // Do the NFFT
  plan_->compute( in, &tmp, dcw_.get(), cuNFFT_plan<REAL,D,ATOMICS>::NFFT_BACKWARDS_NC2C );

  if( !accumulate ){
    clear(out);    
  }
  
  mult_csm_conj_sum( &tmp, out );  
}

template<class REAL, unsigned int D, bool ATOMICS> void
cuNonCartesianSenseOperator<REAL,D,ATOMICS>::setup( _uintd matrix_size, _uintd matrix_size_os, REAL W )
{  
  plan_->setup( matrix_size, matrix_size_os, W );
}

template<class REAL, unsigned int D, bool ATOMICS> void
cuNonCartesianSenseOperator<REAL,D,ATOMICS>::preprocess( cuNDArray<_reald> *trajectory ) 
{
  if( trajectory == 0x0 ){
    BOOST_THROW_EXCEPTION(runtime_error( "cuNonCartesianSenseOperator: cannot preprocess 0x0 trajectory."));
  }
  
  boost::shared_ptr< std::vector<unsigned int> > domain_dims = this->get_domain_dimensions();
  if( domain_dims.get() == 0x0 || domain_dims->size() == 0 ){
    BOOST_THROW_EXCEPTION(runtime_error("cuNonCartesianSenseOperator::preprocess : operator domain dimensions not set"));
  }
  
  {
    std::vector<unsigned int> tmp_dims;
    tmp_dims = *trajectory->get_dimensions();
    tmp_dims.push_back(this->ncoils_);
    this->set_codomain_dimensions(&tmp_dims);
  }
  
  plan_->preprocess( trajectory, cuNFFT_plan<REAL,D,ATOMICS>::NFFT_PREP_ALL );
}

template<class REAL, unsigned int D, bool ATOMICS> void
cuNonCartesianSenseOperator<REAL,D,ATOMICS>::set_dcw( boost::shared_ptr< cuNDArray<REAL> > dcw ) 
{
  dcw_ = dcw;  
}

//
// Instantiations
//

template class EXPORTGPUPMRI cuNonCartesianSenseOperator<float,1,true>;
template class EXPORTGPUPMRI cuNonCartesianSenseOperator<float,1,false>;

template class EXPORTGPUPMRI cuNonCartesianSenseOperator<float,2,true>;
template class EXPORTGPUPMRI cuNonCartesianSenseOperator<float,2,false>;

template class EXPORTGPUPMRI cuNonCartesianSenseOperator<float,3,true>;
template class EXPORTGPUPMRI cuNonCartesianSenseOperator<float,3,false>;

template class EXPORTGPUPMRI cuNonCartesianSenseOperator<float,4,true>;
template class EXPORTGPUPMRI cuNonCartesianSenseOperator<float,4,false>;

template class EXPORTGPUPMRI cuNonCartesianSenseOperator<double,1,false>;
template class EXPORTGPUPMRI cuNonCartesianSenseOperator<double,2,false>;
template class EXPORTGPUPMRI cuNonCartesianSenseOperator<double,3,false>;
template class EXPORTGPUPMRI cuNonCartesianSenseOperator<double,4,false>;
