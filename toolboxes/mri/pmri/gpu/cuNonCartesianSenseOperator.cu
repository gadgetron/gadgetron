#include "cuNonCartesianSenseOperator.h"
#include "vector_td_utilities.h"

using namespace Gadgetron;


template<class REAL,unsigned int D>
cuNonCartesianSenseOperator<REAL,D>::cuNonCartesianSenseOperator(ConvolutionType conv) : cuSenseOperator<REAL,D>() {

    convolutionType = conv;
    is_preprocessed_ = false;
 }

template<class REAL, unsigned int D> void
cuNonCartesianSenseOperator<REAL,D>::mult_M( cuNDArray< complext<REAL> >* in, cuNDArray< complext<REAL> >* out, bool accumulate )
{

  if( !in || !out ){
    throw std::runtime_error("cuNonCartesianSenseOperator::mult_M : 0x0 input/output not accepted");
  }
  if ( !in->dimensions_equal(&this->domain_dims_) || !out->dimensions_equal(&this->codomain_dims_)){
	  throw std::runtime_error("cuNonCartesianSenseOperator::mult_H: input/output arrays do not match specified domain/codomains");
  }

  std::vector<size_t> full_dimensions = *this->get_domain_dimensions();
  full_dimensions.push_back(this->ncoils_);
  cuNDArray< complext<REAL> > tmp(&full_dimensions);  
  this->mult_csm( in, &tmp );
  
  // Forwards NFFT

  if( accumulate ){
    cuNDArray< complext<REAL> > tmp_out(out->get_dimensions());
    plan_->compute( &tmp, tmp_out, dcw_.get(), NFFT_comp_mode::FORWARDS_C2NC );
    *out += tmp_out;
  }
  else
    plan_->compute( tmp, *out, dcw_.get(), NFFT_comp_mode::FORWARDS_C2NC );
}

template<class REAL, unsigned int D> void
cuNonCartesianSenseOperator<REAL,D>::mult_MH( cuNDArray< complext<REAL> >* in, cuNDArray< complext<REAL> >* out, bool accumulate )
{

  if( !in || !out ){
    throw std::runtime_error("cuNonCartesianSenseOperator::mult_MH : 0x0 input/output not accepted");
  }

  if ( !in->dimensions_equal(&this->codomain_dims_) || !out->dimensions_equal(&this->domain_dims_)){
	  throw std::runtime_error("cuNonCartesianSenseOperator::mult_MH: input/output arrays do not match specified domain/codomains");
  }
  std::vector<size_t> tmp_dimensions = *this->get_domain_dimensions();
  std::vector<size_t> tmp_dimensions_data = *this->get_codomain_dimensions();

  
  auto RO = in->get_size(0);
  auto E1E2 = in->get_size(1);
  auto CHA = in->get_size(2);
  
  tmp_dimensions.push_back(this->ncoils_);
  cuNDArray< complext<REAL> > tmp(&tmp_dimensions);

  // Remove channel dimension if the last dimension is the same as the number of coils
   if(tmp_dimensions_data[tmp_dimensions_data.size()-1]==this->ncoils_ && tmp_dimensions_data.size()>2)
   {
     tmp_dimensions_data.pop_back(); 
     cuNDArray<complext<REAL>> temp_ch_data(&tmp_dimensions_data);
     cuNDArray<complext<REAL>> temp_ch_recon(&tmp_dimensions);


      for (int iCHA = 0; iCHA < CHA; iCHA++)
        {
          cudaMemcpy(temp_ch_data.get_data_ptr(),
          in->get_data_ptr() + RO * E1E2 * iCHA,
          RO * E1E2 * sizeof(complext<REAL>), cudaMemcpyDefault);
                  
          plan_->compute(temp_ch_data, temp_ch_recon, dcw_.get(), NFFT_comp_mode::BACKWARDS_NC2C);
          
          cudaMemcpy(tmp.get_data_ptr() + tmp_dimensions[0] * tmp_dimensions[1] * tmp_dimensions[2] * iCHA,
                      temp_ch_recon.get_data_ptr(),
                      tmp_dimensions[0] * tmp_dimensions[1] * tmp_dimensions[2] * sizeof(complext<REAL>), cudaMemcpyDefault);
        }
   }
    else
    {
     //throw std::runtime_error("cuNonCartesianSenseOperator::Last dimension is not the coil dimension");
     plan_->compute(in, tmp, dcw_.get(), NFFT_comp_mode::BACKWARDS_NC2C);

    }

  
  

  if( !accumulate ){
    clear(out);    
  }
  
  this->mult_csm_conj_sum( &tmp, out );
}

template<class REAL, unsigned int D> void
cuNonCartesianSenseOperator<REAL,D>::setup( _uint64d matrix_size, _uint64d matrix_size_os, REAL W )
{  
  plan_ = NFFT<cuNDArray,REAL,D>::make_plan( matrix_size, matrix_size_os, W,convolutionType );
}

template<class REAL, unsigned int D> void
cuNonCartesianSenseOperator<REAL,D>::preprocess( cuNDArray<_reald> *trajectory )
{
  if( trajectory == 0x0 ){
    throw std::runtime_error( "cuNonCartesianSenseOperator: cannot preprocess 0x0 trajectory.");
  }
  
  boost::shared_ptr< std::vector<size_t> > domain_dims = this->get_domain_dimensions();
  if( domain_dims.get() == 0x0 || domain_dims->size() == 0 ){
    throw std::runtime_error("cuNonCartesianSenseOperator::preprocess : operator domain dimensions not set");
  }
  plan_->preprocess( trajectory, NFFT_prep_mode::ALL );
  is_preprocessed_ = true;
}

template<class REAL, unsigned int D> void
cuNonCartesianSenseOperator<REAL,D>::set_dcw( boost::shared_ptr< cuNDArray<REAL> > dcw )
{
  dcw_ = dcw;  
}

//
// Instantiations
//

template class EXPORTGPUPMRI cuNonCartesianSenseOperator<float,1>;
template class EXPORTGPUPMRI cuNonCartesianSenseOperator<float,2>;
template class EXPORTGPUPMRI cuNonCartesianSenseOperator<float,3>;
template class EXPORTGPUPMRI cuNonCartesianSenseOperator<float,4>;

template class EXPORTGPUPMRI cuNonCartesianSenseOperator<double,1>;
template class EXPORTGPUPMRI cuNonCartesianSenseOperator<double,2>;
template class EXPORTGPUPMRI cuNonCartesianSenseOperator<double,3>;
template class EXPORTGPUPMRI cuNonCartesianSenseOperator<double,4>;
