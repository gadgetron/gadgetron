#pragma once

#include "linearOperator.h"
#include "gpupmri_export.h"

#include <boost/smart_ptr.hpp>
#include <iostream>

template<class REAL, unsigned int D, class ARRAY_TYPE> class EXPORTGPUPMRI senseOperator 
	: public linearOperator<REAL, ARRAY_TYPE>
{

public:

  senseOperator() : linearOperator<REAL, ARRAY_TYPE>(), ncoils_(0) {}
  virtual ~senseOperator() {}

  inline unsigned int get_number_of_coils() { return ncoils_; }
  inline boost::shared_ptr<ARRAY_TYPE> get_csm() { return csm_; }

  virtual int set_csm( boost::shared_ptr<ARRAY_TYPE> csm )
  {
    if( csm.get() && csm->get_number_of_dimensions() == D+1 ) {
      csm_ = csm;      
      ncoils_ = csm_->get_size(D);
      dimensionsI_ = *csm_->get_dimensions();
      dimensionsI_.pop_back();
      return 0;
    }
    else{
      std::cerr << "senseOperator::set_csm: dimensionality mismatch " << std::endl;
      return -1;
    }    
  }

  virtual int mult_M( ARRAY_TYPE* in, ARRAY_TYPE* out, bool accumulate = false ) = 0;
  virtual int mult_MH( ARRAY_TYPE* in, ARRAY_TYPE* out, bool accumulate = false ) = 0;

  virtual int mult_MH_M( ARRAY_TYPE* in, ARRAY_TYPE* out, bool accumulate = false )
  {    
    ARRAY_TYPE tmp;
    if( !tmp.create(&dimensionsK_) ) {
      std::cerr << "senseOperator::mult_MH_M: Unable to create temporary storage" << std::endl;
      return -1;
    }
    
    if( mult_M( in, &tmp, false ) < 0 ) {
      std::cerr << "senseOperator::mult_MH_M: Unable to perform mult_M" << std::endl;
      return -2;
    }
    
    if( mult_MH( &tmp, out, accumulate ) < 0 ) {
      std::cerr << "senseOperator::mult_MH_M: Unable to perform mult_MH" << std::endl;
      return -3;
    }
    
    return 0;
  }
    
  virtual int mult_csm( ARRAY_TYPE* in, ARRAY_TYPE* out ) = 0;
  virtual int mult_csm_conj_sum( ARRAY_TYPE* in, ARRAY_TYPE* out) = 0;

protected:

  unsigned int ncoils_;
  boost::shared_ptr< ARRAY_TYPE > csm_;
  std::vector<unsigned int> dimensionsI_;
  std::vector<unsigned int> dimensionsK_;
};
