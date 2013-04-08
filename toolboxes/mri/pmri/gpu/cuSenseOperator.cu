#include "cuSenseOperator.h"
#include "sense_utilities.h"
#include "vector_td_utilities.h"

namespace Gadgetron{

  template<class REAL, unsigned int D> void
  cuSenseOperator<REAL,D>::mult_csm( cuNDArray<complext<REAL> >* in, cuNDArray<complext<REAL> >* out )
  {  
    csm_mult_M<REAL,D>( in, out, this->csm_.get() );
  }
  
  template<class REAL, unsigned int D> void
  cuSenseOperator<REAL,D>::mult_csm_conj_sum( cuNDArray<complext<REAL> > *in, cuNDArray<complext<REAL> > *out )
  {
    csm_mult_MH<REAL,D>( in, out, this->csm_.get() );
  }
  
  //
  // Instantiations
  //
  
  template class EXPORTGPUPMRI cuSenseOperator<float,1>;
  template class EXPORTGPUPMRI cuSenseOperator<float,2>;
  template class EXPORTGPUPMRI cuSenseOperator<float,3>;
  template class EXPORTGPUPMRI cuSenseOperator<float,4>;

  template class EXPORTGPUPMRI cuSenseOperator<double,1>;
  template class EXPORTGPUPMRI cuSenseOperator<double,2>;
  template class EXPORTGPUPMRI cuSenseOperator<double,3>;
  template class EXPORTGPUPMRI cuSenseOperator<double,4>;
}
