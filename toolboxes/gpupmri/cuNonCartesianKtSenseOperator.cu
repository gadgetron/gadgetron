#include "cuNonCartesianKtSenseOperator.h"
#include "cuNDFFT.h"

template<class REAL, unsigned int D> void
cuNonCartesianKtSenseOperator<REAL,D>::mult_M( cuNDArray<_complext>* in, cuNDArray<_complext>* out, bool accumulate )
{
  if( accumulate ){
    throw std::runtime_error( "cuNonCartesianKtSenseOperator::mult_M: accumulation not supported");

  }

  // Make a copy of the input array as the fft transform in-place and we do not want to alter the input
  cuNDArray<_complext> tmp(*in); // TODO: multi-device support;
  cuNDFFT<_complext>().fft( &tmp, D ); // TODO: multi-device support;



	cuNonCartesianSenseOperator<REAL,D>::mult_M( &tmp, out, accumulate );





}

template<class REAL, unsigned int D> void
cuNonCartesianKtSenseOperator<REAL,D>::mult_MH( cuNDArray<_complext>* in, cuNDArray<_complext>* out, bool accumulate )
{  
  if( accumulate ){
    throw std::runtime_error( "cuNonCartesianKtSenseOperator::mult_MH: accumulation not supported");

  }

  cuNonCartesianSenseOperator<REAL,D>::mult_MH( in, out, accumulate );

  cuNDFFT<_complext>().ifft( out, D ); // TODO: multi-device support
    

}

//
// Instantiations
//

template class EXPORTGPUPMRI cuNonCartesianKtSenseOperator<float,2>;
template class EXPORTGPUPMRI cuNonCartesianKtSenseOperator<float,3>;
template class EXPORTGPUPMRI cuNonCartesianKtSenseOperator<float,4>;

template class EXPORTGPUPMRI cuNonCartesianKtSenseOperator<double,2>;
template class EXPORTGPUPMRI cuNonCartesianKtSenseOperator<double,3>;
template class EXPORTGPUPMRI cuNonCartesianKtSenseOperator<double,4>;
