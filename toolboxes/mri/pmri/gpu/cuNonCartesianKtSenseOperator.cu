#include "cuNonCartesianKtSenseOperator.h"
#include "cuNDFFT.h"

using namespace Gadgetron;

template<class REAL, unsigned int D> void
cuNonCartesianKtSenseOperator<REAL,D>::mult_M( cuNDArray< complext<REAL> >* in, cuNDArray< complext<REAL> >* out, bool accumulate )
{
  if( accumulate ){
    throw std::runtime_error( "cuNonCartesianKtSenseOperator::mult_M: accumulation not supported");
  }

  // Make a copy of the input array as the fft transform in-place and we do not want to alter the input
  cuNDArray< complext<REAL> > tmp(*in);
  cuNDFFT<REAL>::instance()->fft( &tmp, D );

  cuNonCartesianSenseOperator<REAL,D>::mult_M( &tmp, out, accumulate );
}

template<class REAL, unsigned int D> void
cuNonCartesianKtSenseOperator<REAL,D>::mult_MH( cuNDArray< complext<REAL> >* in, cuNDArray< complext<REAL> >* out, bool accumulate )
{
  if( accumulate ){
    throw std::runtime_error( "cuNonCartesianKtSenseOperator::mult_MH: accumulation not supported");
  }

  cuNonCartesianSenseOperator<REAL,D>::mult_MH( in, out, accumulate );
  cuNDFFT<REAL>::instance()->ifft( out, D );
}

//
// Instantiations
//

template class cuNonCartesianKtSenseOperator<float,2>;
template class cuNonCartesianKtSenseOperator<float,3>;
template class cuNonCartesianKtSenseOperator<float,4>;

template class cuNonCartesianKtSenseOperator<double,2>;
template class cuNonCartesianKtSenseOperator<double,3>;
template class cuNonCartesianKtSenseOperator<double,4>;
