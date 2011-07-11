#include "cgOperatorNonCartesianKtSense.h"
#include "cuNDFFT.h"

template<class REAL, unsigned int D> int 
cgOperatorNonCartesianKtSense<REAL,D>::mult_M( cuNDArray<_complext>* in, cuNDArray<_complext>* out, bool accumulate )
{
  if( accumulate ){
    std::cerr << "cgOperatorNonCartesianKtSense::mult_M: accumulation not supported" << std::endl;
    exit(1);
  }

  // Make a copy of the input array as the fft transform in-place and we do not want to alter the input
  cuNDArray<_complext> tmp(*in); // TODO: multi-device support;
  int ret = cuNDFFT<_complext>().fft( &tmp, D ); // TODO: multi-device support;

  if( ret == 0 ){

    ret = cgOperatorNonCartesianSense<REAL,D>::mult_M( &tmp, out, accumulate );

    if( ret < 0 ){
      std::cerr << "cgOperatorNonCartesianKtSense::mult_M: parent failed" << std::endl;
    }
  }
  else{
    std::cerr << "cgOperatorNonCartesianKtSense::mult_M: temporal NDFFT failed" << std::endl;
  }

  return ret;
}

template<class REAL, unsigned int D> int 
cgOperatorNonCartesianKtSense<REAL,D>::mult_MH( cuNDArray<_complext>* in, cuNDArray<_complext>* out, bool accumulate )
{  
  if( accumulate ){
    std::cerr << "cgOperatorNonCartesianKtSense::mult_MH: accumulation not supported" << std::endl;
    exit(1);
  }

  int ret = cgOperatorNonCartesianSense<REAL,D>::mult_MH( in, out, accumulate );

  if( ret == 0 ){
    
    ret = cuNDFFT<_complext>().ifft( out, D ); // TODO: multi-device support
    
    if( ret < 0 ){
      std::cerr << "cgOperatorNonCartesianKtSense::mult_MH: temporal NDFFT failed" << std::endl;
    }
  }
  else{
    std::cerr << "cgOperatorNonCartesianKtSense::mult_M: parent failed" << std::endl;
  }

  return ret;
}


//
// Instantiations
//

template class EXPORTGPUPMRI cgOperatorNonCartesianKtSense<float,2>;
