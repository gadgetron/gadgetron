#include "cuKtSenseRHSBuffer.h"
#include "cuNDFFT.h"

template<class REAL, unsigned int D> int 
cuKtSenseRHSBuffer<REAL,D>::mult_MH( cuNDArray<_complext>* in, cuNDArray<_complext>* out, bool accumulate )
{/*
  if( out->get_number_of_dimensions() != D+1 ){
    std::cerr << "cuKtSenseRHSBuffer::mult_MH: unexpected dimensionality of output array: " << std::endl;
    return -1;
  }
  
  int ret = cuSenseRHSBuffer<REAL,D>::mult_MH( in, out, accumulate );
  
  if( ret == 0 )
    return cuNDFFT<_complext>().ifft( out, D ); // TODO: multi-device support
  else
    return ret;*/ return 0;
}

//
// Instantiations
//

template class EXPORTGPUPMRI cuKtSenseRHSBuffer<float,2>;
