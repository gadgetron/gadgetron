#include "cgOperatorKtSenseRHSBuffer.h"
#include "cuNDFFT.h"

template<class REAL, unsigned int D> int 
cgOperatorKtSenseRHSBuffer<REAL,D>::mult_MH( cuNDArray<_complext>* in, cuNDArray<_complext>* out, bool accumulate )
{
  if( out->get_number_of_dimensions() != D+1 ){
    std::cerr << "cgOperatorKtSenseRHSBuffer::mult_MH: unexpected dimensionality of output array: "
	      << D << " " << in->get_number_of_dimensions() << " " << out->get_number_of_dimensions() << std::endl;
    return -1;
  }
  
  int ret = cgOperatorSenseRHSBuffer<REAL,D>::mult_MH( in, out, accumulate );
  
  if( ret == 0 )
    return cuNDFFT<_complext>().ifft( out, D ); // TODO: multi-device support
  else
    return ret;
}

//
// Instantiations
//

template class EXPORTGPUPMRI cgOperatorKtSenseRHSBuffer<float,2>;
