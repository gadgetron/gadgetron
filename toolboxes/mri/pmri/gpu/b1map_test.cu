#include "b1_map.h"
#include "hoNDArray_fileio.h"
#include "cuNDArray.h"
#include "ndarray_vector_td_utilities.hcu"
#include "NFFT.h"
#include "check_CUDA.h"

#include <cutil.h>
#include <iostream>

using namespace std;
using namespace Gadgetron;
int main( int argc, char** argv)
{
  hoNDArray<float_complext::Type> host_data =
    read_nd_array<float_complext::Type>("b1_mapping_data/coil_images.cplx");

  //hoNDArray<float_complext::Type> host_data =
  //read_nd_array<float_complext::Type>("b1_mapping_data/5ch.cplx");

  if( host_data.get_number_of_dimensions() != 3 ){
    printf("\nInput data is not three-dimensional (a series of images). Quitting!\n");
    exit(1);
  }

  // Copy the image data to the device
  cuNDArray<float_complext::Type> device_data(host_data);

  unsigned int timer; cutCreateTimer(&timer); double time;
  printf("\nComputing CSM..."); fflush(stdout);
  cutResetTimer( timer ); cutStartTimer( timer );

  // Compute CSM
  boost::shared_ptr< cuNDArray<float_complext::Type> > csm = estimate_b1_map<float,2>( &device_data );

  cudaDeviceSynchronize(); cutStopTimer( timer );
  time = cutGetTimerValue( timer ); printf("done: %.1f ms.", time ); fflush(stdout);

  // Output result

  hoNDArray<float_complext::Type> host_csm = csm->to_host();
  write_nd_array<float_complext::Type>( host_csm, "csm.cplx" );

  printf("\n", time ); fflush(stdout);

  CHECK_FOR_CUDA_ERROR();
  return 0;
}
