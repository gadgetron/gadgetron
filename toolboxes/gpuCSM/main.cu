#include "b1_map.hcu"
#include "hoNDArray_fileio.h"
#include "cuNDArray.h"
#include "ndarray_device_utilities.hcu"
#include "NFFT.h"
#include "check_CUDA.h"

#include <cutil.h>
#include <cuComplex.h>
#include <iostream>

using namespace std;

int main( int argc, char** argv) 
{
  hoNDArray<float_complex> host_data = read_nd_array<float_complex>("b1_mapping_data/coil_images.cplx");
  
  if( host_data.get_number_of_dimensions() != 3 ){
    printf("\nInput data is not three-dimensional (a series of images). Quitting!\n");
    exit(1);
  }
  
  cuNDArray<float_complex> device_data(host_data); // Use float_complex to ensure alignment

  unsigned int timer; cutCreateTimer(&timer); double time;
  printf("\nComputing CSM..."); fflush(stdout);
  cutResetTimer( timer ); cutStartTimer( timer );

  // Compute CSM
  auto_ptr< cuNDArray<real_complex<float> > > csm = estimate_b1_map<float,2>( (cuNDArray<real_complex<float> >*) &device_data );
  
  cudaThreadSynchronize(); cutStopTimer( timer );
  time = cutGetTimerValue( timer ); printf("done: %.1f ms.", time ); fflush(stdout);

  //
  // Output result
  //

  hoNDArray<real_complex<float> > host_csm = csm->to_host();
  write_nd_array<real_complex<float> >( host_csm, "csm.cplx" );

  printf("\n", time ); fflush(stdout);
  return 0;
}
