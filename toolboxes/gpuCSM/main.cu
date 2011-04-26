#include "b1_map.hcu"
#include "hoNDArray_fileio.h"
#include "cuNDArray.h"
#include "ndarray_vector_td_utilities.hcu"
#include "NFFT.h"
#include "check_CUDA.h"

#include <cutil.h>
#include <iostream>

using namespace std;

int main( int argc, char** argv) 
{
  hoNDArray<float_complext::Type> host_data = read_nd_array<float_complext::Type>("b1_mapping_data/coil_images.cplx");
  //hoNDArray<float_complex> host_data = read_nd_array<float_complex>("b1_mapping_data/5ch.cplx");
  
  if( host_data.get_number_of_dimensions() != 3 ){
    printf("\nInput data is not three-dimensional (a series of images). Quitting!\n");
    exit(1);
  }
  
  // Hold the image data on the device
  cuNDArray<float_complext::Type> device_data(host_data); // Use float_complex to ensure alignment
  
  // But split into two runs (the test data has 32 coils and the csm estimation eats memory)
  vector<unsigned int> reduced_dims(3);
  reduced_dims[0] = (device_data.get_size(0)>>1);
  reduced_dims[1] = device_data.get_size(1);
  reduced_dims[2] = device_data.get_size(2);
  cuNDArray<float_complext::Type> part_data; part_data.create(reduced_dims);

  bool success;
  uintd2 offset = uintd2(0,0);

  // Get reduces size device data array
  success = cuNDA_crop<float_complext::Type,2>( offset, &device_data, &part_data );

  if( !success ){
    printf("\nerror:\n");
    exit(1);
  }

  unsigned int timer; cutCreateTimer(&timer); double time;
  printf("\nComputing CSM (part one)..."); fflush(stdout);
  cutResetTimer( timer ); cutStartTimer( timer );

  auto_ptr< cuNDArray<float_complext::Type> > csm;

  // Compute CSM
  if( success ) 
    csm = estimate_b1_map<float,2>( &part_data );
  
  cudaThreadSynchronize(); cutStopTimer( timer );
  time = cutGetTimerValue( timer ); printf("done: %.1f ms.", time ); fflush(stdout);

  // Output result

  hoNDArray<float_complext::Type> host_csm = csm->to_host();
  write_nd_array<float_complext::Type>( host_csm, "csm_p1.cplx" );

  // Get reduces size device data array - part 2
  if( success ){
    offset = uintd2(reduced_dims[0],0);
    success = cuNDA_crop<float_complext::Type,2>( offset, &device_data, &part_data );
  }

  printf("\nComputing CSM (part two)..."); fflush(stdout);
  cutResetTimer( timer ); cutStartTimer( timer );

  // Compute CSM
  if( success )
    csm = estimate_b1_map<float,2>( (cuNDArray<float_complext::Type>*) &part_data );
  
  cudaThreadSynchronize(); cutStopTimer( timer );
  time = cutGetTimerValue( timer ); printf("done: %.1f ms.", time ); fflush(stdout);

  // Output result

  host_csm = csm->to_host();
  write_nd_array<float_complext::Type>( host_csm, "csm_p2.cplx" );

  if( !success ){
    printf("\nSome error was encountered...");
  }

  printf("\n", time ); fflush(stdout);
  return 0;
}
