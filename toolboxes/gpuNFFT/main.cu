#include "hoNDArray_fileio.h"
#include "cuNDArray.h"
#include "ndarray_device_utilities.hcu"
#include "NFFT.hcu"
#include "check_CUDA.h"

#include <cutil.h>
#include <cuComplex.h>
#include <iostream>

using namespace std;

int main( int argc, char** argv) 
{
  hoNDArray<cuFloatComplex> host_data = read_nd_array<cuFloatComplex>("data/raw_data.cplx");
  hoNDArray<float> host_traj = read_nd_array<float>("data/co.real");
  hoNDArray<float> host_weights = read_nd_array<float>("data/weights.real");

  host_data.squeeze();
  host_traj.squeeze();
  host_weights.squeeze();
  
  if( !(host_data.get_number_of_dimensions() == 1 || host_data.get_number_of_dimensions() == 2) || 
      !(host_weights.get_number_of_dimensions() == 1  || host_weights.get_number_of_dimensions() == 2 ) ||
      host_traj.get_number_of_dimensions() != 2 ){
    
    printf("\nInput data is not two-dimensional. Quitting!\n");
    exit(1);
  }
  
  if( host_data.get_size(0) != host_weights.get_size(0) ){
    printf("\nInput data dimensions mismatch between sample data and weights. Quitting!\n");
    exit(1);
  }
  
  if( host_data.get_size(0) != host_traj.get_size(1) || host_traj.get_size(0) != 2 ){
    printf("\nInput trajectory data mismatch. Quitting!\n");
    exit(1);
  }
  
  // Matrix sizes
  const uint2 matrix_size = make_uint2(128,128);
  const uint2 matrix_size_os = make_uint2(256,256);

  // Kernel width
  const float W = 5.5f;

  // No fixed dimensions
  uint2 fixed_dims = make_uint2(0,0);

  // Get trajectory dimensions 'float2' style
  vector<unsigned int> traj_dims = host_traj.get_dimensions();
  traj_dims.erase(traj_dims.begin());
  cuNDArray<float> _traj(host_traj); 
  cuNDArray<float2> traj; traj.create(traj_dims, (float2*)_traj.get_data_ptr());
  
  cuNDArray<cuFloatComplex> data(host_data);
  cuNDArray<float> weights(host_weights);

  unsigned int num_batches = (data.get_number_of_dimensions() == 2) ? data.get_size(1) : 1;
  
  // Setup result image
  vector<unsigned int> image_dims = cuNDA_toVec(matrix_size); 
  if( num_batches > 1 ) image_dims.push_back(num_batches);
  cuNDArray<cuFloatComplex> image; image.create( image_dims );
  
  // Initialize plan
  NFFT_plan<uint2, float2, float, cuFloatComplex> plan( matrix_size, matrix_size_os, fixed_dims, W );

  // Time calls
  unsigned int timer; cutCreateTimer(&timer);
  double time;

  printf("\npreprocess..."); fflush(stdout);
  cutResetTimer( timer ); cutStartTimer( timer );

  // Preprocess
  bool success = plan.preprocess( &traj, false );
    
  cudaThreadSynchronize(); cutStopTimer( timer );
  time = cutGetTimerValue( timer ); printf("done: %.1f ms.", time ); fflush(stdout);

  printf("\nNFFT_H..."); fflush(stdout);
  cutResetTimer( timer ); cutStartTimer( timer );

  // Gridder
  //if( success )
  //success = plan.compute_iteration( &data, &image, &weights, NFFT_plan<uint2, float2, float, cuFloatComplex>::NFFT_BACKWARDS );

  if( success )
    success = plan.compute( &data, &image, &weights, NFFT_plan<uint2, float2, float, cuFloatComplex>::NFFT_BACKWARDS );

  cudaThreadSynchronize(); cutStopTimer( timer );
  time = cutGetTimerValue( timer ); printf("done: %.1f ms.", time ); fflush(stdout);

  printf("\nNFFT_iteration..."); fflush(stdout);
  cutResetTimer( timer ); cutStartTimer( timer );

  /*
  if( success )
    success = plan.compute_iteration( &data, &image, &weights, NFFT_plan<uint2, float2, float, cuFloatComplex>::NFFT_FORWARDS );
  */

  cudaThreadSynchronize(); cutStopTimer( timer );
  time = cutGetTimerValue( timer ); printf("done: %.1f ms", time ); fflush(stdout);
  
  CHECK_FOR_CUDA_ERROR();

  if( !success ){
    printf("\nNFFT failed. Quitting.\n");
    exit(1);
  }

  //
  // Output result
  //

  hoNDArray<cuFloatComplex> host_image = image.to_host();
  //hoNDArray<float> host_norm = cuNDA_norm<float, cuFloatComplex>(&image)->to_host();

  write_nd_array<cuFloatComplex>( host_image, "result.cplx" );
  //write_nd_array<float>( host_norm, "result.real" );

  if( num_batches > 1 ) {
    //hoNDArray<float> host_rss = cuNDA_rss<float, cuFloatComplex>(&image, 2)->to_host();
    //write_nd_array<float>( host_rss, "result_rss.real" );
  }

  printf("\n", time ); fflush(stdout);
  return 0;
}
