/*
  CUDA implementation of the NFFT - standalone gridding example.

  -----------

  A variant of the implementation publsihed in:

  Accelerating the Non-equispaced Fast Fourier Transform on Commodity Graphics Hardware.
  T.S. Sørensen, T. Schaeffter, K.Ø. Noe, M.S. Hansen. 
  IEEE Transactions on Medical Imaging 2008; 27(4):538-547.

  Real-time Reconstruction of Sensitivity Encoded Radial Magnetic Resonance Imaging Using a Graphics Processing Unit.
  T.S. Sørensen, D. Atkinson, T. Schaeffter, M.S. Hansen.
  IEEE Transactions on Medical Imaging 2009; 28(12):1974-1985. 
*/

#include "hoNDArray_fileio.h"
#include "cuNDArray.h"
#include "vector_td_utilities.h"
#include "ndarray_vector_td_utilities.h"
#include "NFFT.h"
#include "check_CUDA.h"

#include <cutil.h>
#include <iostream>

using namespace std;

int main( int argc, char** argv) 
{
  // TODO: use some data from a public sample data repository
  hoNDArray<float_complext::Type> host_data = read_nd_array<float_complext::Type>("data/data__512.cplx");
  hoNDArray<floatd2::Type> host_traj = read_nd_array<floatd2::Type>("data/traj_512.reald");
  hoNDArray<float> host_weights = read_nd_array<float>("data/dcw_512.real");

  if( (host_data.get_number_of_elements() % host_traj.get_number_of_elements()) ||
      (host_traj.get_number_of_elements() != host_weights.get_number_of_elements()) ){
    
    printf("\nInput data not formatted as expected. Quitting!\n");
    exit(1);
  }

  // Number of images to reconstruct
  unsigned int num_batches = host_data.get_number_of_elements() / host_traj.get_number_of_elements();

  // Matrix sizes
  uintd2 matrix_size = uintd2(192,192);
  uintd2 matrix_size_os = uintd2(256,256);

  // Kernel width
  float W = 5.5f;

  // Upload host arrays to device arrays
  cuNDArray<float_complext::Type> data(host_data);
  cuNDArray<floatd2::Type> traj(host_traj);
  cuNDArray<float> weights(host_weights);

  // Setup result image
  vector<unsigned int> image_dims = uintd_to_vector<2>(matrix_size); image_dims.push_back(num_batches);
  cuNDArray<float_complext::Type> image; image.create(image_dims);
  
  // Initialize plan
  NFFT_plan<float, 2> plan( matrix_size, matrix_size_os, W );

  // Time calls
  unsigned int timer; cutCreateTimer(&timer);
  double time;

  printf("\npreprocess..."); fflush(stdout);
  cutResetTimer( timer ); cutStartTimer( timer );

  // Preprocess
  bool success = plan.preprocess( &traj, NFFT_plan<float,2>::NFFT_PREP_ALL );
    
  cudaThreadSynchronize(); cutStopTimer( timer );
  time = cutGetTimerValue( timer ); printf("done: %.1f ms.", time ); fflush(stdout);

  printf("\nNFFT iteration..."); fflush(stdout);
  cutResetTimer( timer ); cutStartTimer( timer );

  // Gridder
  if( success )
    success = plan.compute_iteration( &data, &image, &weights, NFFT_plan<float,2>::NFFT_BACKWARDS );
  
  cudaThreadSynchronize(); cutStopTimer( timer );
  time = cutGetTimerValue( timer ); printf("done: %.1f ms.", time ); fflush(stdout);

  printf("\nNFFT_H..."); fflush(stdout);
  cutResetTimer( timer ); cutStartTimer( timer );

  if( success )
    success = plan.compute( &data, &image, &weights, NFFT_plan<float,2>::NFFT_BACKWARDS );

  cudaThreadSynchronize(); cutStopTimer( timer );
  time = cutGetTimerValue( timer ); printf("done: %.1f ms.", time ); fflush(stdout);

  printf("\nNFFT_iteration..."); fflush(stdout);
  cutResetTimer( timer ); cutStartTimer( timer );
  
  if( success )
    success = plan.compute_iteration( &data, &image, &weights, NFFT_plan<float,2>::NFFT_FORWARDS );
  
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

  hoNDArray<float_complext::Type> host_image = image.to_host();
  hoNDArray<float> host_norm = cuNDA_norm<float,2>(&image)->to_host();

  write_nd_array<float_complext::Type>( host_image, "result.cplx" );
  write_nd_array<float>( host_norm, "result.real" );

   if( num_batches > 1 ) {
     hoNDArray<float> host_rss = cuNDA_rss<float,float_complext::Type>(&image, 2)->to_host();
     write_nd_array<float>( host_rss, "result_rss.real" );
   }

  return 0;
}
