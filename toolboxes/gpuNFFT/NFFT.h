/*
  CUDA implementation of the NFFT.

  -----------

  Accelerating the Non-equispaced Fast Fourier Transform on Commodity Graphics Hardware.
  T.S. Sørensen, T. Schaeffter, K.Ø. Noe, M.S. Hansen. 
  IEEE Transactions on Medical Imaging 2008; 27(4):538-547.

  Real-time Reconstruction of Sensitivity Encoded Radial Magnetic Resonance Imaging Using a Graphics Processing Unit.
  T.S. Sørensen, D. Atkinson, T. Schaeffter, M.S. Hansen.
  IEEE Transactions on Medical Imaging 2009; 28(12):1974-1985. 
*/

#pragma once

#include "vectord.h"
#include "cuNDArray.h"
#include <thrust/device_vector.h>

template< class REAL, unsigned int D > class NFFT_plan
{
  
public: // Main interface
    
  // Constructors
  NFFT_plan();
  NFFT_plan( uintd<D> matrix_size, uintd<D> matrix_size_os, uintd<D> fixed_dims, REAL W, int device = -1 );
  NFFT_plan( NFFT_plan<REAL,D> *plan );

  // Destructor
  ~NFFT_plan();

  // Clear internal storage in plan
  enum NFFT_wipe_mode { NFFT_WIPE_ALL, NFFT_WIPE_PREPROCESSING };
  void wipe( NFFT_wipe_mode mode );

  // Replan 
  bool setup( uintd<D> matrix_size, uintd<D> matrix_size_os, uintd<D> fixed_dims, REAL W, int device = -1 );
    
  // Preproces trajectory
  enum NFFT_prep_mode { NFFT_PREP_ALL, NFFT_PREP_FORWARDS, NFFT_PREP_BACKWARDS };
  bool preprocess( cuNDArray< vectord<REAL,D> > *trajectory, NFFT_prep_mode mode  );
    
  // Execute NFFT
  enum NFFT_comp_mode { NFFT_FORWARDS, NFFT_BACKWARDS };
  bool compute( cuNDArray< real_complex<REAL> > *samples, cuNDArray< real_complex<REAL> > *image, cuNDArray<REAL> *weights, NFFT_comp_mode mode );
  bool compute_iteration( cuNDArray< real_complex<REAL> > *samples, cuNDArray< real_complex<REAL> > *image, cuNDArray<REAL> *weights, NFFT_comp_mode mode );
  
public: // Utilities
  
  // NFFT convolution
  bool convolve( cuNDArray< real_complex<REAL> > *samples, cuNDArray< real_complex<REAL> > *image, cuNDArray<REAL> *weights, NFFT_comp_mode mode );
    
  // NFFT FFT
  bool FFT( cuNDArray< real_complex<REAL> > *data, NFFT_comp_mode mode, bool do_scale = true );
  
  // NFFT deapodization
  bool deapodize( cuNDArray< real_complex<REAL> > *image );

private:   

  enum NFFT_components { NFFT_CONVOLUTION = 1, NFFT_H_CONVOLUTION = 2, NFFT_FFT = 4, NFFT_DEAPODIZATION = 8 };
  bool check_consistency( cuNDArray< real_complex<REAL> > *samples, cuNDArray< real_complex<REAL> > *image, cuNDArray<REAL> *weights, unsigned char components );

  // Shared barebones constructor
  bool barebones();
    
  // Compute beta control parameter for Kaiser-Bessel kernel
  bool compute_beta();

  // Compute deapodization filter
  bool compute_deapodization_filter();

  // A dedicated compute for each of the two NFFT directions
  bool compute_NFFT( cuNDArray< real_complex<REAL> > *samples, cuNDArray<  real_complex<REAL> > *image );
  bool compute_NFFT_H( cuNDArray< real_complex<REAL> > *samples, cuNDArray< real_complex<REAL> > *image );

  // A dedicated convolution for each of the two NFFT directions
  bool convolve_NFFT( cuNDArray< real_complex<REAL> > *samples, cuNDArray< real_complex<REAL> > *image );
  bool convolve_NFFT_H( cuNDArray< real_complex<REAL> > *samples, cuNDArray<  real_complex<REAL> > *image );
   
  // Internal utility to the NFFT_H convolution
  bool image_wrap( cuNDArray< real_complex<REAL> > *source, cuNDArray< real_complex<REAL> > *target/*, bool accumulate*/ );

private:
    
  vectord<unsigned int,D> matrix_size;                // Matrix size
  vectord<unsigned int,D> matrix_size_os;             // Oversampled matrix size
  vectord<unsigned int,D> matrix_size_wrap;           // Wrap size at border

  vectord<unsigned int,D> fixed_dims;                 // "Boolean" denoting which dimensions are fixed
  vectord<unsigned int,D> non_fixed_dims;             // "Boolean" denoting which dimensions are non-fixed

  REAL alpha;                          // Oversampling factor
  REAL beta;                           // Kaiser-Bessel convolution kernel control parameter
  REAL W;                              // Kernel width in oversampled grid

  unsigned int number_of_samples;      // Number of amples (per batch)
    
  int device;                          // Associated device id

  //
  // Internal data structures for convolution and deapodization
  //

  cuNDArray< real_complex<REAL> > *deapodization_filter; 
   
  thrust::device_vector< vectord<REAL,D> > *trajectory_positions;
  thrust::device_vector<unsigned int> *tuples_last;
  thrust::device_vector<unsigned int> *bucket_begin, *bucket_end;

  //
  // State variables
  //

  bool preprocessed_NFFT, preprocessed_NFFT_H;
  bool initialized;
};
