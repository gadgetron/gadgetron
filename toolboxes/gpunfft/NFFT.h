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

#include "gadgetron_export.h"

#include "vector_td.h"
#include "cuNDArray.h"

#include <thrust/device_vector.h>
#include <boost/shared_ptr.hpp>

template< class REAL, unsigned int D > class EXPORTGPUNFFT NFFT_plan
{
  
public: // Main interface
    
  // Constructors
  NFFT_plan();
  NFFT_plan( typename uintd<D>::Type matrix_size, typename uintd<D>::Type matrix_size_os, REAL W, int device = -1 );

  // Destructor
  virtual ~NFFT_plan();

  // Clear internal storage in plan
  enum NFFT_wipe_mode { NFFT_WIPE_ALL, NFFT_WIPE_PREPROCESSING };
  bool wipe( NFFT_wipe_mode mode );

  // Replan 
  bool setup( typename uintd<D>::Type matrix_size, typename uintd<D>::Type matrix_size_os, REAL W, int device = -1 );
    
  // Preproces trajectory
  enum NFFT_prep_mode { NFFT_PREP_ALL, NFFT_PREP_FORWARDS, NFFT_PREP_BACKWARDS };
  bool preprocess( cuNDArray<typename reald<REAL,D>::Type> *trajectory, NFFT_prep_mode mode  );
    
  // Execute NFFT
  enum NFFT_comp_mode { NFFT_FORWARDS, NFFT_BACKWARDS };
  bool compute( cuNDArray<typename complext<REAL>::Type> *samples, cuNDArray<typename complext<REAL>::Type> *image, cuNDArray<REAL> *weights, NFFT_comp_mode mode );
  bool compute_iteration( cuNDArray<typename complext<REAL>::Type> *samples, cuNDArray<typename complext<REAL>::Type> *image, cuNDArray<REAL> *weights, NFFT_comp_mode mode );
  
public: // Utilities
  
  // NFFT convolution
  bool convolve( cuNDArray<typename complext<REAL>::Type> *samples, cuNDArray<typename complext<REAL>::Type> *image, cuNDArray<REAL> *weights, NFFT_comp_mode mode, bool accumulate = false );
    
  // NFFT FFT
  bool fft( cuNDArray<typename complext<REAL>::Type> *data, NFFT_comp_mode mode, bool do_scale = true );
  
  // NFFT deapodization
  bool deapodize( cuNDArray<typename complext<REAL>::Type> *image );

public: 

  // Custom operators new/delete for windows memory handling across dll boundaries
  void* operator new (size_t bytes) { return ::new char[bytes]; }
  void operator delete (void *ptr) { delete [] static_cast <char *> (ptr); } 
  void * operator new(size_t s, void * p) { return p; }

private:

  enum NFFT_components { NFFT_CONVOLUTION = 1, NFFT_H_CONVOLUTION = 2, NFFT_FFT = 4, NFFT_DEAPODIZATION = 8 };
  bool check_consistency( cuNDArray<typename complext<REAL>::Type> *samples, cuNDArray<typename complext<REAL>::Type> *image, cuNDArray<REAL> *weights, unsigned char components );

  // Shared barebones constructor
  bool barebones();
    
  // Compute beta control parameter for Kaiser-Bessel kernel
  bool compute_beta();

  // Compute deapodization filter
  bool compute_deapodization_filter();

  // A dedicated compute for each of the two NFFT directions
  bool compute_NFFT( cuNDArray<typename complext<REAL>::Type> *samples, cuNDArray<typename complext<REAL>::Type> *image );
  bool compute_NFFT_H( cuNDArray<typename complext<REAL>::Type> *samples, cuNDArray<typename complext<REAL>::Type> *image );

  // A dedicated convolution for each of the two NFFT directions
  bool convolve_NFFT( cuNDArray<typename complext<REAL>::Type> *samples, cuNDArray<typename complext<REAL>::Type> *image, bool accumulate );
  bool convolve_NFFT_H( cuNDArray<typename complext<REAL>::Type> *samples, cuNDArray<typename complext<REAL>::Type> *image, bool accumulate );
   
  // Internal utility to the NFFT_H convolution
  bool image_wrap( cuNDArray<typename complext<REAL>::Type> *source, cuNDArray<typename complext<REAL>::Type> *target, bool accumulate );

private:
    
  typename uintd<D>::Type matrix_size;          // Matrix size
  typename uintd<D>::Type matrix_size_os;       // Oversampled matrix size
  typename uintd<D>::Type matrix_size_wrap;     // Wrap size at border

  REAL alpha;                                   // Oversampling factor
  REAL beta;                                    // Kaiser-Bessel convolution kernel control parameter
  REAL W;                                       // Kernel width in oversampled grid

  unsigned int number_of_samples;               // Number of samples per frame per coil
  unsigned int number_of_frames;                // Number of frames per reconstruction
    
  int device;                                   // Associated device id

  //
  // Internal data structures for convolution and deapodization
  //

  boost::shared_ptr< cuNDArray<typename complext<REAL>::Type> > deapodization_filter; 
   
  thrust::device_vector< typename reald<REAL,D>::Type > *trajectory_positions;
  thrust::device_vector<unsigned int> *tuples_last;
  thrust::device_vector<unsigned int> *bucket_begin, *bucket_end;

  //
  // State variables
  //

  bool preprocessed_NFFT, preprocessed_NFFT_H;
  bool initialized;
};
