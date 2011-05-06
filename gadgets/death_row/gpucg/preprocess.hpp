
//
// Generate the necessary data structures for the CUDA implementations of the NFFT and NFFT_H algorithms.
//


#ifndef _PREPROCESS_HPP_
#define _PREPROCESS_HPP_

#include "types.hcu"

namespace mr_recon
{

  // Class delarations

  template< class UINTd, class FLOATd, char TYPE > class NFFT_plan;
  template< class UINTd, class FLOATd, char TYPE > class NFFT_H_plan;
  template< class UINTd, class FLOATd, char TYPE > class NFFT_iteration_plan;


  /*
    TYPE definition:
    ----------------
    0: generic trajectories
    1: golden ratio radial trajectories computed online
    2: fixed angle radial trajectories computed online
    
  */
  
  /* 
	
     -----------------------------------------------------------------------------------------
     The funcions 'preprocess_NFFT' generates data structures required for the GPU based NFFT.
     -----------------------------------------------------------------------------------------
	   
     Input: 
     ------
     - matrix_size:	     Size of each image/kspace dimension
     - matrix_size_os:	     Size of each oversampled dimension
     - W:		     Window width (in non-oversampled terms, i.e. Jackson paper)
     - domain_size_grid      Dimensions of each convolution domain (per processor)
     - domain_size_samples   Number of samples in each convolution domain (per processor)
     - domain_size_coils     Number of coils to process in each convolution (per processor)
     - trajectory:	     NDArray of dimension 2 with sample trajectory positions, i.e. size(Dim 0) = #samples, size(Dim 1) = d
	   
     Output:
     -------
     - Pointer to the computed 'NFFT_plan', 'NFFT_H_plan', or 'NFFT_iteration_plan' data structure.
	   
  */

  template< class UINTd, class FLOATd, char TYPE > NFFT_plan<UINTd, FLOATd, TYPE>*
  preprocess_NFFT( UINTd matrix_size, UINTd matrix_size_os, UINTd fixed_dims, unsigned int domain_size_samples, unsigned int domain_size_coils, float W, RealFloatArray *trajectory );

  template< class UINTd, class FLOATd, char TYPE > NFFT_H_plan<UINTd, FLOATd, TYPE>*
  preprocess_NFFT( UINTd matrix_size, UINTd matrix_size_os, UINTd fixed_dims, UINTd domain_size_grid, unsigned int domain_size_coils, float W, RealFloatArray *trajectory );

  template< class UINTd, class FLOATd, char TYPE > NFFT_iteration_plan<UINTd, FLOATd, TYPE>*
  preprocess_NFFT( UINTd matrix_size, UINTd matrix_size_os, UINTd fixed_dims, UINTd domain_size_grid, unsigned int domain_size_samples, unsigned int domain_size_coils, float W, RealFloatArray *trajectory );
}

#endif
