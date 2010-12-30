
//
// Generate the necessary data structures for the CUDA implementation of the non-Cartesian Sense reconstruction
//


#ifndef _NSENSE_HPP_
#define _NSENSE_HPP_

#include "preprocess.hpp"

namespace mr_recon
{

	/* 

	---------------------------------------------------------------------------------------------------------------------------
	The funcion 'preprocess_NSense' generates the data structure required for the GPU based non-Cartesian Sense reconstruction.
	---------------------------------------------------------------------------------------------------------------------------

	Input: 
	------
  	- matrix_size:		Size of each image/kspace dimension
	- matrix_size_os:	Size of each oversampled dimension
	- W:				Window width (in non-oversampled terms, i.e. Jackson paper)
	- trajectory:		NDArray of dimension 2 with sample trajectory positions, i.e. size(Dim 0) = #samples, size(Dim 1) = d
	
	Output:
	-------
	- Pointer to the computed 'preprocess_NSense' object.

	*/

  template< class UINTd, class FLOATd, char TYPE > NFFT_iteration_plan<UINTd, FLOATd, TYPE>*
		preprocess_NSense( UINTd matrix_size, UINTd matrix_size_os, UINTd fixed_dims, UINTd domain_size_grid, unsigned int domain_size_samples, unsigned int domain_size_coils, float W, RealFloatArray *trajectory );

  template< class UINTd, class FLOATd, char TYPE > NFFT_iteration_plan<UINTd, FLOATd, TYPE>*
		preprocess_NSense_radial( UINTd matrix_size, UINTd matrix_size_os, UINTd fixed_dims, UINTd domain_size_grid, unsigned int domain_size_samples, unsigned int domain_size_coils, float W, unsigned int num_projections, unsigned int samples_per_projection, unsigned int projections_per_frame, unsigned int angular_offset = 0, unsigned int frames_per_rotation = 1, float gc_factor  = 1.0f );
}

#endif
