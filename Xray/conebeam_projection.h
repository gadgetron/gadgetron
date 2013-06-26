#pragma once

#include "hoNDArray.h"
#include "vector_td.h"

#include <vector>
#include <math_constants.h>

namespace Gadgetron {
  
  // Forwards projection of a 3D volume onto a set of projections.
  // Dependent on the provided binnning indices, 
  // only a subset of the projections may be included.
  
  void 
  conebeam_forwards_projection( hoNDArray<float> *projections,
				hoNDArray<float> *image,
				std::vector<float> angles, 
				std::vector<floatd2> offsets, 
				std::vector<unsigned int> indices,
				int projections_per_batch, 
				int num_samples_per_ray,
				floatd3 is_spacing_in_mm, 
				floatd2 ps_dims_in_mm,
				float SDD, 
				float SAD,
				bool accumulate );
  
  // Backprojection of a set of projections onto a 3D volume.
  // Dependent on the provided binnning indices, 
  // only a subset of the projections may be included.

  void 
  conebeam_backwards_projection( hoNDArray<float> *projections, 
				 hoNDArray<float> *image,
				 std::vector<float> angles, 
				 std::vector<floatd2> offsets, 
				 std::vector<unsigned int> indices,
				 int projections_per_batch,
				 uintd3 is_dims_in_pixels, 
				 floatd3 is_spacing_in_mm, 
				 floatd2 ps_dims_in_mm,
				 float SDD, 
				 float SAD,
				 bool use_fbp, 
				 bool use_oversampling_in_fbp,
				 float maximum_angle,
				 bool accumulate );
}
