#pragma once

#include "hoCuNDArray.h"
#include "cuNDArray.h"
#include "vector_td.h"
#include "gpuxray_export.h"

namespace Gadgetron {



/**
 * Applies the sqrt of the FPB offset correction. Should be used on projection data with offset detector for iterative reconstruction.
 * @param projections
 * @param offsets
 * @param ps_dims_in_mm
 * @param SDD
 * @param SAD
 */
void apply_offset_correct(hoCuNDArray<float>* projections,std::vector<floatd2>& offsets,		floatd2 ps_dims_in_mm, float SDD,	float SAD);
/**
 *
 * @param projections
 * @param image
 * @param angles
 * @param offsets
 * @param indices
 * @param projections_per_batch
 * @param samples_per_pixel
 * @param is_dims_in_mm
 * @param ps_dims_in_mm
 * @param SDD
 * @param SAD
 * @param accumulate
 */
  // Forwards projection of a 3D volume onto a set of projections.
  // - dependening on the provided binnning indices, just a subset of the projections can be targeted.
  //
  
  EXPORTGPUXRAY void conebeam_forwards_projection
    ( hoCuNDArray<float> *projections,
				hoCuNDArray<float> *image,
				std::vector<float> angles, 
				std::vector<floatd2> offsets, 
				std::vector<unsigned int> indices,
				int projections_per_batch, 
				float samples_per_pixel,
				floatd3 is_dims_in_mm, 
				floatd2 ps_dims_in_mm,
				float SDD, 
				float SAD
  );
  
  // Backprojection of a set of projections onto a 3D volume.
  // - depending on the provided binnning indices, just a subset of the projections can be included
  //

  template <bool FBP> EXPORTGPUXRAY void conebeam_backwards_projection( 
        hoCuNDArray<float> *projections,
        hoCuNDArray<float> *image,
        std::vector<float> angles, 
        std::vector<floatd2> offsets, 
        std::vector<unsigned int> indices,
        int projections_per_batch,
        intd3 is_dims_in_pixels, 
        floatd3 is_dims_in_mm, 
        floatd2 ps_dims_in_mm,
        float SDD, 
        float SAD,
        bool short_scan,
        bool use_offset_correction,
        bool accumulate, 
        cuNDArray<float> *cosine_weights = 0x0,
        cuNDArray<float> *frequency_filter = 0x0
  );
}
