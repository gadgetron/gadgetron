#include "hoCudaConebeamProjectionOperator.h"
#include "conebeam_projection.h"
#include "hoNDArray.h"
#include "vector_td_operators.h"
#include "GPUTimer.h"

#include <vector>
#include <stdio.h>

namespace Gadgetron
{
  void hoCudaConebeamProjectionOperator
  ::mult_M( hoNDArray<float> *image, hoNDArray<float> *projections, bool accumulate )
  {
    //
    // Validate the input 
    //
    
    if( image == 0x0 || projections == 0x0 ){
      throw std::runtime_error("Error: hoCudaCobebeamProjectionOperator::mult_M: illegal array pointer provided");
    }
    
    if( image->get_number_of_dimensions() != 4 ){
      throw std::runtime_error("Error: hoCudaCobebeamProjectionOperator::mult_M: image array must be four-dimensional");
    }

    if( projections->get_number_of_dimensions() != 3 ){
      throw std::runtime_error("Error: hoCudaCobebeamProjectionOperator::mult_M: projections array must be three-dimensional");
    }
    
    if( !preprocessed_ ){
      throw std::runtime_error( "Error: hoCudaConebeamProjectionOperator::mult_M: setup not performed");
    }

    if( projections->get_size(2) != acquisition_->get_geometry()->get_angles().size() || 
	projections->get_size(2) != acquisition_->get_geometry()->get_offsets().size() ){
      throw std::runtime_error("Error: hoCudaCobebeamProjectionOperator::mult_M: inconsistent sizes of input arrays/vectors");
    }

    GPUTimer timer("Conebeam forwards projection");

    // Iterate over the temporal dimension.
    // I.e. reconstruct one 3D volume at a time.
    //

    for( int b=0; b<binning_->get_number_of_bins(); b++ ) {

      floatd2 ps_dims_in_pixels_float(projections->get_size(0), projections->get_size(1));
      floatd2 ps_dims_in_mm = acquisition_->get_geometry()->get_FOV();
      floatd2 ps_spacing_in_mm = ps_dims_in_mm / ps_dims_in_pixels_float;

      float SDD = acquisition_->get_geometry()->get_SDD();
      float SAD = acquisition_->get_geometry()->get_SAD();

      std::vector<unsigned int> dims_3d = *image->get_dimensions();
      dims_3d.pop_back();
      
      int num_3d_elements = dims_3d[0]*dims_3d[1]*dims_3d[2];

      hoNDArray<float> image_3d(&dims_3d, image->get_data_ptr()+b*num_3d_elements);

      conebeam_forwards_projection( projections, &image_3d, 
				    acquisition_->get_geometry()->get_angles(), 
				    acquisition_->get_geometry()->get_offsets(),
				    binning_->get_bin(b),
				    projections_per_batch_, samples_per_pixel_,
				    is_dims_in_mm_, ps_dims_in_mm, 
				    SDD, SAD, 
				    accumulate );
    }
  }

  void hoCudaConebeamProjectionOperator
  ::mult_MH( hoNDArray<float> *projections, hoNDArray<float> *image, bool accumulate )
  {
    //
    // Validate the input 
    //
    
    if( image == 0x0 || projections == 0x0 ){
      throw std::runtime_error("Error: hoCudaCobebeamProjectionOperator::mult_MH:: illegal array pointer provided");
    }
    
    if( image->get_number_of_dimensions() != 4 ){
      throw std::runtime_error("Error: hoCudaCobebeamProjectionOperator::mult_MH: image array must be four-dimensional");
    }

    if( projections->get_number_of_dimensions() != 3 ){
      throw std::runtime_error("Error: hoCudaCobebeamProjectionOperator::mult_MH: projections array must be three-dimensional");
    }
    
    if( !preprocessed_ ){
      throw std::runtime_error( "Error: hoCudaConebeamProjectionOperator::mult_MH: setup not performed");
    }

    if( projections->get_size(2) != acquisition_->get_geometry()->get_angles().size() ||
	projections->get_size(2) != acquisition_->get_geometry()->get_offsets().size() ){
      throw std::runtime_error("Error: hoCudaCobebeamProjectionOperator::mult_MH: inconsistent sizes of input arrays/vectors");
    }

    GPUTimer timer("Conebeam backwards projection");

    // Iterate over the temporal dimension.
    // I.e. reconstruct one 3D volume at a time.
    //

    for( int b=0; b<binning_->get_number_of_bins(); b++ ) {

      floatd2 ps_dims_in_pixels_float(projections->get_size(0), projections->get_size(1));
      floatd2 ps_dims_in_mm = acquisition_->get_geometry()->get_FOV();
      floatd2 ps_spacing_in_mm = ps_dims_in_mm / ps_dims_in_pixels_float;

      intd3 is_dims_in_pixels( image->get_size(0), image->get_size(1), image->get_size(2) );

      float SDD = acquisition_->get_geometry()->get_SDD();
      float SAD = acquisition_->get_geometry()->get_SAD();

      std::vector<unsigned int> dims_3d = *image->get_dimensions();
      dims_3d.pop_back();

      int num_3d_elements = dims_3d[0]*dims_3d[1]*dims_3d[2];

      hoNDArray<float> image_3d(&dims_3d, image->get_data_ptr()+b*num_3d_elements);

      conebeam_backwards_projection( projections, &image_3d,
				     acquisition_->get_geometry()->get_angles(), 
				     acquisition_->get_geometry()->get_offsets(),
				     binning_->get_bin(b),
				     projections_per_batch_,
				     is_dims_in_pixels, is_dims_in_mm_,
				     ps_dims_in_mm,
				     SDD, SAD, 
				     use_fbp_, use_oversampling_in_fbp_,
      				     max_angle_, accumulate );
    }
  }
}
