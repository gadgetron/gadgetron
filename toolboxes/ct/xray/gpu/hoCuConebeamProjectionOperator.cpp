#include "hoCuConebeamProjectionOperator.h"
#include "conebeam_projection.h"
#include "vector_td_operators.h"
#include "cuNDArray_operators.h"
#include "cuNDArray_blas.h"
#include "GPUTimer.h"

#include <vector>
#include <stdio.h>

namespace Gadgetron
{
  void hoCuConebeamProjectionOperator
  ::compute_default_frequency_filter()
  {
    // This code computes the default frequency filter used in filtered backprojection
    // _Important_ aspects:
    // - the filter is defined as single precision weights (non-complex)
    // - the filter defines the scalar weights for the positive frequencies only (i.e "one side)
    //   - however, the size of the filter still equals the full size of the 1D dimension to filter +1 ...
    //   - ... due to zero padding and cufft expecting an additional element.
    //

    if( !preprocessed_ )
      throw std::runtime_error("Error: hoCuConebeamProjectionOperator::compute_default_frequency_filter() : setup not performed");
   
    std::vector<size_t> dims;
    dims.push_back(acquisition_->get_projections()->get_size(0)+1);
    
    hoCuNDArray<float> host_weights(&dims);
    float* data = host_weights.get_data_ptr();

    const float A2 = dims[0]*dims[0];

#ifdef USE_OMP
#pragma omp parallel for
#endif    
    for( int i=0; i<dims[0]; i++ ) {
      float k = float(i);
      data[i] = k*A2/(A2-k*k)*std::exp(-A2/(A2-k*k)); // From Guo et al, Journal of X-Ray Science and Technology 2011, doi: 10.3233/XST-2011-0294
    }

    frequency_filter_ = boost::shared_ptr< cuNDArray<float> >(new cuNDArray<float>(&host_weights));
    float sum = asum(frequency_filter_.get());
    *frequency_filter_ /= (dims[0]/sum);
  }
  
  void hoCuConebeamProjectionOperator
  ::compute_cosine_weights()
  {
    if( !preprocessed_ )
      throw std::runtime_error("Error: hoCuConebeamProjectionOperator::compute_cosine_weights() : setup not performed");

    uintd2 ps_dims_in_pixels( acquisition_->get_projections()->get_size(0), acquisition_->get_projections()->get_size(1) );
    floatd2 ps_dims_in_mm = acquisition_->get_geometry()->get_FOV();

    double SAD = double(acquisition_->get_geometry()->get_SAD());
    double SDD = double(acquisition_->get_geometry()->get_SDD());

    std::vector<size_t> dims;
    dims.push_back(ps_dims_in_pixels[0]);
    dims.push_back(ps_dims_in_pixels[1]);

    hoCuNDArray<float> weights(&dims);
    float* data = weights.get_data_ptr();

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for(  int y=0; y<ps_dims_in_pixels[1]; y++ ) {
      for( int x=0; x<ps_dims_in_pixels[0]; x++ ) {

        double xx = (( double(x) / double(ps_dims_in_pixels[0])) - 0.5) * ps_dims_in_mm[0];
        double yy = (( double(y) / double(ps_dims_in_pixels[1])) - 0.5) * ps_dims_in_mm[1];
        double s = SAD * xx/SDD;
        double v = SAD * yy/SDD;

        // Equation 10.1, page 386 in Computed Tomography 2nd edition, Jiang Hsieh
        //

        double value = SAD / std::sqrt( SAD*SAD + s*s + v*v );
        data[x+y*ps_dims_in_pixels[0]] = float(value);
      }
    }
    cosine_weights_ = boost::shared_ptr< cuNDArray<float> >(new cuNDArray<float>(&weights));
  }

  void hoCuConebeamProjectionOperator
  ::mult_M( hoCuNDArray<float> *image, hoCuNDArray<float> *projections, bool accumulate )
  {

    // Validate the input 
    //

    if( image == 0x0 || projections == 0x0 ){
      throw std::runtime_error("Error: hoCuConebeamProjectionOperator::mult_M: illegal array pointer provided");
    }
    
    if( (image->get_number_of_dimensions() != 4) &&  (image->get_number_of_dimensions() != 3) ){
      throw std::runtime_error("Error: hoCuConebeamProjectionOperator::mult_M: image array must be four or three -dimensional");
    }

    if( projections->get_number_of_dimensions() != 3 ){
      throw std::runtime_error("Error: hoCuConebeamProjectionOperator::mult_M: projections array must be three-dimensional");
    }
    
    if( !preprocessed_ ){
      throw std::runtime_error( "Error: hoCuConebeamProjectionOperator::mult_M: setup not performed");
    }

    if( !binning_.get() ){
      throw std::runtime_error( "Error: hoCuConebeamProjectionOperator::mult_M: binning not provided");
    }

    if( projections->get_size(2) != acquisition_->get_geometry()->get_angles().size() || 
        projections->get_size(2) != acquisition_->get_geometry()->get_offsets().size() ){
      throw std::runtime_error("Error: hoCuConebeamProjectionOperator::mult_M: inconsistent sizes of input arrays/vectors");
    }

    // Iterate over the temporal dimension.
    // I.e. reconstruct one 3D volume at a time.
    //

    for( int b=0; b<binning_->get_number_of_bins(); b++ ) {

      floatd2 ps_dims_in_pixels_float(projections->get_size(0), projections->get_size(1));
      floatd2 ps_dims_in_mm = acquisition_->get_geometry()->get_FOV();
      floatd2 ps_spacing_in_mm = ps_dims_in_mm / ps_dims_in_pixels_float;

      float SDD = acquisition_->get_geometry()->get_SDD();
      float SAD = acquisition_->get_geometry()->get_SAD();

      std::vector<size_t> dims_3d = *image->get_dimensions();
      if (dims_3d.size()==4)
      	dims_3d.pop_back();
      
      int num_3d_elements = dims_3d[0]*dims_3d[1]*dims_3d[2];

      hoCuNDArray<float> image_3d(&dims_3d, image->get_data_ptr()+b*num_3d_elements);

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

  void hoCuConebeamProjectionOperator
  ::mult_MH( hoCuNDArray<float> *projections, hoCuNDArray<float> *image, bool accumulate )
  {

    // Validate the input 
    //

    if( image == 0x0 || projections == 0x0 ){
      throw std::runtime_error("Error: hoCuConebeamProjectionOperator::mult_MH:: illegal array pointer provided");
    }
    
    if( (image->get_number_of_dimensions() != 4) &&  (image->get_number_of_dimensions() != 3) ){
      throw std::runtime_error("Error: hoCuConebeamProjectionOperator::mult_MH: image array must be four or three -dimensional");
    }
    
    if( projections->get_number_of_dimensions() != 3 ){
      throw std::runtime_error("Error: hoCuConebeamProjectionOperator::mult_MH: projections array must be three-dimensional");
    }
    
    if( !preprocessed_ ){
      throw std::runtime_error( "Error: hoCuConebeamProjectionOperator::mult_MH: setup not performed");
    }

    if( !binning_.get() ){
      throw std::runtime_error( "Error: hoCuConebeamProjectionOperator::mult_MH: binning not provided");
    }

    if( projections->get_size(2) != acquisition_->get_geometry()->get_angles().size() ||
        projections->get_size(2) != acquisition_->get_geometry()->get_offsets().size() ){
      throw std::runtime_error("Error: hoCuConebeamProjectionOperator::mult_MH: inconsistent sizes of input arrays/vectors");
    }

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

      std::vector<size_t> dims_3d = *image->get_dimensions();
      if (dims_3d.size() ==4)
      	dims_3d.pop_back();

      int num_3d_elements = dims_3d[0]*dims_3d[1]*dims_3d[2];

      hoCuNDArray<float> image_3d(&dims_3d, image->get_data_ptr()+b*num_3d_elements);

      if( use_fbp_ ){

        if( !cosine_weights_.get() )
          compute_cosine_weights();

        if( !frequency_filter_.get() )
          compute_default_frequency_filter();

        conebeam_backwards_projection<true>
          ( projections, &image_3d,
            acquisition_->get_geometry()->get_angles(), 
            acquisition_->get_geometry()->get_offsets(),
            binning_->get_bin(b),
            projections_per_batch_,
            is_dims_in_pixels, is_dims_in_mm_, ps_dims_in_mm,
            SDD, SAD, short_scan_, use_offset_correction_, accumulate,
            cosine_weights_.get(), frequency_filter_.get() );
      }
      else
        conebeam_backwards_projection<false>
          ( projections, &image_3d,
            acquisition_->get_geometry()->get_angles(), 
            acquisition_->get_geometry()->get_offsets(),
            binning_->get_bin(b),
            projections_per_batch_,
            is_dims_in_pixels, is_dims_in_mm_, ps_dims_in_mm,
            SDD, SAD, short_scan_, use_offset_correction_, accumulate );
    }
  }
}
