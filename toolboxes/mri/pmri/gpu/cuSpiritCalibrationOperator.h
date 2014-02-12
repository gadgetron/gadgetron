/** \file cuSpiritCalibrationOperator.h
    \brief Operator to compute Spirit calibration kernels.

    The operator domain size is the image size times the squared number of coils. 
    The codomain size is the image size times the number of coils. 
*/

#pragma once

#include "cuDiagonalSumOperator.h"

namespace Gadgetron {

  template<class REAL> class cuSpirit2DCalibrationOperator : public linearOperator< cuNDArray< complext<REAL> > >
  {
  public:
    
    cuSpirit2DCalibrationOperator() : linearOperator< cuNDArray< complext<REAL> > >() {}
    virtual ~cuSpirit2DCalibrationOperator() {}
    
    virtual void set_accumulated_kspace( boost::shared_ptr< cuNDArray< complext<REAL> > > coil_data )
    { 
      if( coil_data->get_number_of_dimensions() != 3 ){
        throw std::runtime_error("cuSpirit2DCalibrationOperator::set_accumulated_kspace: coil data must be three-dimensionsal");
      }
    
      D_ = boost::shared_ptr< cuDiagonalSumOperator< complext<REAL> > >(new cuDiagonalSumOperator< complext<REAL> >());
      D_->set_diagonal(coil_data);
    }
  
    virtual void mult_M( cuNDArray< complext<REAL> > *in, cuNDArray< complext<REAL> > *out, bool accumulate = false )
    {
      if( !D_.get() ){
        throw std::runtime_error("cuSpiritCalibrationOperator::mult_M failed: accumulated kspace data not set");
      }
      
      if( in->get_number_of_dimensions() != 3 || out->get_number_of_dimensions() != 3 ){
        throw std::runtime_error("cuSpiritCalibrationOperator::mult_M failed: expected exactly 3 dimensions in input and output images");
      }

      boost::shared_ptr< cuNDArray< complext<REAL> > > diagonal = D_->get_diagonal();

      const unsigned int num_coils = diagonal->get_size(2);
      const unsigned int num_phases_in = in->get_size(2);
      const unsigned int num_phases_out = out->get_size(2);
      
      if( num_phases_in != num_coils*num_coils ){
        throw std::runtime_error("cuSpirit2DCalibrationOperator::mult_M failed: the input image array's last dimension must have size #coils squared");
      }

      if( num_phases_out != num_coils ){
        throw std::runtime_error("cuSpirit2DCalibrationOperator::mult_M failed: the output image size must match that of the provided coil array");
      }
      
      std::vector<size_t> dim_coils = *diagonal->get_dimensions();
      std::vector<size_t> dim_image = dim_coils; dim_image.pop_back();

      size_t num_elements_image = dim_coils[0]*dim_coils[1];
      size_t num_elements_coils = num_elements_image*dim_coils[2];

      // Iterate over the coils
      //
      
      for( unsigned int i=0; i<num_coils; i++ ){
        
        cuNDArray< complext<REAL> > tmp_in( &dim_coils, in->get_data_ptr()+i*num_elements_coils );
        cuNDArray< complext<REAL> > tmp_out( &dim_image, out->get_data_ptr()+i*num_elements_image );
        
        D_->mult_M( &tmp_in, &tmp_out, accumulate );
      }
    }
    
    virtual void mult_MH( cuNDArray< complext<REAL> > *in, cuNDArray< complext<REAL> > *out, bool accumulate = false )
    {
      if( !D_.get() ){
        throw std::runtime_error("cuSpiritCalibrationOperator::mult_MH failed: accumulated kspace data not set");
      }
      
      if( in->get_number_of_dimensions() != 3 || out->get_number_of_dimensions() != 3 ){
        throw std::runtime_error("cuSpiritCalibrationOperator::mult_MH failed: expected exactly 3 dimensions in input and output images");
      }

      boost::shared_ptr< cuNDArray< complext<REAL> > > diagonal = D_->get_diagonal();

      const unsigned int num_coils = diagonal->get_size(2);
      const unsigned int num_phases_in = in->get_size(2);
      const unsigned int num_phases_out = out->get_size(2);
      
      if( num_phases_in != num_coils ){
        throw std::runtime_error("cuSpirit2DCalibrationOperator::mult_MH failed: the input image size must match that of the provided coil array");
      }
      
      if( num_phases_out != num_coils*num_coils ){
        throw std::runtime_error("cuSpirit2DCalibrationOperator::mult_MH failed: the output image array's last dimension must have size #coils squared");
      }      
      
      std::vector<size_t> dim_coils = *diagonal->get_dimensions();
      std::vector<size_t> dim_image = dim_coils; dim_image.pop_back();

      size_t num_elements_image = dim_coils[0]*dim_coils[1];
      size_t num_elements_coils = num_elements_image*dim_coils[2];

      // Iterate over the coils
      //
      
      for( unsigned int i=0; i<num_coils; i++ ){
        
        cuNDArray< complext<REAL> > tmp_in( &dim_image, in->get_data_ptr()+i*num_elements_image );
        cuNDArray< complext<REAL> > tmp_out( &dim_coils, out->get_data_ptr()+i*num_elements_coils );
        
        D_->mult_MH( &tmp_in, &tmp_out, accumulate );
      }
    }
    
    virtual boost::shared_ptr< linearOperator< cuNDArray< complext<REAL> > > > clone() {
      return linearOperator< cuNDArray< complext<REAL> > >::clone(this);
    }

  protected:    
    boost::shared_ptr< cuDiagonalSumOperator< complext<REAL> > > D_;
  };
}
