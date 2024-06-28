/** \file cuSpiritOperator.h
    \brief Spirit regularization operator.

    The operator domain and codomain sizes are the image size times the number of coils. 
 */

#pragma once

#include "cuDiagonalSumOperator.h"
#include <numeric>

namespace Gadgetron {

template<class REAL> class cuSpirit2DOperator : public linearOperator< cuNDArray< complext<REAL> > >
{
public:

	cuSpirit2DOperator() : linearOperator< cuNDArray< complext<REAL> > >() {
		D_ = boost::shared_ptr< cuDiagonalSumOperator< complext<REAL> > >(new cuDiagonalSumOperator< complext<REAL> >());
	}

	virtual ~cuSpirit2DOperator() {}

	virtual void set_calibration_kernels( boost::shared_ptr< cuNDArray< complext<REAL> > > kernels )
	{
		if( kernels->get_number_of_dimensions() !=3 ){
			throw std::runtime_error("cuSpirit2DOperator::set_calibration kernels: kernels array must be three-dimensionsal (x,y,squared num coils)");
		}
		kernels_ = kernels;
	}

	virtual void mult_M( cuNDArray< complext<REAL> > *in, cuNDArray< complext<REAL> > *out, bool accumulate = false )
	{
		if( !kernels_.get() ){
			throw std::runtime_error("cuSpiritCalibrationOperator::mult_M failed: calibration kernels not set");
		}

		if( in->get_number_of_dimensions() < 3 || out->get_number_of_dimensions() < 3 ){
			throw std::runtime_error("cuSpiritCalibrationOperator::mult_M failed: expected at least 3 dimensions in input and output images");
		}

		const unsigned int num_coils_squared = kernels_->get_size(2);
		const unsigned int num_phases_in = in->get_size(in->get_number_of_dimensions()-1);
		const unsigned int num_phases_out = out->get_size(out->get_number_of_dimensions()-1);
		const unsigned int num_frames = in->get_number_of_dimensions() == 3 ? 1 : in->get_size(2); //Number of frames 3rd dimension

		if( num_phases_out != num_phases_in ){
			throw std::runtime_error("cuSpirit2DOperator::mult_M failed: array size mismatch between input/output images");
		}

		if( num_phases_in*num_phases_out != num_coils_squared ){
			throw std::runtime_error("cuSpirit2DOperator::mult_M failed: the calibration kernels do not correspond to the squared number of coils");
		}

		std::vector<size_t> dim_image = {in->get_size(0),in->get_size(1)};
		std::vector<size_t> dim_coils = dim_image;
		dim_coils.push_back(num_phases_in);


		size_t num_elements_image = std::accumulate(dim_image.begin(),dim_image.end(),1,std::multiplies<size_t>());
		size_t num_elements_coils = num_elements_image*dim_coils.back();

		// Iterate over the coils
		//
		for( unsigned int i=0; i<num_phases_out; i++ ){
			boost::shared_ptr< cuNDArray< complext<REAL> > > tmp_kernel( new cuNDArray< complext<REAL> >(dim_coils, kernels_->get_data_ptr()+i*num_elements_coils ));
			D_->set_diagonal(tmp_kernel);
			for (unsigned int k=0; k < num_frames; k++){
				cuNDArray<complext<REAL>> tmp_in(dim_image,in->get_data_ptr()+k*num_elements_image);
				cuNDArray< complext<REAL> > tmp_out(dim_image, out->get_data_ptr()+i*num_elements_image*num_frames+k*num_elements_image );
				D_->mult_M( &tmp_in, &tmp_out, accumulate );
			}
		}

		// Subtract identity
		//

		*out -= *in;
	}

	virtual void mult_MH( cuNDArray< complext<REAL> > *in, cuNDArray< complext<REAL> > *out, bool accumulate = false )
	{
		if( !kernels_.get() ){
			throw std::runtime_error("cuSpiritCalibrationOperator::mult_MH failed: calibration kernels not set");
		}

		if( in->get_number_of_dimensions() != 3 || out->get_number_of_dimensions() != 3 ){
			throw std::runtime_error("cuSpiritCalibrationOperator::mult_MH failed: expected exactly 3 dimensions in input and output images");
		}
const unsigned int num_coils_squared = kernels_->get_size(2);
		const unsigned int num_phases_in = in->get_size(in->get_number_of_dimensions()-1);
		const unsigned int num_phases_out = out->get_size(out->get_number_of_dimensions()-1);
		const unsigned int num_frames = in->get_number_of_dimensions() == 3 ? 1 : in->get_size(2); //Number of frames 3rd dimension

		if( num_phases_out != num_phases_in ){
			throw std::runtime_error("cuSpirit2DOperator::mult_M failed: array size mismatch between input/output images");
		}

		if( num_phases_in*num_phases_out != num_coils_squared ){
			throw std::runtime_error("cuSpirit2DOperator::mult_M failed: the calibration kernels do not correspond to the squared number of coils");
		}

		std::vector<size_t> dim_image = {in->get_size(0),in->get_size(1)};
		std::vector<size_t> dim_coils = dim_image;
		dim_coils.push_back(num_phases_in);


		size_t num_elements_image = std::accumulate(dim_image.begin(),dim_image.end(),1,std::multiplies<size_t>());
		size_t num_elements_coils = num_elements_image*dim_coils.back();


		// Iterate over the coils
		//

		for( unsigned int i=0; i<num_phases_in; i++ ){
			boost::shared_ptr< cuNDArray< complext<REAL> > > tmp_kernel( new cuNDArray< complext<REAL> >(dim_coils, kernels_->get_data_ptr()+i*num_elements_coils ));

			D_->set_diagonal(tmp_kernel);
			for (unsigned int k=0; i<num_frames; k++){
			cuNDArray< complext<REAL> > tmp_in(dim_image, in->get_data_ptr()+i*num_elements_image*num_frames+k*num_elements_image );
			cuNDArray<complext<REAL> > tmp_out(dim_image,out->get_data_ptr()+k*num_elements_image);
			if( i==0 && !accumulate )
				D_->mult_MH( &tmp_in, &tmp_out, false );
			else
				D_->mult_MH( &tmp_in, &tmp_out, true );
		}
		}

		// Subtract identity
		//

		*out -= *in;
	}



protected:
	boost::shared_ptr< cuNDArray< complext<REAL> > > kernels_;
	boost::shared_ptr< cuDiagonalSumOperator< complext<REAL> > > D_;
};
}
