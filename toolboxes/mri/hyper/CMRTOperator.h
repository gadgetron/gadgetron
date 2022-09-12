/*
 * CMRTOperator.h
 *
 *  Created on: Apr 15, 2014
 *      Author: u051747
 */
#pragma once

#include "linearOperator.h"
#include "cuNDArray.h"
#include "radial_utilities.h"
#include "cuNFFT.h"
#include "../../nfft/NFFTOperator.h"
#include "cuNDFFT.h"
#include "vector_td_operators.h"

#include "hoNDArray_fileio.h"
#include "cudaDeviceManager.h"

#include <boost/make_shared.hpp>

namespace Gadgetron {

template<class REAL> class CMRTOperator: public linearOperator<cuNDArray< complext<REAL > > > {
	typedef complext<REAL> COMPLEX;
public:
	CMRTOperator(): W_(5.5), alpha_(2),readout_oversampling_factor_(1){};
	virtual ~CMRTOperator(){};

	virtual void mult_MH(cuNDArray<COMPLEX>* in, cuNDArray<COMPLEX>* out, bool accumulate = false){
		cuNDArray<COMPLEX> projections(&projection_dims);
		E_.mult_MH(in,&projections);
		std::vector<size_t> permute_dims;
		permute_dims.push_back(0);
		permute_dims.push_back(2);
		permute_dims.push_back(1);
		permute_dims.push_back(3);
		projections = permute(projections, permute_dims);
		cuNDFFT<REAL>::instance()->fft(&projections,0u);

		COMPLEX* proj_ptr = projections.get_data_ptr();
		std::vector<size_t> proj_dim3d(projection_dims.begin(),projection_dims.end()-1);
		std::vector<size_t> out_dim3d({out->get_size(0),out->get_size(1),out->get_size(2)});
		COMPLEX* out_ptr = out->get_data_ptr();
		for (size_t t = 0; t< out->get_size(3); t++){
			cuNDArray<COMPLEX> proj_view(proj_dim3d,proj_ptr);
			cuNDArray<COMPLEX> out_view(out_dim3d,out_ptr);
			backprojections[t]->mult_MH(&proj_view,&out_view,accumulate);
			proj_ptr += proj_view.get_number_of_elements();
			out_ptr += out_view.get_number_of_elements();
		}
	}


	virtual void mult_M(cuNDArray<COMPLEX>* in, cuNDArray<COMPLEX>* out, bool accumulate = false){

		cuNDArray<COMPLEX> projections(&projection_dims_permuted);

		COMPLEX* proj_ptr = projections.get_data_ptr();
		std::vector<size_t> proj_dim3d(projection_dims.begin(),projection_dims.end()-1);
		std::vector<size_t> in_dim3d({in->get_size(0),in->get_size(1),in->get_size(2)});
		COMPLEX* in_ptr = in->get_data_ptr();
		for (size_t t = 0; t< in->get_size(3); t++){
			cuNDArray<COMPLEX> proj_view(proj_dim3d,proj_ptr);
			cuNDArray<COMPLEX> in_view(in_dim3d,in_ptr);
			backprojections[t]->mult_M(&in_view,&proj_view,accumulate);
			proj_ptr += proj_view.get_number_of_elements();
			in_ptr += in_view.get_number_of_elements();
		}


		cuNDFFT<REAL>::instance()->ifft(&projections,0u);
		std::vector<size_t> permute_dims = {0,2,1,3};
		projections = permute(projections,permute_dims);


		E_.mult_M(&projections,out,accumulate);
	}


	void setup( boost::shared_ptr<cuNDArray<vector_td<REAL,2> > > traj, std::vector<size_t>& dims,std::vector<size_t>& projection_dims, unsigned int offset, bool golden_ratio ){

		E_.setup( uint64d2(projection_dims[0], projection_dims[1]),
				uint64d2(projection_dims[0], projection_dims[1])*size_t(2), // !! <-- alpha_
				W_ );
		E_.preprocess(traj.get());

		this->projection_dims = projection_dims;
		projection_dims_permuted = projection_dims;
		projection_dims_permuted[1] = projection_dims[2];
		projection_dims_permuted[2] = projection_dims[1];

		size_t ntimeframes = projection_dims.size() > 3 ? projection_dims[3] : 1;
		/*
		boost::shared_ptr< cuNDArray<REAL> > b_dcw = compute_radial_dcw_fixed_angle_2d
						( dims[0], projection_dims[2], alpha_, 1.0f/readout_oversampling_factor_ );
		sqrt_inplace(b_dcw.get());

		//backprojection.set_dcw(b_dcw);
		 */

		size_t time_offset =offset;
		backprojections.clear();
		for (size_t t = 0; t < ntimeframes; t++){
			auto backprojection = boost::make_shared<NFFTOperator<cuNDArray,REAL,2>>();
			backprojection->setup( uint64d2(dims[0], dims[1]),
					uint64d2(dims[0], dims[1])*size_t(2), // !! <-- alpha_
					W_ );

			boost::shared_ptr< cuNDArray<floatd2> > traj2;

			if (golden_ratio){
				traj2= compute_radial_trajectory_golden_ratio_2d<REAL>
				( projection_dims[0], projection_dims[2],1, time_offset, GR_ORIGINAL );

			}
			else{
				traj2= compute_radial_trajectory_fixed_angle_2d<REAL>
				( projection_dims[0], projection_dims[2], 1/*number of frames*/ );
			}


			backprojection->preprocess(traj2.get());
			backprojections.push_back(backprojection);
			time_offset += projection_dims[2];
		}

	}


protected:

	NFFTOperator<cuNDArray,REAL,2> E_; //cuNFFTOperator reconstructing the 2d projections
	std::vector< boost::shared_ptr< NFFTOperator<cuNDArray,REAL,2>>> backprojections; //cuNFFTOperator doing the equivalent of backprojection

	std::vector<size_t> projection_dims;
	std::vector<size_t> projection_dims_permuted;
	REAL W_;
	REAL readout_oversampling_factor_;
	REAL alpha_;
};

} /* namespace Gadgetron */

