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
#include "cuNFFTOperator.h"
#include "cuNDFFT.h"
#include "vector_td_operators.h"

#include "hoNDArray_fileio.h"
#include "cudaDeviceManager.h"
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
		projections = *permute(&projections,&permute_dims);
		cuNDFFT<REAL>::instance()->fft(&projections,0u);
		backprojection.mult_MH(&projections,out,accumulate);
	}


	virtual void mult_M(cuNDArray<COMPLEX>* in, cuNDArray<COMPLEX>* out, bool accumulate = false){

		cuNDArray<COMPLEX> projections(&projection_dims_permuted);
		backprojection.mult_M(in,&projections);

		cuNDFFT<REAL>::instance()->ifft(&projections,0u);
		std::vector<size_t> permute_dims;
		permute_dims.push_back(0);
		permute_dims.push_back(2);
		permute_dims.push_back(1);
		projections = *permute(&projections,&permute_dims);


		E_.mult_M(&projections,out,accumulate);
	}


	void setup(boost::shared_ptr<cuNDArray<COMPLEX> > data, boost::shared_ptr<cuNDArray<REAL> > dcw, boost::shared_ptr<cuNDArray<vector_td<REAL,2> > > traj, std::vector<size_t>& dims,std::vector<size_t>& projection_dims, bool golden_ratio ){
		E_.set_dcw(dcw);

		E_.setup( uint64d2(projection_dims[0], projection_dims[1]),
		               uint64d2(projection_dims[0], projection_dims[1])*size_t(2), // !! <-- alpha_
		               W_ );
		E_.preprocess(traj.get());

		this->projection_dims = projection_dims;
		projection_dims_permuted = projection_dims;
		projection_dims_permuted[1] = projection_dims[2];
		projection_dims_permuted[2] = projection_dims[1];
		projection_dims_permuted[0] = projection_dims[0];

/*
		boost::shared_ptr< cuNDArray<REAL> > b_dcw = compute_radial_dcw_fixed_angle_2d
						( dims[0], projection_dims[2], alpha_, 1.0f/readout_oversampling_factor_ );
		sqrt_inplace(b_dcw.get());

		//backprojection.set_dcw(b_dcw);
*/

		*data *= *dcw;

		backprojection.setup( uint64d2(dims[0], dims[1]),
				               uint64d2(dims[0], dims[1])*size_t(2), // !! <-- alpha_
				               W_ );

		boost::shared_ptr< cuNDArray<floatd2> > traj2;

		if (golden_ratio)
			traj2= compute_radial_trajectory_golden_ratio_2d<REAL>
							( dims[0], projection_dims[2], 1, 0, GR_ORIGINAL );
		else
			traj2= compute_radial_trajectory_fixed_angle_2d<REAL>
							( dims[0], projection_dims[2], 1 /*number of frames*/ );
		backprojection.preprocess(traj2.get());

	}

  virtual boost::shared_ptr< linearOperator< cuNDArray<complext<REAL> > > > clone(){
      return linearOperator<cuNDArray<complext<REAL> > >::clone(this);
    }
protected:

	cuNFFTOperator<REAL,2> E_; //cuNFFTOperator reconstructing the 2d projections
	cuNFFTOperator<REAL,2> backprojection; //cuNFFTOperator doing the equivalent of backprojection

	std::vector<size_t> projection_dims;
	std::vector<size_t> projection_dims_permuted;
	REAL W_;
	REAL readout_oversampling_factor_;
	REAL alpha_;
};

} /* namespace Gadgetron */

