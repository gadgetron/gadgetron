/*
 * gpuCSICoilEstimationGadget.cpp
 *
 *  Created on: Nov 12, 2014
 *      Author: dch
 */

#include "gpuCSICoilEstimationGadget.h"
#include <utility>
#include <ismrmrd/xml.h>
#include "cuNDArray.h"
#include "cuNFFT.h"
#include "b1_map.h"
#include "vector_td_utilities.h"
#include "../../toolboxes/nfft/NFFTOperator.h"
#include "cuCgSolver.h"
#include "cuNDArray_fileio.h"
#include "trajectory_utils.h"
#include <boost/make_shared.hpp>
#include "hoNDArray_math.h"

namespace Gadgetron {

gpuCSICoilEstimationGadget::gpuCSICoilEstimationGadget() {

}

gpuCSICoilEstimationGadget::~gpuCSICoilEstimationGadget() {
	// TODO Auto-generated destructor stub
}

int gpuCSICoilEstimationGadget::process_config(ACE_Message_Block* mb) {
	ISMRMRD::IsmrmrdHeader h;
	ISMRMRD::deserialize(mb->rd_ptr(),h);

	if (h.encoding.size() != 1){
		GDEBUG("Coil estimation gadget only supports encoding spaces of size 1\n");
		return GADGET_FAIL;
	}
	img_size = {h.encoding[0].reconSpace.matrixSize.x,h.encoding[0].reconSpace.matrixSize.y};
	kernel_width_ = kernel_width.value();

	coils = h.acquisitionSystemInformation->receiverChannels;
	skip_lines_ = skip_lines.value();
	return GADGET_OK;
}

int gpuCSICoilEstimationGadget::process(
		GadgetContainerMessage<IsmrmrdReconData>* m1) {
	IsmrmrdReconData* bucket = m1->getObjectPtr();

	auto cm1 = new GadgetContainerMessage<cuSenseData>();
	auto senseData = cm1->getObjectPtr();

	coils = bucket->rbit_.front().data_.headers_[0].active_channels;
	GDEBUG("Active channels %i \n",coils);


	{

		hoNDArray<std::complex<float>> * ho_data = &bucket->rbit_.front().data_.data_;
		hoNDArray<float>* ho_traj = &bucket->rbit_.front().data_.trajectory_.value();


		if (skip_lines_ > 0){



			auto seperated = split_calibration_lines<std::complex<float>>(*ho_data,skip_lines_,1);

			senseData->freq_calibration = boost::make_shared<cuNDArray<float_complext>>((hoNDArray<float_complext>*)std::get<0>(seperated).get());
			senseData->freq_calibration->squeeze();
			senseData->data = boost::make_shared<cuNDArray<float_complext>>((hoNDArray<float_complext>*)std::get<1>(seperated).get());
		} else {

			senseData->data = boost::make_shared<cuNDArray<float_complext>>(reinterpret_cast<hoNDArray<float_complext>*>(ho_data));
		}


		if (ho_traj->get_size(0) > 2){ //We have dcw
			auto traj_dcw = separate_traj_and_dcw(ho_traj);
			senseData->traj = boost::make_shared<cuNDArray<floatd2>>(*std::get<0>(traj_dcw));
			senseData->dcw = boost::make_shared<cuNDArray<float>>(*std::get<1>(traj_dcw));
		} else {
			std::vector<size_t> tdims = *ho_traj->get_dimensions();
			std::vector<size_t> tmp_dim(tdims.begin()+1,tdims.end());
			hoNDArray<floatd2> tmp(tmp_dim,reinterpret_cast<floatd2*>(ho_traj->get_data_ptr()));
			senseData->traj = boost::make_shared<cuNDArray<floatd2>>(tmp);
		}

	}


	//Remove Initial Spirals

	boost::shared_ptr< cuNDArray<float_complext> >  ref_data;
	boost::shared_ptr< cuNDArray<floatd2> > ref_traj;
	boost::shared_ptr<cuNDArray<float> > ref_dcw;



	if (!bucket->rbit_.front().ref_){
		GDEBUG("Setting reference data to real data\n");
		ref_data = senseData->data;
		ref_traj = senseData->traj;
		ref_dcw = senseData->dcw;
	} else {

		auto & ref = *bucket->rbit_.front().ref_;

		hoNDArray<std::complex<float>> * ho_data = &ref.data_;
		hoNDArray<float>* ho_traj = &ref.trajectory_.value();

		ref_data = boost::make_shared<cuNDArray<float_complext>>(reinterpret_cast<hoNDArray<float_complext>*>(ho_data));
		if (ho_traj->get_size(0) > 2){
			auto traj_dcw = separate_traj_and_dcw(ho_traj);
			ref_traj =boost::make_shared<cuNDArray<floatd2>>(*std::get<0>(traj_dcw));
			ref_dcw = boost::make_shared<cuNDArray<float>>(*std::get<1>(traj_dcw));
		} else {
			std::vector<size_t> tdims = *ho_traj->get_dimensions();
			std::vector<size_t> tmp_dim(tdims.begin()+1,tdims.end());
			hoNDArray<floatd2> tmp(tmp_dim,reinterpret_cast<floatd2*>(ho_traj->get_data_ptr()));
			ref_traj = boost::make_shared<cuNDArray<floatd2>>(tmp);
		}
	}


	senseData->csm = calculate_CSM(ref_data.get(),ref_traj.get(),ref_dcw.get());



	if (this->next()->putq(cm1) == GADGET_FAIL){
		GERROR("Failed to put message on que\n");
		return GADGET_FAIL;
	}

	return GADGET_OK;


}


boost::shared_ptr<cuNDArray<float_complext> > gpuCSICoilEstimationGadget::calculate_CSM(
		cuNDArray<float_complext>* data, cuNDArray<floatd2>* traj, cuNDArray<float>* dcw ) {

	std::vector<size_t> csm_dims = img_size;
	csm_dims.push_back(coils);
	std::vector<size_t> flat_dims = {traj->get_number_of_elements()};
	cuNDArray<floatd2> flat_traj(flat_dims,traj->get_data_ptr());
	std::vector<size_t> spiral_dims{data->get_size(0)*data->get_size(1)*data->get_size(2),data->get_size(3)}; //Trajectories, coils
	cuNDArray<complext<float>> second_spiral(spiral_dims,data->get_data_ptr());
	std::vector<size_t> spiral_traj_dims{spiral_dims[0]};
	cuNDArray<floatd2> spiral_traj(spiral_traj_dims,traj->get_data_ptr());
	if (dcw) { //We have density compensation, so we can get away with gridding

		cuNFFT_impl<float,2> plan(from_std_vector<size_t,2>(img_size),from_std_vector<size_t,2>(img_size)*size_t(2),kernel_width_);
		cuNDArray<float_complext> tmp(csm_dims);
		GDEBUG("Coils %i \n\n",tmp.get_size(2));

		cuNDArray<float> spiral_dcw(spiral_traj_dims,dcw->get_data_ptr());

		GDEBUG("Preprocessing\n\n");
		plan.preprocess(spiral_traj,NFFT_prep_mode::NC2C);
		GDEBUG("Computing\n\n");

		auto tmp_abs = abs(&tmp);

		return boost::make_shared<cuNDArray<float_complext>>(estimate_b1_map<float,2>(tmp));

	} else { //No density compensation, we have to do iterative reconstruction.


		auto E = boost::make_shared<NFFTOperator<cuNDArray,float,2>>();

		E->setup(from_std_vector<size_t,2>(img_size),from_std_vector<size_t,2>(img_size)*size_t(2),kernel_width_);

		E->set_codomain_dimensions(&spiral_dims);
		E->set_domain_dimensions(&csm_dims);
		cuCgSolver<float_complext> solver;
		solver.set_max_iterations(20);
		solver.set_encoding_operator(E);

		E->preprocess(&spiral_traj);

		auto tmp = solver.solve(&second_spiral);
		auto tmp_abs = abs(tmp.get());



		auto res = estimate_b1_map<float,2>(*tmp);
		//fill(res.get(),float_complext(1,0));
		//auto res= boost::make_shared<cuNDArray<float_complext>>(csm_dims);
		//fill(res.get(),float_complext(1,0));
		return boost::make_shared<cuNDArray<float_complext>>(res);

	}

}

std::tuple<boost::shared_ptr<hoNDArray<floatd2 > >, boost::shared_ptr<hoNDArray<float >>> gpuCSICoilEstimationGadget::separate_traj_and_dcw(
		hoNDArray<float >* traj_dcw) {
	std::vector<size_t> dims = *traj_dcw->get_dimensions();
	std::vector<size_t> reduced_dims(dims.begin()+1,dims.end()); //Copy vector, but leave out first dim
	auto  dcw = boost::make_shared<hoNDArray<float>>(reduced_dims);

	auto traj = boost::make_shared<hoNDArray<floatd2>>(reduced_dims);

	auto dcw_ptr = dcw->get_data_ptr();
	auto traj_ptr = traj->get_data_ptr();
	auto ptr = traj_dcw->get_data_ptr();

	for (size_t i = 0; i < traj_dcw->get_number_of_elements()/3; i++){
		traj_ptr[i][0] = ptr[i*3];
		traj_ptr[i][1] = ptr[i*3+1];
		dcw_ptr[i] = ptr[i*3+2];
	}

	return std::make_tuple(traj,dcw);


}


template<class T> std::tuple<boost::shared_ptr<hoNDArray<T>>, boost::shared_ptr<hoNDArray<T>> > gpuCSICoilEstimationGadget::split_calibration_lines(hoNDArray<T>& data, int ncal_lines, int dim){

	auto dimensions = *data.get_dimensions();

	size_t cal_dim = dimensions[dim];
	if (ncal_lines >= dimensions[dim]){
		throw std::runtime_error("Number of calibration lines is longer or equal to data size. No data left for reconstruction!");
	}
	auto cal_dims = dimensions;
	cal_dims[dim] = ncal_lines;
	dimensions[dim] -= ncal_lines;

	auto calibration_data = boost::make_shared<hoNDArray<T>>(cal_dims);
	auto recon_data = boost::make_shared<hoNDArray<T>>(dimensions);

	size_t smaller_dims = 1;
	for (auto i = 0u; i < dim;i++)
		smaller_dims *= dimensions[i];

	size_t larger_dims = 1;
	for (auto i = dim+1; i < dimensions.size(); i++)
		larger_dims *= dimensions[i];


	auto cal_ptr = calibration_data->get_data_ptr();
	auto recon_ptr = recon_data->get_data_ptr();
	auto data_ptr = data.get_data_ptr();

	for (auto k = 0u; k < larger_dims; k++){
		for (auto n = 0u; n < ncal_lines; n++){
			for (auto m = 0u; m < smaller_dims; m++){
				cal_ptr[m+n*smaller_dims+k*ncal_lines*smaller_dims] = data_ptr[m+n*smaller_dims+k*smaller_dims*cal_dim];
			}
		}
		for (auto n = 0u; n < dimensions[dim]; n++){
			for (auto m = 0u; m < smaller_dims; m++){
				recon_ptr[m+n*smaller_dims+k*dimensions[dim]*smaller_dims] = data_ptr[m+(n+ncal_lines)*smaller_dims+k*smaller_dims*cal_dim];
			}
		}
	}


	return std::make_tuple(calibration_data,recon_data);




}


GADGET_FACTORY_DECLARE(gpuCSICoilEstimationGadget)

} /* namespace Gadgetron */
