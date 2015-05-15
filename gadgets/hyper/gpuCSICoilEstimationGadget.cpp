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
#include "cuNFFTOperator.h"
#include "cuCgSolver.h"
#include "cuNDArray_fileio.h"
#include "trajectory_utils.h"
#include <boost/make_shared.hpp>

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
		GadgetContainerMessage<IsmrmrdAcquisitionBucket>* m1) {
	IsmrmrdAcquisitionBucket* bucket = m1->getObjectPtr();

	auto cm1 = new GadgetContainerMessage<cuSenseData>();
	auto senseData = cm1->getObjectPtr();

	coils = bucket->data_.front().head_->getObjectPtr()->active_channels;
	GDEBUG("Active channels %i \n",coils);


	{

		hoNDArray<std::complex<float>> * ho_data;
		hoNDArray<float>* ho_traj;

		std::tie(ho_data,ho_traj) = combine_data(bucket->data_);

		if (skip_lines_ > 0){
			auto cal_dims = *ho_data->get_dimensions();
			cal_dims.back() = skip_lines_;
			auto data_dims = *ho_data->get_dimensions();
			data_dims.back() -= skip_lines_;


			hoNDArray<float_complext> cal_view(cal_dims,(float_complext*) ho_data->get_data_ptr());
			senseData->freq_calibration = boost::make_shared<cuNDArray<float_complext>>(cal_view);
			senseData->freq_calibration->squeeze();
			hoNDArray<float_complext> data_view(data_dims,(float_complext*)ho_data->get_data_ptr()+cal_view.get_number_of_elements());
			senseData->data = boost::make_shared<cuNDArray<float_complext>>(data_view);
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

		delete ho_data;
		delete ho_traj;
	}


	//Remove Initial Spirals

	boost::shared_ptr< cuNDArray<float_complext> >  ref_data;
	boost::shared_ptr< cuNDArray<floatd2> > ref_traj;
	boost::shared_ptr<cuNDArray<float> > ref_dcw;


	if (bucket->ref_.empty()){
		ref_data = senseData->data;
		ref_traj = senseData->traj;
		ref_dcw = senseData->dcw;
	} else {

		hoNDArray<std::complex<float>> * ho_data;
		hoNDArray<float>* ho_traj;
		std::tie(ho_data,ho_traj) = combine_data(bucket->ref_);

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
		delete ho_data;
		delete ho_traj;


	}

	senseData->csm = calculate_CSM(ref_data.get(),ref_traj.get(),ref_dcw.get());


	if (this->next()->putq(cm1) == GADGET_FAIL){
		GERROR("Failed to put message on que\n");
		return GADGET_FAIL;
	}

	return GADGET_OK;


}

std::tuple<hoNDArray<std::complex<float>>*, hoNDArray<float>*> gpuCSICoilEstimationGadget::combine_data(
		std::vector<IsmrmrdAcquisitionData>& acquisitions) {


	std::vector<size_t> data_dims = *acquisitions.front().data_->getObjectPtr()->get_dimensions();
	std::vector<size_t> traj_dims = *acquisitions.front().traj_->getObjectPtr()->get_dimensions();
	std::vector<size_t> base_dim = data_dims;
	data_dims.push_back(acquisitions.size());
	if (acquisitions.size() == 1 ||  acquisitions[2].traj_) //Trajectory present on all acquisitions
		traj_dims.push_back(acquisitions.size());


	auto result = new hoNDArray<std::complex<float> >(data_dims);
	auto traj = new hoNDArray<float>(traj_dims);

	std::complex<float>* ptr = result->get_data_ptr();
	float* traj_ptr = traj->get_data_ptr();
	for (const IsmrmrdAcquisitionData & data : acquisitions){
		hoNDArray<std::complex<float>>* array = data.data_->getObjectPtr();
		if (data.traj_) { //Only copy if trajectory is present
			hoNDArray<float>* array_traj = data.traj_->getObjectPtr();
			memcpy(traj_ptr,array_traj->get_data_ptr(),array_traj->get_number_of_bytes());
			traj_ptr += array_traj->get_number_of_elements();
		}
		if (!array->dimensions_equal(&base_dim)){
			return std::tuple<hoNDArray<std::complex<float>>*, hoNDArray<float>*>(nullptr,nullptr);
		}
		memcpy(ptr,array->get_data_ptr(),array->get_number_of_bytes());
		ptr += array->get_number_of_elements();

	}

	return std::make_tuple(result,traj);

}

boost::shared_ptr<cuNDArray<float_complext> > gpuCSICoilEstimationGadget::calculate_CSM(
		cuNDArray<float_complext>* data, cuNDArray<floatd2>* traj, cuNDArray<float>* dcw ) {


	if (dcw) { //We have density compensation, so we can get away with gridding

		cuNFFT_plan<float,2> plan(from_std_vector<size_t,2>(img_size),from_std_vector<size_t,2>(img_size)*size_t(2),kernel_width_);
		std::vector<size_t> csm_dims = img_size;
		csm_dims.push_back(coils);
		cuNDArray<float_complext> tmp(csm_dims);
		GDEBUG("Coils %i \n\n",tmp.get_size(2));
		std::vector<size_t> flat_dims = {traj->get_number_of_elements()};
		cuNDArray<floatd2> flat_traj(flat_dims,traj->get_data_ptr());

		std::vector<size_t> spiral_dims{data->get_size(0),data->get_size(1)}; //Trajectories, coils
		cuNDArray<complext<float>> second_spiral(spiral_dims,data->get_data_ptr()+spiral_dims[0]*spiral_dims[1]*0);
		std::vector<size_t> spiral_traj_dims{spiral_dims[0]};
		cuNDArray<floatd2> spiral_traj(spiral_traj_dims,traj->get_data_ptr()+spiral_dims[0]*0);
		cuNDArray<float> spiral_dcw(spiral_traj_dims,dcw->get_data_ptr()+spiral_dims[0]*0);

		GDEBUG("Preprocessing\n\n");
		plan.preprocess(&spiral_traj,cuNFFT_plan<float,2>::NFFT_PREP_NC2C);
		GDEBUG("Computing\n\n");
		plan.compute(&second_spiral,&tmp,&spiral_dcw,cuNFFT_plan<float,2>::NFFT_BACKWARDS_NC2C);
		auto tmp_abs = abs(&tmp);

		return estimate_b1_map<float,2>(&tmp);

	} else { //No density compensation, we have to do iterative reconstruction.
		std::vector<size_t> csm_dims = img_size;
		csm_dims.push_back(coils);

		auto E = boost::make_shared<cuNFFTOperator<float,2>>();

		E->setup(from_std_vector<size_t,2>(img_size),from_std_vector<size_t,2>(img_size)*size_t(2),kernel_width_);
		std::vector<size_t> flat_dims = {traj->get_number_of_elements()};
		cuNDArray<floatd2> flat_traj(flat_dims,traj->get_data_ptr());

		E->set_domain_dimensions(&csm_dims);
		cuCgSolver<float_complext> solver;
		solver.set_max_iterations(20);
		solver.set_encoding_operator(E);
		std::vector<size_t> spiral_dims{data->get_size(0),data->get_size(1)}; //Trajectories, coils
		cuNDArray<complext<float>> second_spiral(spiral_dims,data->get_data_ptr()+spiral_dims[0]*spiral_dims[1]*0);
		E->set_codomain_dimensions(&spiral_dims);
		std::vector<size_t> spiral_traj_dims{spiral_dims[0]};
		cuNDArray<floatd2> spiral_traj(spiral_traj_dims,traj->get_data_ptr()+spiral_dims[0]*0);
		E->preprocess(&spiral_traj);
		auto tmp = solver.solve(&second_spiral);
		auto tmp_abs = abs(tmp.get());


		auto res = estimate_b1_map<float,2>(tmp.get());
		//fill(res.get(),float_complext(1,0));
		//auto res= boost::make_shared<cuNDArray<float_complext>>(csm_dims);
		//fill(res.get(),float_complext(1,0));
		return res;

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



GADGET_FACTORY_DECLARE(gpuCSICoilEstimationGadget)

} /* namespace Gadgetron */
