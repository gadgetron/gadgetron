/*
 * gpuCoilEstimationGadget.cpp
 *
 *  Created on: Nov 12, 2014
 *      Author: dch
 */

#include "gpuCoilEstimationGadget.h"
#include <utility>
#include <ismrmrd/xml.h>
#include "cuNDArray.h"
#include "cuNFFT.h"
#include "b1_map.h"
#include "vector_td_utilities.h"
#include "cuNFFTOperator.h"
#include "cuCgSolver.h"
#include "cuNDArray_fileio.h"

namespace Gadgetron {

gpuCoilEstimationGadget::gpuCoilEstimationGadget() {
	set_parameter("kernel_width","5.5");

}

gpuCoilEstimationGadget::~gpuCoilEstimationGadget() {
	// TODO Auto-generated destructor stub
}

int gpuCoilEstimationGadget::process_config(ACE_Message_Block* mb) {
	ISMRMRD::IsmrmrdHeader h;
	ISMRMRD::deserialize(mb->rd_ptr(),h);

	if (h.encoding.size() != 1){
		GADGET_DEBUG1("Coil estimation gadget only supports encoding spaces of size 1\n");
		return GADGET_FAIL;
	}
	img_size = {h.encoding[0].reconSpace.matrixSize.x,h.encoding[0].reconSpace.matrixSize.y};
	kernel_width = get_double_value("kernel_width");

	coils = h.acquisitionSystemInformation->receiverChannels;
}

int gpuCoilEstimationGadget::process(
		GadgetContainerMessage<IsmrmrdAcquisitionBucket>* m1) {
	IsmrmrdAcquisitionBucket* bucket = m1->getObjectPtr();

	auto cm1 = new GadgetContainerMessage<cuSenseData>();
	auto senseData = cm1->getObjectPtr();

	coils = bucket->data_.front().head_->getObjectPtr()->active_channels;
	GADGET_DEBUG2("Active channels %i \n",coils);


	{

		hoNDArray<std::complex<float>> * ho_data;
		hoNDArray<float>* ho_traj;

		std::tie(ho_data,ho_traj) = combine_data(bucket->data_);
		senseData->data = boost::make_shared<cuNDArray<float_complext>>(reinterpret_cast<hoNDArray<float_complext>*>(ho_data));

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


	cuNDArray<float_complext> * ref_data;
	cuNDArray<floatd2>* ref_traj;
	cuNDArray<float>* ref_dcw;


	if (bucket->ref_.empty()){
		ref_data = senseData->data.get();
		ref_traj = senseData->traj.get();
		ref_dcw = senseData->dcw.get();
	} else {

		hoNDArray<std::complex<float>> * ho_data;
		hoNDArray<float>* ho_traj;
		std::tie(ho_data,ho_traj) = combine_data(bucket->ref_);

		ref_data = new cuNDArray<float_complext>(reinterpret_cast<hoNDArray<float_complext>*>(ho_data));
		if (ho_traj->get_size(0) > 2){
			auto traj_dcw = separate_traj_and_dcw(ho_traj);
			ref_traj = new cuNDArray<floatd2>(*std::get<0>(traj_dcw));
			ref_dcw = new cuNDArray<float>(*std::get<1>(traj_dcw));
		} else {
			std::vector<size_t> tdims = *ho_traj->get_dimensions();
			std::vector<size_t> tmp_dim(tdims.begin()+1,tdims.end());
			hoNDArray<floatd2> tmp(tmp_dim,reinterpret_cast<floatd2*>(ho_traj->get_data_ptr()));
			ref_traj = new cuNDArray<floatd2>(tmp);
		}
		delete ho_data;
		delete ho_traj;


	}


	senseData->csm = calculate_CSM(ref_data,ref_traj,ref_dcw);

	this->next()->putq(cm1);
	//All important stuff has been taken from the bucket. Free it.
	m1->release();




}

std::tuple<hoNDArray<std::complex<float>>*, hoNDArray<float>*> gpuCoilEstimationGadget::combine_data(
		std::vector<IsmrmrdAcquisitionData>& acquisitions) {


	std::vector<size_t> data_dims = *acquisitions.front().data_->getObjectPtr()->get_dimensions();
	std::vector<size_t> traj_dims = *acquisitions.front().traj_->getObjectPtr()->get_dimensions();
	std::vector<size_t> base_dim = data_dims;
	data_dims.push_back(acquisitions.size());
	traj_dims.push_back(acquisitions.size());
	auto result = new hoNDArray<std::complex<float> >(data_dims);
	auto traj = new hoNDArray<float>(traj_dims);

	std::complex<float>* ptr = result->get_data_ptr();
	float* traj_ptr = traj->get_data_ptr();
	for (const IsmrmrdAcquisitionData & data : acquisitions){
		hoNDArray<std::complex<float>>* array = data.data_->getObjectPtr();
		hoNDArray<float>* array_traj = data.traj_->getObjectPtr();
		if (!array->dimensions_equal(&base_dim)){
			return std::tuple<hoNDArray<std::complex<float>>*, hoNDArray<float>*>(nullptr,nullptr);
		}
		memcpy(ptr,array->get_data_ptr(),array->get_number_of_bytes());
		ptr += array->get_number_of_elements();
		memcpy(traj_ptr,array_traj->get_data_ptr(),array_traj->get_number_of_bytes());
		traj_ptr += array_traj->get_number_of_elements();

	}

	return std::make_tuple(result,traj);

}

boost::shared_ptr<cuNDArray<float_complext> > gpuCoilEstimationGadget::calculate_CSM(
		cuNDArray<float_complext>* data, cuNDArray<floatd2>* traj, cuNDArray<float>* dcw ) {


	if (dcw) { //We have density compensation, so we can get away with gridding

		cuNFFT_plan<float,2> plan(from_std_vector<size_t,2>(img_size),from_std_vector<size_t,2>(img_size)*size_t(2),kernel_width);
		std::vector<size_t> csm_dims = img_size;
		csm_dims.push_back(coils);
		cuNDArray<float_complext> tmp(csm_dims);
		GADGET_DEBUG2("Coils %i \n",tmp.get_size(2));
		std::vector<size_t> flat_dims = {traj->get_number_of_elements()};
		cuNDArray<floatd2> flat_traj(flat_dims,traj->get_data_ptr());
		plan.preprocess(&flat_traj,cuNFFT_plan<float,2>::NFFT_PREP_NC2C);
		plan.compute(data,&tmp,dcw,cuNFFT_plan<float,2>::NFFT_BACKWARDS_NC2C);
		auto tmp_abs = abs(&tmp);
		write_nd_array(tmp_abs.get(),"images.real");

		return estimate_b1_map<float,2>(&tmp);

	} else { //No density compensation, we have to do iterative reconstruction.
		std::vector<size_t> csm_dims = img_size;
		csm_dims.push_back(coils);

		auto E = boost::make_shared<cuNFFTOperator<float,2>>();
		E->setup(from_std_vector<size_t,2>(img_size),from_std_vector<size_t,2>(img_size)*size_t(2),kernel_width);
		E->preprocess(traj);
		E->set_domain_dimensions(&csm_dims);
		cuCgSolver<float_complext> solver;
		solver.set_max_iterations(20);
		solver.set_encoding_operator(E);
		auto tmp = solver.solve(data);
		return estimate_b1_map<float,2>(tmp.get());


	}

}

std::tuple<boost::shared_ptr<hoNDArray<floatd2 > >, boost::shared_ptr<hoNDArray<float >>> gpuCoilEstimationGadget::separate_traj_and_dcw(
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

GADGET_FACTORY_DECLARE(gpuCoilEstimationGadget)

} /* namespace Gadgetron */
