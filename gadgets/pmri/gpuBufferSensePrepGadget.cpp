/*
 * gpuBufferSensePrepGadget.cpp
 *
 *  Created on: Dec 10, 2014
 *      Author: dch
 */

#include "gpuBufferSensePrepGadget.h"
#include <ismrmrd/xml.h>
#include "GenericReconJob.h"
#include "cuNFFT.h"
#include "NFFTOperator.h"
#include "cuNDArray_math.h"
#include "vector_td_utilities.h"
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include "b1_map.h"
#include "cuCgSolver.h"
#include "hoNDArray_math.h"
#include "hoNDArray_utils.h"
#include "cuNDArray_fileio.h"
#include "cudaDeviceManager.h"
#include <numeric>

namespace Gadgetron {

gpuBufferSensePrepGadget::gpuBufferSensePrepGadget() {

}

gpuBufferSensePrepGadget::~gpuBufferSensePrepGadget() {

}

int gpuBufferSensePrepGadget::process_config(ACE_Message_Block* mb) {
	ISMRMRD::IsmrmrdHeader h;
	ISMRMRD::deserialize(mb->rd_ptr(),h);

	auto matrixsize = h.encoding.front().encodedSpace.matrixSize;


	profiles_per_frame_ = profiles_per_frame.value();
	kernel_width_ = kernel_width.value();
	oversampling_factor_ = buffer_convolution_oversampling_factor.value();

	unsigned int warp_size = cudaDeviceManager::Instance()->warp_size();
	image_dims_.push_back(((matrixsize.x+warp_size-1)/warp_size)*warp_size);
	image_dims_.push_back(((matrixsize.y+warp_size-1)/warp_size)*warp_size);

	image_dims_recon_.push_back(((static_cast<size_t>(std::ceil(matrixsize.x*reconstruction_os_factor.value()))+warp_size-1)/warp_size)*warp_size);
	image_dims_recon_.push_back(((static_cast<size_t>(std::ceil(matrixsize.y*reconstruction_os_factor.value()))+warp_size-1)/warp_size)*warp_size);

	image_dims_recon_os_ = uint64d2
			(((static_cast<size_t>(std::ceil(image_dims_recon_[0]*oversampling_factor_))+warp_size-1)/warp_size)*warp_size,
					((static_cast<size_t>(std::ceil(image_dims_recon_[1]*oversampling_factor_))+warp_size-1)/warp_size)*warp_size);

	// In case the warp_size constraint kicked in
	oversampling_factor_ = float(image_dims_recon_os_[0])/float(image_dims_recon_[0]);

	return GADGET_OK;

}

int gpuBufferSensePrepGadget::process(
		GadgetContainerMessage<IsmrmrdReconData>* m1) {

	IsmrmrdReconData* recondata= m1->getObjectPtr();

	if (recondata->rbit_.size() != 1){
		throw std::runtime_error("gpuBufferSensePrepGadget only support a single encoding space");
	}

	IsmrmrdReconBit& reconbit = recondata->rbit_[0];

	GenericReconJob job;

	IsmrmrdDataBuffered* buffer = &reconbit.data_;

	//Use reference data if available.
	if (reconbit.ref_){
		GDEBUG("Using Reference data for CSM estimation\n");
		buffer = &reconbit.ref_.value();
	}

	size_t ncoils = buffer->headers_[0].active_channels;


	std::vector<size_t> new_order = {0,1,2,4,5,6,3};

	boost::shared_ptr<cuNDArray<float>> dcw;
	boost::shared_ptr<cuNDArray<floatd2>> traj;
	if (buffer->trajectory_){
		auto & trajectory = *buffer->trajectory_;

		if (buffer->headers_[0].trajectory_dimensions == 3){
			auto traj_dcw = separate_traj_and_dcw(&trajectory);
			dcw = boost::make_shared<cuNDArray<float>>(std::get<1>(traj_dcw).get());
			traj = boost::make_shared<cuNDArray<floatd2>>(std::get<0>(traj_dcw).get());
		} else if (buffer->headers_[0].trajectory_dimensions == 2){
			auto old_traj_dims = *trajectory.get_dimensions();
			std::vector<size_t> traj_dims (old_traj_dims.begin()+1,old_traj_dims.end()); //Remove first element
			hoNDArray<floatd2> tmp_traj(traj_dims,(floatd2*)trajectory.get_data_ptr());
			traj = boost::make_shared<cuNDArray<floatd2>>(tmp_traj);
		} else {
			throw std::runtime_error("Unsupported number of trajectory dimensions");
		}
	}
	{
		auto tmpdim = *buffer->data_.get_dimensions();
		std::stringstream stream; 
		for (auto dim : tmpdim)
			stream << dim << " ";
		stream << "\n";
		GINFO(stream.str().c_str());
		auto permuted = permute(*(hoNDArray<float_complext>*)&buffer->data_,new_order);
		cuNDArray<float_complext> data(permuted);
		if (dcw){
			float scale_factor = float(prod(image_dims_recon_os_))/asum(dcw.get());
			*dcw *= scale_factor;
		}

		auto reg_images = reconstruct_regularization(&data,traj.get(),dcw.get(),ncoils);
		//reg_images->squeeze();

		auto csm = estimate_b1_map<float,2>(reg_images.get());
		*reg_images *= csm;
		auto combined = sum(reg_images.get(),reg_images->get_number_of_dimensions()-1);

		auto tmp_combined = abs(reg_images.get());
		auto tmpcsm = abs(&csm);
		job.csm_host_ = csm.to_host();
		job.reg_host_ = combined->to_host();
	}


	IsmrmrdDataBuffered* mainbuffer = &reconbit.data_;

	//Permute as Sensegadgets expect last dimension to be coils. *Sigh*
	job.dat_host_ = boost::make_shared<hoNDArray<complext<float>>>(permute(*(hoNDArray<float_complext>*)&mainbuffer->data_,new_order));

	if (mainbuffer->trajectory_){
		auto & trajectory = *mainbuffer->trajectory_;
		if (mainbuffer->headers_[0].trajectory_dimensions >2 ){
			auto traj_dcw = separate_traj_and_dcw(&trajectory);
			job.tra_host_ = std::get<0>(traj_dcw);
			job.dcw_host_ = std::get<1>(traj_dcw);
		} else if (mainbuffer->headers_[0].trajectory_dimensions == 2){
			auto old_traj_dims = *trajectory.get_dimensions();
			std::vector<size_t> traj_dims (old_traj_dims.begin()+1,old_traj_dims.end()); //Remove first element
			hoNDArray<floatd2> tmp_traj(traj_dims,(floatd2*)trajectory.get_data_ptr());
			job.tra_host_ = boost::make_shared<hoNDArray<floatd2>>(tmp_traj);
			auto host_dcw = boost::make_shared<hoNDArray<float>>(traj_dims);
			fill(host_dcw.get(),1.0f);
			job.dcw_host_ = host_dcw;

		} else {
			throw std::runtime_error("Unsupported number of trajectory dimensions");
		}
	}
	{
		float scale_factor = float(prod(image_dims_recon_os_))/asum(job.dcw_host_.get());
		*job.dcw_host_  *= scale_factor;
	}

	auto data_dims = *job.dat_host_->get_dimensions();
	//Sense gadgets expect only 1 dimension for encoding, so collapse the first
	size_t elements = std::accumulate(data_dims.begin(),data_dims.end()-1,1,std::multiplies<size_t>());
	std::vector<size_t> new_data_dims = {elements,data_dims.back()};
	job.dat_host_->reshape(&new_data_dims);

	size_t traj_elements = job.tra_host_->get_number_of_elements();
	auto traj_dims = *job.tra_host_->get_dimensions();

	size_t kpoints_per_frame = traj_dims[0]*profiles_per_frame_;
	if (traj_elements%kpoints_per_frame){
		std::stringstream ss;
		ss << "Profiles per frame (" << profiles_per_frame_ << ") must be a divisor of total number of profiles (" << traj_elements/traj_dims[0] << ")";
		throw std::runtime_error(ss.str());
	}
	std::vector<size_t> new_traj_dims ={kpoints_per_frame,traj_elements/kpoints_per_frame};

	job.tra_host_->reshape(&new_traj_dims);
	job.dcw_host_->reshape(&new_traj_dims);


	//Let's invent some image headers!
	size_t total_frames = profiles_per_frame_ > 0 ? mainbuffer->headers_.get_number_of_elements()/profiles_per_frame_ : 1 ;
	job.image_headers_ = boost::shared_array<ISMRMRD::ImageHeader>(new ISMRMRD::ImageHeader[total_frames]);
	for (size_t i = 0; i < total_frames; i++){
		job.image_headers_[i] = create_image_header(mainbuffer->headers_[i*profiles_per_frame_],mainbuffer->sampling_,i,total_frames);
	}


	m1->release(); //We be done with everything now.

	auto header_message = new GadgetContainerMessage<ISMRMRD::ImageHeader>(job.image_headers_[0]);

	auto job_message = new GadgetContainerMessage<GenericReconJob>(job);

	header_message->cont(job_message);

	if (!this->next()->putq(header_message)){
		GDEBUG("Failed to put message on que");
		return GADGET_FAIL;
	} else
		return GADGET_OK;



	//cuNDArray<float_complext> reg_images = reconstruct_regularization(reconbit.data_);
}

boost::shared_ptr<cuNDArray<float_complext> > gpuBufferSensePrepGadget::reconstruct_regularization(
		cuNDArray<float_complext>* data, cuNDArray<floatd2>* traj, cuNDArray<float>* dcw, size_t ncoils ) {

	if (dcw) { //We have density compensation, so we can get away with gridding

		cuNFFT_impl<float,2> plan(from_std_vector<size_t,2>(image_dims_recon_),image_dims_recon_os_,kernel_width_);
		std::vector<size_t> csm_dims = image_dims_recon_;
		csm_dims.push_back(ncoils);
		auto result = new cuNDArray<float_complext>(csm_dims);
		GDEBUG("Coils %i \n\n",ncoils);

		std::vector<size_t> flat_dims = {traj->get_number_of_elements()};
		cuNDArray<floatd2> flat_traj(flat_dims,traj->get_data_ptr());
		GDEBUG("traj: %i data %i\n",traj->get_number_of_elements(),data->get_number_of_elements());
		GDEBUG("Preprocessing\n\n");
		plan.preprocess(&flat_traj,NFFT_prep_mode::NC2C);
		GDEBUG("Computing\n\n");
		plan.compute(data,*result,dcw,NFFT_comp_mode::BACKWARDS_NC2C);

		return boost::shared_ptr<cuNDArray<float_complext>>(result);

	} else { //No density compensation, we have to do iterative reconstruction.
		std::vector<size_t> csm_dims = image_dims_recon_;
		csm_dims.push_back(ncoils);

		auto E = boost::make_shared<NFFTOperator<cuNDArray,float,2>>();

		E->setup(from_std_vector<size_t,2>(image_dims_recon_),image_dims_recon_os_,kernel_width_);
		std::vector<size_t> flat_dims = {traj->get_number_of_elements()};
		cuNDArray<floatd2> flat_traj(flat_dims,traj->get_data_ptr());

		E->set_domain_dimensions(&csm_dims);
		cuCgSolver<float_complext> solver;
		solver.set_max_iterations(0);
		solver.set_encoding_operator(E);
		solver.set_output_mode(cuCgSolver<float_complext>::OUTPUT_VERBOSE);
		E->set_codomain_dimensions(data->get_dimensions().get());
		E->preprocess(&flat_traj);
		auto res = solver.solve(data);
		return res;
	}
}

std::tuple<boost::shared_ptr<hoNDArray<floatd2 > >, boost::shared_ptr<hoNDArray<float >>> gpuBufferSensePrepGadget::separate_traj_and_dcw(
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

ISMRMRD::ImageHeader gpuBufferSensePrepGadget::create_image_header(
		ISMRMRD::AcquisitionHeader& base_head, const SamplingDescription& samp, size_t idx, size_t num_frames) {

	ISMRMRD::ImageHeader header;
	header.version = base_head.version;

	header.matrix_size[0] = image_dims_recon_[0];
	header.matrix_size[1] = image_dims_recon_[1];
	header.matrix_size[2] = num_frames;


	header.field_of_view[0] = samp.recon_FOV_[0];
	header.field_of_view[1] = samp.recon_FOV_[1];
	header.field_of_view[2] = samp.recon_FOV_[2];

	header.channels = 1;
	header.slice = base_head.idx.slice;
	header.set = base_head.idx.set;

	header.acquisition_time_stamp = base_head.acquisition_time_stamp;
	memcpy(header.physiology_time_stamp, base_head.physiology_time_stamp, sizeof(uint32_t)*ISMRMRD::ISMRMRD_PHYS_STAMPS);

	memcpy(header.position, base_head.position, sizeof(float)*3);
	memcpy(header.read_dir, base_head.read_dir, sizeof(float)*3);
	memcpy(header.phase_dir, base_head.phase_dir, sizeof(float)*3);
	memcpy(header.slice_dir, base_head.slice_dir, sizeof(float)*3);
	memcpy(header.patient_table_position, base_head.patient_table_position, sizeof(float)*3);

	header.data_type = ISMRMRD::ISMRMRD_CXFLOAT;
	header.image_index = idx;
	header.image_series_index = 0;

	return header;



}

GADGET_FACTORY_DECLARE(gpuBufferSensePrepGadget)

} /* namespace Gadgetron */
