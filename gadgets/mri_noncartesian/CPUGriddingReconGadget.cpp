#include "CPUGriddingReconGadget.h"
#include "hoNFFT.h"
#include "vector_td_utilities.h"
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include "hoNDArray_math.h"
#include "b1_map.h"

#include <numeric>
#include <random>

namespace Gadgetron{
	CPUGriddingReconGadget::CPUGriddingReconGadget() : BaseClass(){}
	CPUGriddingReconGadget::~CPUGriddingReconGadget() {}

	int CPUGriddingReconGadget::process_config(ACE_Message_Block* mb){
		GADGET_CHECK_RETURN(BaseClass::process_config(mb) == GADGET_OK, GADGET_FAIL);
		ISMRMRD::IsmrmrdHeader h;
		try{
			deserialize(mb->rd_ptr(), h);
		}catch(...){
			GDEBUG("Error parsing ISMRMRD Header");
		}
		auto matrixsize = h.encoding.front().encodedSpace.matrixSize;
		kernel_width_ = kernel_width.value();
		oversampling_factor_ = gridding_oversampling_factor.value();
		
		image_dims_.push_back(matrixsize.x);
		image_dims_.push_back(matrixsize.y);

		image_dims_os_ = uint64d2(
			static_cast<size_t>(std::ceil(image_dims_[0]*oversampling_factor_)),
			static_cast<size_t>(std::ceil(image_dims_[1]*oversampling_factor_))
		);
		return GADGET_OK;
	}

	int CPUGriddingReconGadget::process(Gadgetron::GadgetContainerMessage<IsmrmrdReconData> *m1){
		IsmrmrdReconData* recon_bit_ = m1->getObjectPtr();
		if(recon_bit_->rbit_.size() > num_encoding_spaces_){
			GWARN_STREAM("Incoming recon+bit has more encoding spaces than the protocol");
		}
		for(size_t e = 0; e < recon_bit_->rbit_.size(); e++){
			IsmrmrdDataBuffered *buffer = &(recon_bit_->rbit_[e].data_);
			IsmrmrdImageArray imarray;

			size_t RO = buffer->data_.get_size(0);
			size_t E1 = buffer->data_.get_size(1);
			size_t E2 = buffer->data_.get_size(2);
			size_t CHA = buffer->data_.get_size(3);
			size_t N = buffer->data_.get_size(4);
			size_t S = buffer->data_.get_size(5);
			size_t SLC = buffer->data_.get_size(6);
			
			imarray.data_.create(image_dims_[0], image_dims_[1], 1, 1, N, S, SLC);
			std::vector<size_t> new_order = {0, 1, 2, 4, 5, 6, 3};
			boost::shared_ptr<hoNDArray<float>> dcw;
			boost::shared_ptr<hoNDArray<floatd2>> traj;

			auto & trajectory = *buffer->trajectory_;
			if(buffer->headers_[0].trajectory_dimensions == 3){
				auto traj_dcw = separate_traj_and_dcw(&trajectory);
				dcw = boost::make_shared<hoNDArray<float>>(std::get<1>(traj_dcw).get());
				traj = boost::make_shared<hoNDArray<floatd2>>(std::get<0>(traj_dcw).get());
			}else if(buffer->headers_[0].trajectory_dimensions == 2){
				auto old_traj_dims = *trajectory.get_dimensions();
				std::vector<size_t> traj_dims (old_traj_dims.begin()+1, old_traj_dims.end());
				hoNDArray<floatd2> tmp_traj(traj_dims, (floatd2*) trajectory.get_data_ptr());
				traj = boost::make_shared<hoNDArray<floatd2>>(tmp_traj);
			}else{
				throw std::runtime_error("unsupported number of trajectory dimensions");
			}

			auto permuted = permute((hoNDArray<float_complext>*)&buffer->data_, &new_order);
			hoNDArray<float_complext> data(*permuted);

			if(dcw){
				float scale_factor = float(prod(image_dims_os_))/asum(dcw.get());
				*dcw *= scale_factor;
			}

			auto images = reconstruct(&data, traj.get(), dcw.get(), CHA);
			// auto csm = estimate_b1_map<float,2>(images.get());
			// *images *= *conj(csm.get());
			auto combined  = sum(images.get(), images->get_number_of_dimensions()-1);
			auto host_img = combined;


			auto elements = imarray.data_.get_number_of_elements();
			memcpy(imarray.data_.get_data_ptr(), host_img->get_data_ptr(), sizeof(float)*2*elements);

			this->compute_image_header(recon_bit_->rbit_[e], imarray, e);
			this->send_out_image_array(recon_bit_->rbit_[e], imarray, e, ((int)e+1), GADGETRON_IMAGE_REGULAR);

			if(replicas.value() > 0 && snr_frame.value() == process_called_times_){
				hoNDArray<std::complex<float>> rep_array(image_dims_[0], image_dims_[1], replicas.value());
				std::mt19937 engine;
				std::normal_distribution<float> distribution;
				for(size_t r = 0; r < replicas.value(); ++r){
					hoNDArray<std::complex<float>> dtmp = buffer->data_;
					auto permuted_rep = permute((hoNDArray<float_complext>*)&dtmp, &new_order);
					auto dataptr = permuted_rep->get_data_ptr();
					for(size_t k = 0; k < permuted_rep->get_number_of_elements(); k++){
						dataptr[k] += std::complex<float>(distribution(engine), distribution(engine));
					}
					hoNDArray<float_complext> data_rep(*permuted_rep);
					images = reconstruct(&data_rep, traj.get(), dcw.get(), CHA);
					// *images *= *conj(csm.get());
					auto combined = sum(images.get(), images->get_number_of_dimensions()-1);
					auto host_img = combined;
					auto elements = imarray.data_.get_number_of_elements();
					size_t offset = image_dims_[0]*image_dims_[1]*r;

					memcpy(rep_array.get_data_ptr()+offset, host_img->get_data_ptr(), sizeof(float)*2*elements);
				}
				
				hoNDArray<float> msg(rep_array.get_dimensions());
				hoNDArray<float> mean(image_dims_[0], image_dims_[1]);
				hoNDArray<float> std(image_dims_[0], image_dims_[1]);

				Gadgetron::abs(rep_array, msg);
				Gadgetron::sum_over_dimension(msg, mean, 2);
				Gadgetron::scal(1.0f/replicas.value(), mean);

				msg -= mean;
				msg *= msg;

				Gadgetron::sum_over_dimension(msg, std, 2);
				Gadgetron::scal(1.0f/(replicas.value()-1),std);
				Gadgetron::sqrt_inplace(&std);

				mean /= std;

				imarray.data_ = *real_to_complex<std::complex<float>>(&mean);
				this->compute_image_header(recon_bit_->rbit_[e], imarray, e);
				this->send_out_image_array(recon_bit_->rbit_[e], imarray, e, image_series.value() + 100 * ((int)e+3),GADGETRON_IMAGE_SNR_MAP);
			}
		}
		m1->release();
		return GADGET_OK;
	}


	boost::shared_ptr<hoNDArray<float_complext>> CPUGriddingReconGadget::reconstruct(
		hoNDArray<float_complext> *data,
		hoNDArray<floatd2> *traj,
		hoNDArray<float> *dcw,
		size_t ncoils
	){
		if(!iterate.value() && dcw){
			hoNFFT_plan<float, 2> plan(from_std_vector<size_t,2>(image_dims_),image_dims_os_,kernel_width_);
			std::vector<size_t> recon_dims = image_dims_;
			recon_dims.push_back(ncoils);
			auto result = new hoNDArray<float_complext>(recon_dims);

			std::vector<size_t> flat_dims = {traj->get_number_of_elements()};
			hoNDArray<floatd2> flat_traj(flat_dims, traj->get_data_ptr());

			plan.preprocess(&flat_traj, hoNFFT_plan<float,2>::NFFT_PREP_NC2C);
			plan.compute(data, result, dcw, hoNFFT_plan<float,2>::NFFT_BACKWARDS_NC2C);
			return boost::shared_ptr<hoNDArray<float_complext>>(result);
		}else{
			// todo 
		}
	}

	std::tuple<boost::shared_ptr<hoNDArray<floatd2>>, boost::shared_ptr<hoNDArray<float>>>
	CPUGriddingReconGadget::separate_traj_and_dcw(
		hoNDArray<float> *traj_dcw
	){
		std::vector<size_t> dims = *traj_dcw->get_dimensions();
		std::vector<size_t> reduced_dims(dims.begin()+1, dims.end());
		auto dcw = boost::make_shared<hoNDArray<float>>(reduced_dims);
		auto traj = boost::make_shared<hoNDArray<floatd2>>(reduced_dims);
		auto dcw_ptr = dcw->get_data_ptr();
		auto traj_ptr = traj->get_data_ptr();
		auto ptr = traj_dcw->get_data_ptr();
		for(size_t i = 0; i < traj_dcw->get_number_of_elements()/3; i++){
			traj_ptr[i][0] = ptr[i*3];
			traj_ptr[i][1] = ptr[i*3+1];
			dcw_ptr[i] = ptr[i*3+2];
		}
		return std::make_tuple(traj,dcw);
	}

	GADGET_FACTORY_DECLARE(CPUGriddingReconGadget);
}

