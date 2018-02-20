#include "CPUGriddingReconGadget.h"
#include "mri_core_grappa.h"
#include "vector_td_utilities.h"
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include "hoNFFT.h"
#include "hoNDArray.h"
#include "hoNDArray_elemwise.h"
#include "hoNDArray_math.h"
#include "hoCgSolver.h"
#include <time.h>
#include <numeric>

namespace Gadgetron{
	CPUGriddingReconGadget::CPUGriddingReconGadget(){}

	CPUGriddingReconGadget::~CPUGriddingReconGadget(){}

	int CPUGriddingReconGadget::process_config(ACE_Message_Block *mb){
		GADGET_CHECK_RETURN(GenericReconGadget::process_config(mb) == GADGET_OK, GADGET_FAIL);
		ISMRMRD::IsmrmrdHeader h;
		deserialize(mb->rd_ptr(), h);
		auto matrixSize = h.encoding.front().encodedSpace.matrixSize;

		kernelWidth = kernelWidthProperty.value();
		oversamplingFactor = oversamplingFactorProperty.value();

		imageDims.push_back(matrixSize.x); 
		imageDims.push_back(matrixSize.y);
		
		imageDimsOs.push_back(matrixSize.x*oversamplingFactor);
		imageDimsOs.push_back(matrixSize.y*oversamplingFactor);

		return GADGET_OK;
	}

	int CPUGriddingReconGadget::process(GadgetContainerMessage<IsmrmrdReconData> *m1){
		IsmrmrdReconData *recon_bit_ = m1->getObjectPtr();
		process_called_times_++;
		for(size_t e = 0; e < recon_bit_->rbit_.size(); e++){
			IsmrmrdDataBuffered* buffer = &(recon_bit_->rbit_[e].data_);
			IsmrmrdImageArray imarray;
			
			size_t RO = buffer->data_.get_size(0);
			size_t E1 = buffer->data_.get_size(1);
			size_t E2 = buffer->data_.get_size(2);
			size_t CHA = buffer->data_.get_size(3);
			size_t N = buffer->data_.get_size(4);
			size_t S = buffer->data_.get_size(5);
			size_t SLC = buffer->data_.get_size(6);

			imarray.data_.create(imageDims[0], imageDims[1], 1, 1, N, S, SLC);			

			auto &trajectory = *buffer->trajectory_;
			auto trajDcw = separateDcwAndTraj(&trajectory);

			boost::shared_ptr<hoNDArray<float>> dcw = 
				boost::make_shared<hoNDArray<float>>(std::get<1>(trajDcw).get());
			boost::shared_ptr<hoNDArray<floatd2>> traj = 
				boost::make_shared<hoNDArray<floatd2>>(std::get<0>(trajDcw).get());
			
			std::vector<size_t> newOrder = {0, 1, 2, 4, 5, 6, 3};
			auto permuted = permute((hoNDArray<float_complext>*)&buffer->data_,&newOrder);
			hoNDArray<float_complext> data(*permuted);

			auto image = reconstruct(&data, traj.get(), dcw.get(), CHA);
			auto img = *image;
			hoNDArray<float_complext> finalImage; finalImage.create(imageDims[0], imageDims[1]);
			
			// Crop the image
			size_t halfImageDims = (imageDimsOs[0]-imageDims[0])/2;
			for(size_t i = halfImageDims; i < imageDims[0]+halfImageDims; i++)
				for(size_t j = halfImageDims; j < imageDims[1]+halfImageDims; j++)
					finalImage[(i-halfImageDims)+(j-halfImageDims)*imageDims[0]] = img[i+j*imageDimsOs[0]]; 

			auto elements = imarray.data_.get_number_of_elements();
		 	memcpy(imarray.data_.get_data_ptr(), finalImage.get_data_ptr(), sizeof(float)*2*elements);			
			this->compute_image_header(recon_bit_->rbit_[e], imarray, e);
			this->send_out_image_array(recon_bit_->rbit_[e], imarray, e, ((int)e + 1), GADGETRON_IMAGE_REGULAR);		
		}

		m1->release();

		return GADGET_OK;
	}

	boost::shared_ptr<hoNDArray<float_complext>> CPUGriddingReconGadget::reconstruct(
		hoNDArray<float_complext> *data,
		hoNDArray<floatd2> *traj,
		hoNDArray<float> *dcw,
		size_t nCoils
	){
		hoNDArray<float_complext> arg;
		arg.create(imageDimsOs[0], imageDimsOs[0]);
		for(unsigned int i = 0; i < nCoils; ++i){
			hoNDArray<float_complext> channelData(data->get_number_of_elements()/nCoils);
			std::copy(data->begin()+i*(data->get_number_of_elements()/nCoils), data->begin()+(i+1)*(data->get_number_of_elements()/nCoils), channelData.begin());
			
			hoNDArray<float_complext> channelRecon = *reconstructChannel(&channelData, traj, dcw);
			multiplyConj(channelRecon, channelRecon, channelRecon);
			add(arg, channelRecon, arg);	
		}
		sqrt_inplace(&arg);
		return boost::make_shared<hoNDArray<float_complext>>(arg);
	}	

	boost::shared_ptr<hoNDArray<float_complext>> CPUGriddingReconGadget::reconstructChannel(
		hoNDArray<float_complext> *data,
		hoNDArray<floatd2> *traj,
		hoNDArray<float> *dcw
	){	
		if(!iterateProperty.value()){
			hoNFFT_plan<float, 2> plan(
				from_std_vector<size_t, 2>(imageDims),
				oversamplingFactor,
				kernelWidth
			);	
			hoNDArray<float_complext> result; result.create(imageDimsOs[0], imageDimsOs[1]);
			plan.preprocess(*traj);
			plan.compute(*data, result, *dcw, hoNFFT_plan<float, 2>::NFFT_BACKWARDS_NC2C); 

			return boost::make_shared<hoNDArray<float_complext>>(result);
		}else{	
			// do iterative reconstruction
		}
	}

	std::tuple<boost::shared_ptr<hoNDArray<floatd2>>, boost::shared_ptr<hoNDArray<float>>>
	CPUGriddingReconGadget::separateDcwAndTraj(
		hoNDArray<float> *dcwTraj
	){
		std::vector<size_t> dims = *dcwTraj->get_dimensions();
		std::vector<size_t> reducedDims(dims.begin()+1, dims.end());
		auto dcw = boost::make_shared<hoNDArray<float>>(reducedDims);
		auto traj = boost::make_shared<hoNDArray<floatd2>>(reducedDims);

		auto dcwPtr = dcw->get_data_ptr();
		auto trajPtr = traj->get_data_ptr();
		auto ptr = dcwTraj->get_data_ptr();
		for(unsigned int i = 0; i != dcwTraj->get_number_of_elements()/3; i++){
			trajPtr[i][0] = ptr[i*3];
			trajPtr[i][1] = ptr[i*3+1];
			dcwPtr[i] = ptr[i*3+2];
		}
		return std::make_tuple(traj, dcw);
	}

	GADGET_FACTORY_DECLARE(CPUGriddingReconGadget);
}
