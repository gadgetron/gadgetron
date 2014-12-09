#include "MatlabBufferGadget.h"
#include "MatlabUtils.h"

namespace Gadgetron{

int MatlabBufferGadget::process(GadgetContainerMessage<IsmrmrdReconData>* m1)
{
	// Initialize a string for matlab commands
	std::string cmd;

	auto recon_data = m1->getObjectPtr();
	mwSize nencoding_spaces = recon_data->rbit_.size();


	auto reconArray = mxCreateCellArray(1,&nencoding_spaces);

	for (int i = 0; i <  recon_data->rbit_.size(); i++){
		auto mxrecon = BufferToMatlabStruct(&recon_data->rbit_[i].data_);
		mxSetCell(reconArray,i,mxrecon);

	}
	engPutVariable(engine_, "recon_data", reconArray);

	cmd = "[imageQ,bufferQ] = matgadget.run_process(recon_data); matgadget.emptyQ(); whos()";
	send_matlab_command(cmd);

	// Get the size of the gadget's queue

	mxArray *imageQ = engGetVariable(engine_, "imageQ");
	if (imageQ == NULL) {
		GADGET_DEBUG1("Failed to get the imageQ from matgadget\n");
		return GADGET_FAIL;
	}

	size_t qlen = mxGetNumberOfElements(imageQ);
	GADGET_DEBUG2("Image Queue size: %d \n", qlen);

	const mwSize* dims = mxGetDimensions(imageQ);
	mwSize ndims = mxGetNumberOfDimensions(imageQ);



	//Read all Image bytes
	for (mwIndex idx = 0; idx < qlen; idx++) {
		mxArray *res_hdr  = mxGetField(imageQ, idx, "bytes");
		mxArray *res_data = mxGetField(imageQ, idx, "image");

		GadgetContainerMessage<ISMRMRD::ImageHeader>* m3 =
				new GadgetContainerMessage<ISMRMRD::ImageHeader>();
		ISMRMRD::ImageHeader *hdr_new = m3->getObjectPtr();
		memcpy(hdr_new, mxGetData(res_hdr), sizeof(ISMRMRD::ImageHeader));

		auto image= MatlabToHoNDArray<std::complex<float>>(res_data);
		auto m4 = new GadgetContainerMessage< hoNDArray< std::complex<float> > >(image);
		auto dims = *image->get_dimensions();

		delete image;
		m3->cont(m4);
		if (this->next()->putq(m3) < 0) {
			GADGET_DEBUG1("Failed to put Image message on queue\n");
			return GADGET_FAIL;
		}

	}
	//Match engGetVariable with mxDestroy___s


	mxArray* bufferQ = engGetVariable(engine_,"bufferQ");

	qlen = mxGetNumberOfElements(bufferQ);
	GADGET_DEBUG2("Buffer Queue size: %d \n", qlen);

	IsmrmrdReconData output_data;
	for (mwIndex idx = 0; idx <qlen; idx++){

		IsmrmrdReconBit bit;
		bit.data_ = MatlabStructToBuffer(mxGetCell(bufferQ,idx));
		output_data.rbit_.push_back(bit);
	}

	if (!output_data.rbit_.empty()){
		auto m3 = new GadgetContainerMessage<IsmrmrdReconData>(output_data.rbit_);
		if (this->next()->putq(m3) < 0){
			GADGET_DEBUG1("Failed to put Buffer message on queue\n");
			return GADGET_FAIL;
		}
	}
	mxDestroyArray(bufferQ);
	mxDestroyArray(imageQ);
	//mxDestroyArray(reconArray); //We're not supposed to delete this?

	// We are finished with the incoming messages m1 and m2
	m1->release();

	return GADGET_OK;
}
}


