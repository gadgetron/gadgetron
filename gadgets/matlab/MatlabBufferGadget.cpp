#include "MatlabBufferGadget.h"
#include "MatlabUtils.h"

std::mutex mutex_;

namespace Gadgetron{

int MatlabBufferGadget::process(GadgetContainerMessage<IsmrmrdReconData>* m1)
{
    GDEBUG("Starting MatlabBufferGadget::process\n");
    
    std::lock_guard<std::mutex> lock(mutex_);   

	// Initialize a string for matlab commands
	std::string cmd;

	auto recon_data = m1->getObjectPtr();
	mwSize nencoding_spaces = recon_data->rbit_.size();

	const char* fieldnames[2] = {"data","reference"};
	auto reconArray = mxCreateStructArray(1,&nencoding_spaces,2,fieldnames); // what is this mysterious encoding_spaces ?

    ///////////////////////////////////
                 
    // send the structure to matlab with the data
    for (int i = 0; i <  recon_data->rbit_.size(); i++)
    {
        auto mxrecon = BufferToMatlabStruct(&recon_data->rbit_[i].data_);
        mxSetField(reconArray,i,"data",mxrecon);
        if (recon_data->rbit_[i].ref_)
        {
            auto mxref = BufferToMatlabStruct(recon_data->rbit_[i].ref_.get_ptr());
            mxSetField(reconArray,i,"reference",mxref);
        }
    }
    engPutVariable(engine_, "recon_data", reconArray);
    
    
    /*
    // 2e9 bytes data is the published (as of 2017a) hardcoded limit that MALTAB can load
    // Empirically, it seems that variables up to 2^32 bytes (~4.3 GB) can be loaded.
    // Above that, the error message displayed is: Error using load: cannot read file stdio.
    size_t max_data_size = 2e9;
    
    // compute the data size in bytes
    size_t data_bytes = 0;
    for (int i=0; i < recon_data->rbit_.size(); i++)
    {
        data_bytes += recon_data->rbit_[i].data_.data_.get_number_of_bytes();
        if (recon_data->rbit_[i].ref_)
            data_bytes += recon_data->rbit_[i].ref_.get_ptr()->data_.get_number_of_bytes();
    }
    
    GDEBUG("Bucket size: %lu bytes\n", data_bytes);
    
    bool split_data = data_bytes >= max_data_size;
    
    // send the structure to matlab, with or without the data
    for (int i = 0; i <  recon_data->rbit_.size(); i++)
    {
        auto mxrecon = BufferToMatlabStruct(&recon_data->rbit_[i].data_, split_data);
        mxSetField(reconArray,i,"data",mxrecon);
        if (recon_data->rbit_[i].ref_)
        {
            auto mxref = BufferToMatlabStruct(recon_data->rbit_[i].ref_.get_ptr(), split_data);
            mxSetField(reconArray,i,"reference",mxref);
        }
    }
    engPutVariable(engine_, "recon_data", reconArray);
    
    if(split_data)
    {
        // the dataset needs to be sent in multiple packets
        // The algorithm here splits the multidimensional arrays (data.data
        // and reference.data) into n_packets in the RO (1st) dimension. After all
        // packets are sent, MATLAB reconcatenates everything.
        
        int n_packets = ceil( float(data_bytes) / float(max_data_size) );
        
        GDEBUG("Bucket size limit reached, parsing it into %i packets.\n", n_packets);
        
        for (int i = 0; i < recon_data->rbit_.size(); i++)
        {
            // Allocate memory in MATLAB for faster concatenation
            // Extra dimensions = 1 are automatically discarded in MATLAB
            // command sent: "recon_data(i).data.data = zeros(dim1,dim2,...,dimn)";
            size_t n_dims = recon_data->rbit_[i].data_.data_.get_number_of_dimensions();
            std::string allocate_cmd = "recon_data(" + std::to_string(i+1) + ").data.data = zeros(";
            for(size_t j = 0; j < n_dims; j++)
                allocate_cmd += std::to_string(recon_data->rbit_[i].data_.data_.get_size(j)) + ((j == n_dims - 1 ) ? ");" : ", ");
            std::string dbstring_mcmd1 = allocate_cmd + "\n"; GDEBUG(dbstring_mcmd1.c_str());
            send_matlab_command(allocate_cmd);
            
            if (recon_data->rbit_[i].ref_)
            {
                size_t n_dims_ref = recon_data->rbit_[i].ref_.get_ptr()->data_.get_number_of_dimensions();
                std::string allocate_cmd_ref = "recon_data(" + std::to_string(i+1) + ").reference.data = zeros(";
                for(size_t j = 0; j < n_dims_ref; j++)
                    allocate_cmd_ref += std::to_string(recon_data->rbit_[i].ref_.get_ptr()->data_.get_size(j)) + ((j == n_dims_ref - 1 ) ? ");" : ", ");
                std::string dbstring_mcmd1_ref = allocate_cmd_ref + "\n"; GDEBUG(dbstring_mcmd1_ref.c_str());
                send_matlab_command(allocate_cmd_ref);
            }
            
            GDEBUG("Starting to process packets for data index %i:\n", i+1);
            
            float step = float(recon_data->rbit_[i].data_.data_.get_size(0))/float(n_packets); // step MUST be in float, because the decimal information is used
            for(int p = 0; p < n_packets; p++)
            {
                // (RO) indexes of data to be split
                size_t beg = roundf(float(p  )*step       );
                size_t end = roundf(float(p+1)*step - 1.0f);
                
                // get the split data
                GDEBUG("Creating data packet #%i: from index %lu to %lu...\n", p+1, beg, end);
                mxArray* mxdata = GetSplitReconData(&recon_data->rbit_[i].data_, beg, end);
                
                // send it to MATLAB
                GDEBUG("Sending data packet #%i...\n", p+1);
                std::string packet_name = "data_" + std::to_string(i) + "_" + std::to_string(p);
                engPutVariable(engine_, packet_name.c_str(), mxdata);
                
                // tell MATLAB to concatenate it to recon_data and clear the packet
                // command sent: "recon_data(i).data.data(beg:end,:,:,:,:,:,:) = data_i_p; clear data_i_p;
                std::string concat_cmd = "recon_data(" + std::to_string(i+1) + ").data.data(" + 
                                         std::to_string(beg+1) + ":" + std::to_string(end+1) + 
                                         ",:,:,:,:,:,:) = " + packet_name + "; " +
                                         "clear " + packet_name + ";";
                std::string dbstring_mcmd2 = concat_cmd + "\n"; GDEBUG(dbstring_mcmd2.c_str());
                send_matlab_command(concat_cmd);
                mxDestroyArray(mxdata);
                
                
                if (recon_data->rbit_[i].ref_)
                {
                    // get the split reference data
                    GDEBUG("Creating ref packet #%i: from index %lu to %lu...\n", p+1, beg, end);
                    mxArray* mxdata_ref = GetSplitReconData(recon_data->rbit_[i].ref_.get_ptr(), beg, end);

                    // send it to MATLAB
                    GDEBUG("Sending ref packet #%i...\n", p+1);
                    std::string packet_name_ref = "ref_" + std::to_string(i) + "_" + std::to_string(p);
                    engPutVariable(engine_, packet_name_ref.c_str(), mxdata_ref);

                    // tell MATLAB to concatenate it to recon_data and clear the packet
                    // command sent: "recon_data(i).reference.data(beg:end,:,:,:,:,:,:) = ref_i_p; clear ref_i_p;
                    std::string concat_cmd_ref = "recon_data(" + std::to_string(i+1) + ").reference.data(" + 
                                             std::to_string(beg+1) + ":" + std::to_string(end+1) + 
                                             ",:,:,:,:,:,:) = " + packet_name_ref + "; " +
                                             "clear " + packet_name_ref + ";";
                    std::string dbstring_mcmd2_ref = concat_cmd_ref + "\n"; GDEBUG(dbstring_mcmd2_ref.c_str());
                    send_matlab_command(concat_cmd_ref);

                    mxDestroyArray(mxdata_ref);
                }
            }
        }
    }
    */
    
    GDEBUG("Sending cmd...\n");
    cmd = "[imageQ,bufferQ] = matgadget.run_process(recon_data); matgadget.emptyQ();";
    send_matlab_command(cmd);
    GDEBUG("done.\n");

	// Get the size of the gadget's queue
	mxArray *imageQ = engGetVariable(engine_, "imageQ");
	if (imageQ == NULL) {
		GERROR("Failed to get the imageQ from matgadget\n");
		return GADGET_FAIL;
	}

	size_t qlen = mxGetNumberOfElements(imageQ);
	if (debug_mode_) {
	GDEBUG("Image Queue size: %d \n", qlen);
	}
	const mwSize* dims = mxGetDimensions(imageQ);
	mwSize ndims = mxGetNumberOfDimensions(imageQ);

	GDEBUG("Number of ndims %i \n",ndims);

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
		auto dims = *image.get_dimensions();

		m3->cont(m4);
		if (this->next()->putq(m3) < 0) {
			GDEBUG("Failed to put Image message on queue\n");
			return GADGET_FAIL;
		}

	}
	//Match engGetVariable with mxDestroy___s


	mxArray* bufferQ = engGetVariable(engine_,"bufferQ");

	qlen = mxGetNumberOfElements(bufferQ);
	if (debug_mode_) {
		GDEBUG("Buffer Queue size: %d \n", qlen);
		}

	for (mwIndex idx = 0; idx <qlen; idx++){

		IsmrmrdReconData output_data;
		IsmrmrdReconBit bit;
		bit.data_ = MatlabStructToBuffer(mxGetField(bufferQ,idx,"data"));

		auto ref = mxGetField(bufferQ,idx,"reference");
		if (ref){
			GDEBUG("Adding reference");
			bit.ref_ = MatlabStructToBuffer(ref);
		}
		output_data.rbit_.push_back(bit);
		auto m3 = new GadgetContainerMessage<IsmrmrdReconData>(output_data);
		if (this->next()->putq(m3) < 0){
			GDEBUG("Failed to put Buffer message on queue\n");
			return GADGET_FAIL;
		}

	}

	mxDestroyArray(bufferQ);
	mxDestroyArray(imageQ);
	mxDestroyArray(reconArray); //We're not supposed to delete this? //LA: apparent memory leak if not done

	// We are finished with the incoming messages m1 and m2
	m1->release();

	return GADGET_OK;
}
}


