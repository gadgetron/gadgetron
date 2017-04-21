#include "MatlabBufferGadget.h"
#include "MatlabUtils.h"

std::mutex mutex_;

namespace Gadgetron{

int MatlabBufferGadget::process(GadgetContainerMessage<IsmrmrdReconData>* m1)
{
    std::lock_guard<std::mutex> lock(mutex_);   

	// Initialize a string for matlab commands
	std::string cmd;

	auto recon_data = m1->getObjectPtr();
	mwSize nencoding_spaces = recon_data->rbit_.size();

	const char* fieldnames[2] = {"data","reference"};
	auto reconArray = mxCreateStructArray(1,&nencoding_spaces,2,fieldnames);
	//auto reconArray = mxCreateCellArray(1,&nencoding_spaces);

    // 2e9 bytes data is the published (as of 2017a) hardcoded limit that engPutVariable can transfer.
    // Empirically, it seems that variables up to 2^32 bytes (~4.3 GB) can be sent.
    size_t max_data_size = 2e9;
    if(sizeof(recon_data->rbit_) < max_data_size) 
    {
        // the dataset is small enough to be sent all at once (original code)
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
        
        GDEBUG("Sending the whole buckets... ");
        engPutVariable(engine_, "recon_data", reconArray);
        GDEBUG("done\n");
    }
    else
    {
        // the dataset needs to be sent in multiple packets
        // The algorithm here splits the multidimensional arrays (data.data
        // and reference.data) into n_packets in the RO dimension. After all
        // packets are sent, MATLAB reconcatenates everything.
        
        int n_packets = ceil( float(sizeof(recon_data->rbit_)) / float(max_data_size) );

        
        GDEBUG("Bucket size limit reached, parsing it into %i packets.\n", n_packets);
        
        for (int i = 0; i <  recon_data->rbit_.size(); i++)
        {            
            // Create the regular MATLAB structure, but omits the data for the fields "data" and "reference".
            auto mxrecon = BufferToMatlabStruct(&recon_data->rbit_[i].data_, true);
            mxSetField(reconArray,i,"data",mxrecon);
            if (recon_data->rbit_[i].ref_)
            {
                auto mxref = BufferToMatlabStruct(recon_data->rbit_[i].ref_.get_ptr(), true);
                mxSetField(reconArray,i,"reference",mxref);
            }
            
            /*
            // send the packets
            size_t n_RO = sizeof(recon_data->rbit_[i].data_.data_) / sizeof(recon_data->rbit_[i].data_.data_[0]);
            float step = float(n_RO)/float(n_packets);
            
            GDEBUG("Starting to process packets for index %i:\n", i+1);
            for(int p = 0; p < n_packets; p++)
            {
                // indexes of data to be split
                size_t beg = roundf(float(p  )*step       );
                size_t end = roundf(float(p+1)*step - 1.0f);
                
                // create the packet. A copy of the data is being done here,
                // which overall increase the RAM usage if packets are needed.
                // There may be a more efficient way to do this.
                GDEBUG("Creating data packet #%i...\n", p+1);
                
                //void *packet = malloc( (end-beg)*sizeof(recon_data->rbit_[i].data_.data_[0]) );
                //std::copy( &(recon_data->rbit_[i].data_.data_[beg]),
                //           &(recon_data->rbit_[i].data_.data_[end]),
                //           &(packet[0]));
                //auto mxdata = hoNDArrayToMatlab(&packet);
                decltype(recon_data->rbit_[i].data_.data_) packet [end-beg]
                            [sizeof(recon_data->rbit_[i].data_.data_[0])                / sizeof(recon_data->rbit_[i].data_.data_[0][0])]
                            [sizeof(recon_data->rbit_[i].data_.data_[0][0])             / sizeof(recon_data->rbit_[i].data_.data_[0][0][0])]
                            [sizeof(recon_data->rbit_[i].data_.data_[0][0][0])          / sizeof(recon_data->rbit_[i].data_.data_[0][0][0][0])]
                            [sizeof(recon_data->rbit_[i].data_.data_[0][0][0][0])       / sizeof(recon_data->rbit_[i].data_.data_[0][0][0][0][0])]
                            [sizeof(recon_data->rbit_[i].data_.data_[0][0][0][0][0])    / sizeof(recon_data->rbit_[i].data_.data_[0][0][0][0][0][0])]
                            [sizeof(recon_data->rbit_[i].data_.data_[0][0][0][0][0][0]) / sizeof(recon_data->rbit_[i].data_.data_[0][0][0][0][0][0][0])]
                
                // convert the packet to MALTAB array and send it
                GDEBUG("Sending data packet #%i...\n", p+1);
                engPutVariable(engine_, "data_" + to_sring(i) + "_" + to_sring(p), mxdata);
                free(packet);
                
                // do the same for the reference
                if (recon_data->rbit_[i].ref_)
                {
                    GDEBUG("Creating reference packet #%i...\n", p+1);
                    // create the ref packet
                    auto packet_ref = malloc( (end-beg)*sizeof(recon_data->rbit_[i].ref_.data_[0]) );
                    std::copy( &(recon_data->rbit_[i].ref_.data_[beg]),
                               &(recon_data->rbit_[i].ref_.data_[end]),
                               &(packet_ref[0]));
                    auto mxdata_ref = hoNDArrayToMatlab(&packet_ref);

                    // convert the ref packet to MALTAB array and send it
                    GDEBUG("Sending reference packet #%i...\n", p+1);
                    engPutVariable(engine_, "ref_" + to_sring(i) + "_" + to_sring(p), mxdata_ref);
                    free(packet_ref);
                }
            }*/
        }
        engPutVariable(engine_, "recon_data", reconArray);
        
        /*
        //send the command to reconcatenate the data and ref
        for (int i = 0; i <  recon_data->rbit_.size(); i++)
        {
            GDEBUG("MATLAB concatenation for index %i...\n", i+1);
            
            // create a concatenation MATLAB command
            string concat_data = "[";
            string concat_ref  = "[";
            for(int p = 0; p < n_packets; p++)
            {
                concat_data += "data_" + to_sring(i) + "_" + to_sring(p) + "; ";
                if (recon_data->rbit_[i].ref_)
                    concat_data += "ref_" + to_sring(i) + "_" + to_sring(p) + "; ";
            }
            
            // send the concatenation command to MATLAB
            send_matlab_command("recon_data.data(" + to_string(i) + ").data  = " + concat_data + "];");
            if (recon_data->rbit_[i].ref_)
                send_matlab_command("recon_data.ref(" + to_string(i) + ").data  = " + concat_ref + "];");
            
            // clear the MATLAB data copies
            for(int p = 0; p < n_packets; p++)
            {
                send_matlab_command("clear " + "data_" + to_sring(i) + "_" + to_sring(p) + "; ");
                if (recon_data->rbit_[i].ref_)
                    send_matlab_command("clear " + "ref_" + to_sring(i) + "_" + to_sring(p) + "; ");
            }
        }        */
    }
    
    cmd = "[imageQ,bufferQ] = matgadget.run_process(recon_data); matgadget.emptyQ();";
    send_matlab_command(cmd);


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
	mxDestroyArray(reconArray); //We're not supposed to delete this?

	// We are finished with the incoming messages m1 and m2
	m1->release();

	return GADGET_OK;
}
}


