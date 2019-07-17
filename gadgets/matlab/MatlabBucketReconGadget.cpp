#include "MatlabBucketReconGadget.h"
#include "../mri_core/GadgetIsmrmrdReadWrite.h" //LA: added ../mri_core/, is that the correct way ?
#include "mri_core_data.h"
#include "mri_core_def.h"
#include "hoNDArray_elemwise.h"
#include "hoNDArray_reductions.h"
#include "MatlabUtils.h"


std::mutex mutex_MBRG_;

namespace Gadgetron{

MatlabBucketReconGadget::MatlabBucketReconGadget()
{
}

MatlabBucketReconGadget::~MatlabBucketReconGadget()
{
    std::lock_guard<std::mutex> lock(mutex_MBRG_);   
    // Close the Matlab engine
    GDEBUG("Closing down Matlab\n");
    engClose(engine_);
}

int MatlabBucketReconGadget::process_config(ACE_Message_Block* mb)
{
    if (N_dimension.value().size() == 0) {
        N_ = NONE;
    } else if (N_dimension.value().compare("average") == 0) {
        N_ = AVERAGE;
    } else if (N_dimension.value().compare("contrast") == 0) {
        N_ = CONTRAST;
    } else if (N_dimension.value().compare("phase") == 0) {
        N_ = PHASE;
    } else if (N_dimension.value().compare("repetition") == 0) {
        N_ = REPETITION;
    } else if (N_dimension.value().compare("set") == 0) {
        N_ = SET;
    } else if (N_dimension.value().compare("segment") == 0) {
        N_ = SEGMENT;
    } else if (N_dimension.value().compare("slice") == 0){
        N_ = SLICE;
    } else {
        GDEBUG("WARNING: Unknown N dimension (%s), N set to NONE", N_dimension.value().c_str());
        N_ = NONE;
    }

    GDEBUG("N DIMENSION IS: %s (%d)\n", N_dimension.value().c_str(), N_);

    if (S_dimension.value().size() == 0) {
        S_ = NONE;
    } else if (S_dimension.value().compare("average") == 0) {
        S_ = AVERAGE;
    } else if (S_dimension.value().compare("contrast") == 0) {
        S_ = CONTRAST;
    } else if (S_dimension.value().compare("phase") == 0) {
        S_ = PHASE;
    } else if (S_dimension.value().compare("repetition") == 0) {
        S_ = REPETITION;
    } else if (S_dimension.value().compare("set") == 0) {
        S_ = SET;
    } else if (S_dimension.value().compare("segment") == 0) {
        S_ = SEGMENT;
    } else if (S_dimension.value().compare("slice") == 0){
        S_ = SLICE;
    } else {
        GDEBUG("WARNING: Unknown sort dimension (%s), sorting set to NONE\n", S_dimension.value().c_str());
        S_ = NONE;
    }

    GDEBUG("S DIMENSION IS: %s (%d)\n", S_dimension.value().c_str(), S_);

    split_slices_  = split_slices.value();
    GDEBUG("SPLIT SLICES IS: %b\n", split_slices_);

    ignore_segment_  = ignore_segment.value();
    GDEBUG("IGNORE SEGMENT IS: %b\n", ignore_segment_);

    // keep a copy of the deserialized ismrmrd xml header for runtime
    ISMRMRD::deserialize(mb->rd_ptr(), hdr_);
    
    
    std::lock_guard<std::mutex> lock(mutex_MBRG_);   
    std::string cmd;

    debug_mode_  = debug_mode.value();
    path_        = matlab_path.value();
    classname_   = matlab_classname.value();
    startcmd_    = matlab_startcmd.value();

    if (classname_.empty()) {
        GERROR("Missing Matlab Gadget classname in config!");
        return GADGET_FAIL;
    }

    GDEBUG("MATLAB Class Name : %s\n", classname_.c_str());


    // Open the Matlab Engine on the current host
    GDEBUG("Starting MATLAB engine with command: %s\n", startcmd_.c_str());
    if (!(engine_ = engOpen(startcmd_.c_str()))) {
        // TODO: error checking!
        GDEBUG("Can't start MATLAB engine\n");
    } else {
        // Prepare a buffer for collecting Matlab's output
        char matlab_buffer_[2049] = "\0";
        engOutputBuffer(engine_, matlab_buffer_, 2048);

        // Add the necessary paths to the matlab environment
        std::string gadgetron_matlab_path = get_gadgetron_home().string() + "/share/gadgetron/matlab";
        std::string add_path_cmd = std::string("addpath('") + gadgetron_matlab_path + std::string("');");
        // Gadgetron matlab scripts
        engEvalString(engine_, add_path_cmd.c_str());
        // ISMRMRD matlab library
        engEvalString(engine_, "addpath(fullfile(getenv('ISMRMRD_HOME'), '/share/ismrmrd/matlab'));");

        GDEBUG("%s", matlab_buffer_);
    }

    char matlab_buffer_[20481] = "\0";
    engOutputBuffer(engine_, matlab_buffer_, 20480);

    // add user specified path for this gadget
    if (!path_.empty()) {
        cmd = "addpath('" + path_ + "');";
        send_matlab_command(cmd);
    }

    // Put the XML Header into the matlab workspace
    std::string xmlConfig = std::string(mb->rd_ptr());
    mxArray *xmlstring = mxCreateString(xmlConfig.c_str());
    engPutVariable(engine_, "xmlstring", xmlstring);

    // Instantiate the Matlab gadget object from the user specified class
    // Call matlab gadget's init method with the XML Header
    // and the user defined config method
    cmd = "matgadget = " + classname_ + "();";
    cmd += "matgadget.init(xmlstring); matgadget.config();";
    if (send_matlab_command(cmd) != GADGET_OK) {
        GDEBUG("Failed to send matlab command.\n");
        return GADGET_FAIL;
    }

    mxDestroyArray(xmlstring);
    return GADGET_OK;
}

int MatlabBucketReconGadget::send_matlab_command(std::string& command)
{
    char matlab_buffer_[8193] = "\0";
    engOutputBuffer(engine_, matlab_buffer_, 8192);
    engEvalString(engine_, command.c_str());
    if (debug_mode_) {
        GDEBUG("%s\n", matlab_buffer_);
    }
    return GADGET_OK;
}
  
  
int MatlabBucketReconGadget::process(GadgetContainerMessage<IsmrmrdAcquisitionBucket>* m1)
{
    
       std::lock_guard<std::mutex> lock(mutex_MBRG_);
 
    high_resolution_clock::time_point time1 = high_resolution_clock::now();
    
    IsmrmrdAcquisitionBucket* bucket = m1->getObjectPtr();
    
    // Number of RO data existing in this bucket.
    long RO_counter = bucket->data_.size()+bucket->ref_.size();

    bool init = false;
    uint16_t NE0;
    uint16_t NE1;
    uint16_t NE2;
    uint16_t NN;
    uint16_t NS;
    uint16_t NCHA;
    std::complex<float>* raw_data;
    uint32_t* phase_coordinates; //cannot go under 32 bit data, as it would be very quickly reached
    
    size_t key;
    std::map<size_t, GadgetContainerMessage<IsmrmrdReconData>* > recon_data_buffers;    
    
    // boolean for packet labelling
    uint8_t isLastPacket=0;
    uint8_t isFirstPacket=1;

    // Iterate over the RO lines of the bucket, copy them into raw_data
//     IsmrmrdDataBuffered* pCurrDataBuffer = NULL;
    Gadgetron::GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* headersToMatlab = NULL;

    //GDEBUG("xxx1\n");

    for (std::vector<IsmrmrdAcquisitionData>::iterator it = bucket->data_.begin(); it != bucket->data_.end(); ++it)
    {
        //Get a reference to the header and data for this acquisition
        ISMRMRD::AcquisitionHeader       & acqhdr  = *it->head_->getObjectPtr();
        hoNDArray< std::complex<float> > & acqdata = *it->data_->getObjectPtr();
        
        if(!init)
        {
            NCHA = acqhdr.active_channels;
            NE0 = acqhdr.number_of_samples - acqhdr.discard_pre - acqhdr.discard_post;
            
            key = getKey(acqhdr.idx);
            uint16_t espace = acqhdr.encoding_space_ref;
            IsmrmrdReconBit & rbit = getRBit(recon_data_buffers, key, espace);
            IsmrmrdDataBuffered & dataBuffer = rbit.data_;
            ISMRMRD::Encoding & encoding = hdr_.encoding[espace];
            IsmrmrdAcquisitionBucketStats & stats = m1->getObjectPtr()->datastats_[espace];
            fillSamplingDescription(dataBuffer.sampling_, encoding, stats, acqhdr, false);
            uint16_t* dims = getEncodingDimensions(acqhdr, encoding, stats, false);
            
            NE1 = dims[1];
            NE2 = dims[2];
            NN  = dims[5];
            NS  = dims[6];
            delete[] dims;
            
            raw_data = (std::complex<float>*) malloc(RO_counter*NCHA*NE0*sizeof(std::complex<float>));
            phase_coordinates = (uint32_t*) mxCalloc(RO_counter*NCHA,    sizeof(uint32_t));
            
            RO_counter = 0; //recycle this counter so that it can be used to track the index within this loop
            init=true;
        }
        
        // Here we're accessing at the 24th bit of flags. 24th is the bit index of ACQ_LAST_IN_MEASUREMENT
        uint64_t flags = it->head_->getObjectPtr()->flags;
	if((flags & ( 1 << 24 )) >> 24)  
        {
	   isLastPacket = 1;
        }
	if(isFirstPacket)
	{
	   headersToMatlab = it->head_;
        }
        
        //Copy this RO line into raw_data
        for (uint16_t cha = 0; cha < NCHA; cha++)
        {
            memcpy(raw_data + RO_counter*NCHA*NE0 + cha*NE0, &acqdata(acqhdr.discard_pre, cha), sizeof(std::complex<float>)*NE0);
            
            phase_coordinates[RO_counter*NCHA + cha] = (uint32_t)    1 + /* Matlab indices start from 1 */
                                                                     it->head_->getObjectPtr()->idx.kspace_encode_step_1 +
                                                       NE1          *it->head_->getObjectPtr()->idx.kspace_encode_step_2 +
                                                       NE1*NE2      *cha +
                                                       NE1*NE2*NCHA *it->head_->getObjectPtr()->idx.contrast;
            //I thought this would take less time to extract, but it doesn't improve anything, so better keep the standard format.
//             phase_coordinates[RO_counter*NCHA + cha] = (uint32_t)   cha +
//                                                        NCHA        *it->head_->getObjectPtr()->idx.contrast +
//                                                        NCHA*NN     *it->head_->getObjectPtr()->idx.kspace_encode_step_2 +
//                                                        NCHA*NN*NE2 *it->head_->getObjectPtr()->idx.kspace_encode_step_1; 
        }

        ++RO_counter;
    }

    //GDEBUG("xxx2\n");
    
    for (std::vector<IsmrmrdAcquisitionData>::iterator it = bucket->ref_.begin(); it != bucket->ref_.end(); ++it)
    {
        //Get a reference to the header and data for this acquisition
        ISMRMRD::AcquisitionHeader       & acqhdr  = *it->head_->getObjectPtr();
        hoNDArray< std::complex<float> > & acqdata = *it->data_->getObjectPtr();
        
        if(!init)
        {
            NCHA = acqhdr.active_channels;
            NE0 = acqhdr.number_of_samples - acqhdr.discard_pre - acqhdr.discard_post;
            
            key = getKey(acqhdr.idx);
            uint16_t espace = acqhdr.encoding_space_ref;
            IsmrmrdReconBit & rbit = getRBit(recon_data_buffers, key, espace);
            if (!rbit.ref_)
                rbit.ref_ = IsmrmrdDataBuffered();
            IsmrmrdDataBuffered & dataBuffer = *rbit.ref_;
            ISMRMRD::Encoding & encoding = hdr_.encoding[espace];
            IsmrmrdAcquisitionBucketStats & stats = m1->getObjectPtr()->refstats_[espace];
            
            fillSamplingDescription(dataBuffer.sampling_, encoding, stats, acqhdr, true);
            uint16_t* dims = getEncodingDimensions(acqhdr, encoding, stats, true);

            NE1 = dims[1];
            NE2 = dims[2];
            NN  = dims[5];
            NS  = dims[6];
            delete[] dims;

            raw_data = (std::complex<float>*) malloc(RO_counter*NCHA*NE0*sizeof(std::complex<float>));
            phase_coordinates = (uint32_t*) mxCalloc(RO_counter*NCHA,    sizeof(uint32_t));

            RO_counter = 0; //recycle this counter so that it can be used to track the index within this loop
            init=true;
        }
        
        // Here we're accessing at the 24th bit of flags. 24th is the bit index of ACQ_LAST_IN_MEASUREMENT
        uint64_t flags = it->head_->getObjectPtr()->flags;
        if((flags & ( 1 << 24 )) >> 24)  // JAC: we need to find a workaround for (RO_counter == 0) to prevent setting isLastPacket=1 too soon for non-EPI data
        {
            isLastPacket = 1;
        }
        //headersToMatlab = it->head_;
	if(isFirstPacket)
	{
	   headersToMatlab = it->head_;
        }
        
        //Copy this RO line into raw_data
        for (uint16_t cha = 0; cha < NCHA; cha++)
        {
            memcpy(raw_data + RO_counter*NCHA*NE0 + cha*NE0, &acqdata(acqhdr.discard_pre, cha), sizeof(std::complex<float>)*NE0);
            
            phase_coordinates[RO_counter*NCHA + cha] = (uint32_t)    1 + /* Matlab indices start from 1 */
                                                                     it->head_->getObjectPtr()->idx.kspace_encode_step_1 +
                                                       NE1          *it->head_->getObjectPtr()->idx.kspace_encode_step_2 +
                                                       NE1*NE2      *cha +
                                                       NE1*NE2*NCHA *it->head_->getObjectPtr()->idx.contrast;
        }

        ++RO_counter;
    }
    
    high_resolution_clock::time_point time2 = high_resolution_clock::now();
    
    //GDEBUG("xxx3\n");
    
    ///////////////////////// RAWDATA TO MXSTRUCT //////////////////////////
//     const char * field_names[] = {"data","trajectory","headers","samplingdescription", "kspace_encode_step"};
    const char * field_names[] = {"data", "headers", "kspace_encode_step"};
	auto mxstruct = mxCreateStructMatrix(1,1,3,field_names); // always check the number of fields, otherwise segfault
	if (!mxstruct)
        throw std::runtime_error("Failed to allocate Matlab struct");

    // create the packet. A copy of the data is being done here
    size_t packet_n_elem = RO_counter * NE0 * NCHA;
    size_t packet_ndim = 2;
    mwSize* packet_dims = new mwSize[packet_ndim]; //no need to delete ?

    packet_dims[0] = NE0;
    packet_dims[1] = RO_counter*NCHA; // this is thus echo x line x partition

    //GDEBUG("xxx4\n");

    // recon data
#if MX_HAS_INTERLEAVED_COMPLEX
    GDEBUG("----------------- Interleaved Complex -----------------\n");
#else
    GDEBUG("--------------- Non-interleaved Complex ---------------\n");
#endif

    auto mxdata =  mxCreateNumericArray(packet_ndim, packet_dims, mxSINGLE_CLASS, mxCOMPLEX);

    float* real_data=(float *)mxGetData(mxdata);
#if MX_HAS_INTERLEAVED_COMPLEX
    float* imag_data=real_data+1;
#else
    float* imag_data=(float *)mxGetImagData(mxdata);
#endif

    //GDEBUG("xxx5\n");

    // Copy from C++ to matlab
    for(size_t i = 0; i<packet_n_elem; ++i)
    {
#if MX_HAS_INTERLEAVED_COMPLEX
        real_data[i*2] = real(raw_data[i]);
        imag_data[i*2] = imag(raw_data[i]);
#else
        real_data[i] = real(raw_data[i]);
        imag_data[i] = imag(raw_data[i]);
#endif
    }
    delete[] raw_data;

    mxSetField(mxstruct,0,"data",mxdata);
    
    //GDEBUG("xxx6\n");

    // encode
    auto mxstep = mxCreateNumericMatrix((mwSize) RO_counter*NCHA, 1, mxUINT32_CLASS, mxREAL);
    mxFree(mxGetData(mxstep));
    mxSetData(mxstep, phase_coordinates);
    mxSetField(mxstruct,0,"kspace_encode_step",mxstep);
    
    //GDEBUG("xxx7\n");

    /*
	//Add trajectory if available
	if (pCurrDataBuffer->trajectory_){
		auto & trajectory = *pCurrDataBuffer->trajectory_;
		int traj_fieldnumber = mxAddField(mxstruct,"trajectory");
		auto mxtraj = hoNDArrayToMatlab(&trajectory);
		mxSetFieldByNumber(mxstruct,0,traj_fieldnumber,mxtraj);
	}
    
    
	//Add headers
	std::cout << "Adding headers...";
	mwSize num_headers = pCurrDataBuffer->headers_.get_number_of_elements();
	auto mxheaders = mxCreateNumericMatrix(sizeof(ISMRMRD::AcquisitionHeader),num_headers,mxUINT8_CLASS,mxREAL);
	memcpy(mxGetData(mxheaders),pCurrDataBuffer->headers_.get_data_ptr(),sizeof(ISMRMRD::AcquisitionHeader)*num_headers);
	mxSetField(mxstruct,0,"headers",mxheaders);
    
    
	auto samplingdescription = samplingdescriptionToMatlabStruct(&pCurrDataBuffer->sampling_);
	mxSetField(mxstruct,0,"samplingdescription",samplingdescription);
    */
//     for(int i=0;i<pCurrDataBuffer->headers_.get_number_of_elements(); ++i)
//         std::cout << (*(pCurrDataBuffer->headers_.get_data_ptr())).patient_table_position[2] << std::endl;
//     
//     std::cout << "nelem: "<< headersToMatlab->getObjectPtr()->get_number_of_elements() << std::endl;
    
    if(isFirstPacket)
	{   isFirstPacket = 0;
        mwSize num_headers = 1; //headersToMatlab->size(); //= 1;
        auto mxheaders = mxCreateNumericMatrix(sizeof(ISMRMRD::AcquisitionHeader),num_headers,mxUINT8_CLASS,mxREAL);
        memcpy(mxGetData(mxheaders),headersToMatlab->getObjectPtr(),sizeof(ISMRMRD::AcquisitionHeader)*num_headers);
	mxSetField(mxstruct,0,"headers",mxheaders);
    }

    //GDEBUG("xxx8\n");

    // Create matlab boolean and set value to C++
    auto mxIsLastPacket = mxCreateNumericMatrix(1, 1, mxINT8_CLASS, mxREAL); // destroyed by recon_array destroy later
    *(uint8_t *)mxGetData(mxIsLastPacket)=isLastPacket;
        
    //GDEBUG("xxx9\n");

	// Initialize a string for matlab commands
	std::string cmd;

	mwSize nencoding_spaces = 1;
	const char* fieldnames[3] = {"data","reference", "isLastPacket"};
	auto reconArray = mxCreateStructArray(1,&nencoding_spaces,3,fieldnames);

    //GDEBUG("xxx10\n");

    ///////////////////////////////////

    mxSetField(reconArray,0,"data",mxstruct);
    mxSetField(reconArray,0,"isLastPacket",mxIsLastPacket);
    
    high_resolution_clock::time_point time3 = high_resolution_clock::now();
    
    engPutVariable(engine_, "recon_data", reconArray);
    
    high_resolution_clock::time_point time4 = high_resolution_clock::now();
    
    cmd = "[imageQ,bufferQ] = matgadget.run_process(recon_data); matgadget.emptyQ();";
    send_matlab_command(cmd);
    
    high_resolution_clock::time_point time5 = high_resolution_clock::now();

    ///////////////////////// FINITION //////////////////////////
	// Get the size of the gadget's queue
	mxArray *imageQ = engGetVariable(engine_, "imageQ");
	if (imageQ == NULL) {
		GERROR("Failed to get the imageQ from matgadget\n");
		return GADGET_FAIL;
	}

	size_t qlen = mxGetNumberOfElements(imageQ);

	const mwSize* dims = mxGetDimensions(imageQ);
	mwSize ndims = mxGetNumberOfDimensions(imageQ);

	

	//Read all Image bytes
	for (mwIndex idx = 0; idx < qlen; idx++) {
		mxArray *res_hdr  = mxGetField(imageQ, idx, "bytes");
		mxArray *res_data = mxGetField(imageQ, idx, "image");
        mxArray *res_comment = mxGetField(imageQ, idx, "image_comment"); //LA: no idea when this is freed

		GadgetContainerMessage<ISMRMRD::ImageHeader>* m3 = new GadgetContainerMessage<ISMRMRD::ImageHeader>();
		ISMRMRD::ImageHeader *hdr_new = m3->getObjectPtr();
		memcpy(hdr_new, mxGetData(res_hdr), sizeof(ISMRMRD::ImageHeader));

		auto image= MatlabToHoNDArray<std::complex<float>>(res_data);
		auto m4 = new GadgetContainerMessage< hoNDArray< std::complex<float> > >(image);
		auto dims = *image.get_dimensions();

		m3->cont(m4);
        
        char* image_comment = mxArrayToString(res_comment);
        
        m3->getObjectPtr()->image_type = ISMRMRD::ISMRMRD_IMTYPE_MAGNITUDE;
        
        GadgetContainerMessage< ISMRMRD::MetaContainer >* cm3 = new GadgetContainerMessage< ISMRMRD::MetaContainer >();

        cm3->getObjectPtr()[0].append(GADGETRON_IMAGECOMMENT, image_comment);

        m4->cont(cm3);
        
		if (this->next()->putq(m3) < 0) {
			GDEBUG("Failed to put Image message on queue\n");
			return GADGET_FAIL;
		}

	}
	//Match engGetVariable with mxDestroy___s
	mxArray* bufferQ = engGetVariable(engine_,"bufferQ");
    
	qlen = mxGetNumberOfElements(bufferQ);
// 	if (debug_mode_) {
// 		GDEBUG("Buffer Queue size: %d \n", qlen);
//         GDEBUG("Image Queue size: %d \n", qlen);
//         GDEBUG("Number of ndims %i \n"  ,ndims);
//     }

	for (mwIndex idx = 0; idx <qlen; idx++){

		IsmrmrdReconData output_data;
		IsmrmrdReconBit bit;
		bit.data_ = MatlabStructToBuffer(mxGetField(bufferQ,idx,"data"));

        
		auto ref = mxGetField(bufferQ,idx,"reference");
		if (ref){
			GDEBUG("Adding reference");
			bit.ref_ = MatlabStructToBuffer(ref);
			
		}
		/*else {
		isLastPacket = 0;
		}*/

		output_data.rbit_.push_back(bit);
		auto m3 = new GadgetContainerMessage<IsmrmrdReconData>(output_data);
		if (this->next()->putq(m3) < 0){
			GDEBUG("Failed to put Buffer message on queue\n");
			return GADGET_FAIL;
		}

	}

	mxDestroyArray(bufferQ);
	mxDestroyArray(imageQ);
	mxDestroyArray(reconArray);

	m1->release();
    
    high_resolution_clock::time_point time6 = high_resolution_clock::now();

    
//     clock_t CPS = CLOCKS_PER_SEC;
    duration<double> time_span1 = duration_cast<duration<double>>(time1 - exitTime_);
    duration<double> time_span2 = duration_cast<duration<double>>(time2 - time1);
    duration<double> time_span3 = duration_cast<duration<double>>(time3 - time2);
    duration<double> time_span4 = duration_cast<duration<double>>(time4 - time3);
    duration<double> time_span5 = duration_cast<duration<double>>(time5 - time4);
    duration<double> time_span6 = duration_cast<duration<double>>(time6 - time5);
    
    //     exitTime_ = std::clock();
    exitTime_ = high_resolution_clock::now();
    
    std::cout   <<   "----------------- Execution times [s] -----------------"
                << "\nOutside process   : " << time_span1.count()
                << "\nBucket to raw data: " << time_span2.count()
                << "\nRaw data to mxdata: " << time_span3.count()
                << "\nmxdata transfer   : " << time_span4.count()
                << "\nMATLAB extraction : " << time_span5.count()
                << "\nFinition          : " << time_span6.count()
                << "\n-------------------------------------------------------\n";

    return GADGET_OK;
  }
  
  void MatlabBucketReconGadget::stuff(std::vector<IsmrmrdAcquisitionData>::iterator it, IsmrmrdDataBuffered & dataBuffer, ISMRMRD::Encoding encoding, IsmrmrdAcquisitionBucketStats & stats, bool forref)
  {

    // The acquisition header and data
    ISMRMRD::AcquisitionHeader & acqhdr = *it->head_->getObjectPtr();
    hoNDArray< std::complex<float> > & acqdata = *it->data_->getObjectPtr();
    // we make one for the trajectory down below if we need it

    uint16_t NE0  = (uint16_t)dataBuffer.data_.get_size(0);
    uint16_t NE1  = (uint16_t)dataBuffer.data_.get_size(1);
    uint16_t NE2  = (uint16_t)dataBuffer.data_.get_size(2);
    uint16_t NCHA = (uint16_t)dataBuffer.data_.get_size(3);
    uint16_t NN   = (uint16_t)dataBuffer.data_.get_size(4);
    uint16_t NS   = (uint16_t)dataBuffer.data_.get_size(5);
    uint16_t NLOC = (uint16_t)dataBuffer.data_.get_size(6);

    size_t slice_loc;
    if (split_slices_ || NLOC==1) {
        slice_loc = 0;
    }
    else {
        slice_loc = acqhdr.idx.slice;
    }

    //Stuff the data
    uint16_t npts_to_copy = acqhdr.number_of_samples - acqhdr.discard_pre - acqhdr.discard_post;
    long long offset;
    if ( ((encoding.trajectory == ISMRMRD::TrajectoryType::CARTESIAN)) || (encoding.trajectory == ISMRMRD::TrajectoryType::EPI) ) {
        if ((acqhdr.number_of_samples == dataBuffer.data_.get_size(0)) && (acqhdr.center_sample == acqhdr.number_of_samples/2)) // acq has been corrected for center , e.g. by asymmetric handling
        {
            offset = acqhdr.discard_pre;
        }
        else
        {
            offset = (long long)dataBuffer.sampling_.sampling_limits_[0].center_ - (long long)acqhdr.center_sample;
        }
    } else {
        //TODO what about EPI with asymmetric readouts?
        //TODO any other sort of trajectory?
        offset = 0;
    }
    long long roffset = (long long) dataBuffer.data_.get_size(0) - npts_to_copy - offset;


    if ((offset < 0) | (roffset < 0) )
      {
        throw std::runtime_error("Acquired reference data does not fit into the reference data buffer.\n");
      }

    std::complex<float> *dataptr;

    uint16_t NUsed = (uint16_t)getN(acqhdr.idx);
    if (NUsed >= NN) NUsed = NN - 1;

    uint16_t SUsed = (uint16_t)getS(acqhdr.idx);
    if (SUsed >= NS) SUsed = NS - 1;

    int16_t e1 = (int16_t)acqhdr.idx.kspace_encode_step_1;
    int16_t e2 = (int16_t)acqhdr.idx.kspace_encode_step_2;

    bool is_cartesian_sampling = ((encoding.trajectory == ISMRMRD::TrajectoryType::CARTESIAN));
    bool is_epi_sampling = (encoding.trajectory == ISMRMRD::TrajectoryType::EPI);
    if(is_cartesian_sampling || is_epi_sampling)
    {
        if (!forref || (forref && (encoding.parallelImaging.get().calibrationMode.get() == "embedded")))
        {
            // compute the center offset for E1 and E2
            int16_t space_matrix_offset_E1 = 0;
            if (encoding.encodingLimits.kspace_encoding_step_1.is_present())
            {
                space_matrix_offset_E1 = (int16_t)encoding.encodedSpace.matrixSize.y / 2 - (int16_t)encoding.encodingLimits.kspace_encoding_step_1->center;
            }

            int16_t space_matrix_offset_E2 = 0;
            if (encoding.encodingLimits.kspace_encoding_step_2.is_present() && encoding.encodedSpace.matrixSize.z > 1)
            {
                space_matrix_offset_E2 = (int16_t)encoding.encodedSpace.matrixSize.z / 2 - (int16_t)encoding.encodingLimits.kspace_encoding_step_2->center;
            }

            // compute the used e1 and e2 indices and make sure they are in the valid range
            e1 = (int16_t)acqhdr.idx.kspace_encode_step_1 + space_matrix_offset_E1;
            e2 = (int16_t)acqhdr.idx.kspace_encode_step_2 + space_matrix_offset_E2;
        }

        // for external or separate mode, it is possible the starting numbers of ref lines are not zero, therefore it is needed to subtract the staring ref line number
        // because the ref array size is set up by the actual number of lines acquired
        // only assumption for external or separate ref line mode is that all ref lines are numbered sequentially
        // the acquisition order of ref line can be arbitrary
        if (forref && ( (encoding.parallelImaging.get().calibrationMode.get() == "separate") || (encoding.parallelImaging.get().calibrationMode.get() == "external") ) )
        {
            if(*stats.kspace_encode_step_1.begin()>0)
            {
                e1 = acqhdr.idx.kspace_encode_step_1 - *stats.kspace_encode_step_1.begin();
            }

            if(*stats.kspace_encode_step_2.begin()>0)
            {
                e2 = acqhdr.idx.kspace_encode_step_2 - *stats.kspace_encode_step_2.begin();
            }
        }

        if (e1 < 0 || e1 >= (int16_t)NE1)
        {
            // if the incoming line is outside the encoding limits, something is wrong
            GADGET_CHECK_THROW(acqhdr.idx.kspace_encode_step_1>=encoding.encodingLimits.kspace_encoding_step_1->minimum && acqhdr.idx.kspace_encode_step_1 <= encoding.encodingLimits.kspace_encoding_step_1->maximum);

            // if the incoming line is inside encoding limits but outside the encoded matrix, do not include the data
            GWARN_STREAM("incoming readout " << acqhdr.scan_counter << " is inside the encoding limits, but outside the encoded matrix for kspace_encode_step_1 : " << e1 << " out of " << NE1);
            return;
        }

        if (e2 < 0 || e2 >= (int16_t)NE2)
        {
            GADGET_CHECK_THROW(acqhdr.idx.kspace_encode_step_2 >= encoding.encodingLimits.kspace_encoding_step_2->minimum && acqhdr.idx.kspace_encode_step_2 <= encoding.encodingLimits.kspace_encoding_step_2->maximum);

            GWARN_STREAM("incoming readout " << acqhdr.scan_counter << " is inside the encoding limits, but outside the encoded matrix for kspace_encode_step_2 : " << e2 << " out of " << NE2);
            return;
        }
    }

    std::complex<float>* pData = &dataBuffer.data_(offset, e1, e2, 0, NUsed, SUsed, slice_loc);

    // DATA COPY
    for (uint16_t cha = 0; cha < NCHA; cha++)
    {
        dataptr = pData + cha*NE0*NE1*NE2;
        memcpy(dataptr, &acqdata(acqhdr.discard_pre, cha), sizeof(std::complex<float>)*npts_to_copy);
    }

    dataBuffer.headers_(e1, e2, NUsed, SUsed, slice_loc) = acqhdr;

    if (acqhdr.trajectory_dimensions > 0)
    {

        hoNDArray< float > & acqtraj = *it->traj_->getObjectPtr();  // TODO do we need to check this?

        float * trajptr;

        trajptr = &(*dataBuffer.trajectory_)(0, offset, e1, e2, NUsed, SUsed, slice_loc);

        memcpy(trajptr, &acqtraj(0, acqhdr.discard_pre), sizeof(float)*npts_to_copy*acqhdr.trajectory_dimensions);

    }
  }
  
  

  int MatlabBucketReconGadget::close(unsigned long flags)
  {

    int ret = Gadget::close(flags);
    GDEBUG("MatlabBucketReconGadget::close\n");

    return ret;
  }

  size_t MatlabBucketReconGadget::getSlice(ISMRMRD::ISMRMRD_EncodingCounters idx)
  {
    size_t index;

    if( split_slices_ ) {
        index = idx.slice;
    } else {
        index = 0;
    }

    return index;
  }

  size_t MatlabBucketReconGadget::getN(ISMRMRD::ISMRMRD_EncodingCounters idx)
  {
    size_t index;

    if (N_ == AVERAGE) {
        index = idx.average;
    } else if (N_ == CONTRAST) {
        index = idx.contrast;
    } else if (N_ == PHASE) {
        index = idx.phase;
    } else if (N_ == REPETITION) {
        index = idx.repetition;
    } else if (N_ == SET) {
        index = idx.set;
    } else if (N_ == SEGMENT) {
        index = idx.segment;
    } else {
        index = 0;
    }

    return index;
  }

  size_t MatlabBucketReconGadget::getS(ISMRMRD::ISMRMRD_EncodingCounters idx)
  {
    size_t index;

    if (S_ == AVERAGE) {
        index = idx.average;
    } else if (S_ == CONTRAST) {
        index = idx.contrast;
    } else if (S_ == PHASE) {
        index = idx.phase;
    } else if (S_ == REPETITION) {
        index = idx.repetition;
    } else if (S_ == SET) {
        index = idx.set;
    } else if (S_ == SEGMENT) {
        index = idx.segment;
    } else {
        index = 0;
    }

    return index;
  }

  size_t MatlabBucketReconGadget::getKey(ISMRMRD::ISMRMRD_EncodingCounters idx)
  {
    //[SLC, PHS, CON, REP, SET, SEG, AVE]
    //collapse across two of them (N and S)

    size_t slice, phase, contrast, repetition, set, segment, average;

    if (split_slices_) {
        slice = idx.slice;
    } else {
        slice = 0;
    }

    if ((N_ == PHASE) || (S_ == PHASE)) {
        phase = 0;
    } else {
        phase = idx.phase;
    }

    if ((N_ == CONTRAST) || (S_ == CONTRAST)) {
        contrast = 0;
    } else {
        contrast = idx.contrast;
    }

    if ((N_ == REPETITION) || (S_ == REPETITION)) {
        repetition = 0;
    } else {
        repetition = idx.repetition;
    }

    if ((N_ == SET) || (S_ == SET)) {
        set = 0;
    } else {
        set = idx.set;
    }

    if ((N_ == SEGMENT) || (S_ == SEGMENT) || ignore_segment_) {
        segment = 0;
    } else {
        segment = idx.segment;
    }

    if ((S_ == AVERAGE) || (N_ == AVERAGE)) {
        average = 0;
    } else {
        average = idx.average;
    }

    size_t key = 0;
    key += slice      * 0x1;
    key += phase      * 0x100;
    key += contrast   * 0x10000;
    key += repetition * 0x1000000;
    key += set        * 0x100000000;
    key += segment    * 0x10000000000;
    key += average    * 0x1000000000000;

    return key;
  }

  IsmrmrdReconBit & MatlabBucketReconGadget::getRBit(std::map<size_t, GadgetContainerMessage<IsmrmrdReconData>* > & recon_data_buffers, size_t key, uint16_t espace)
  {
    //Look up the corresponding ReconData buffer
    if (recon_data_buffers.find(key) == recon_data_buffers.end())
      {
        //ReconData buffer does not exist, create it
        recon_data_buffers[key] = new GadgetContainerMessage<IsmrmrdReconData>;
      }

    //Look up the DataBuffered entry corresponding to this encoding space
    // create if needed and set the fields of view and matrix size
    if ( recon_data_buffers[key]->getObjectPtr()->rbit_.size() < (espace+1) )
      {
        recon_data_buffers[key]->getObjectPtr()->rbit_.resize(espace+1);
      }

    return recon_data_buffers[key]->getObjectPtr()->rbit_[espace];

  }

  
uint16_t* MatlabBucketReconGadget::getEncodingDimensions(ISMRMRD::AcquisitionHeader & acqhdr, ISMRMRD::Encoding encoding, IsmrmrdAcquisitionBucketStats & stats, bool forref)
{
    //// NE0 ////
    uint16_t NE0;
    if ( ((encoding.trajectory == ISMRMRD::TrajectoryType::CARTESIAN)) || (encoding.trajectory == ISMRMRD::TrajectoryType::EPI) ) {
        // if seperate or external calibration mode, using the acq length for NE0
        if (encoding.parallelImaging)
        {
            NE0 = acqhdr.number_of_samples;
        }
        else
        {
            NE0 = acqhdr.number_of_samples - acqhdr.discard_pre - acqhdr.discard_post;
        }
    } else {
        NE0 = acqhdr.number_of_samples - acqhdr.discard_pre - acqhdr.discard_post;
    }
    
    //// NE1 ////
    uint16_t NE1;
    if ( ((encoding.trajectory == ISMRMRD::TrajectoryType::CARTESIAN)) || (encoding.trajectory == ISMRMRD::TrajectoryType::EPI) )
    {
        if (encoding.parallelImaging)
        {
            if (forref && (encoding.parallelImaging.get().calibrationMode.get() == "separate"
                || encoding.parallelImaging.get().calibrationMode.get() == "external"))
            {
                NE1 = *stats.kspace_encode_step_1.rbegin() - *stats.kspace_encode_step_1.begin() + 1;
            }
            else
            {
                NE1 = encoding.encodedSpace.matrixSize.y;
            }
        }
        else
        {
            if (encoding.encodingLimits.kspace_encoding_step_1.is_present())
            {
                NE1 = encoding.encodingLimits.kspace_encoding_step_1->maximum - encoding.encodingLimits.kspace_encoding_step_1->minimum + 1;
            }
            else
            {
                NE1 = encoding.encodedSpace.matrixSize.y;
            }
        }
    }
    else {
        if (encoding.encodingLimits.kspace_encoding_step_1.is_present()) {
            NE1 = encoding.encodingLimits.kspace_encoding_step_1->maximum - encoding.encodingLimits.kspace_encoding_step_1->minimum + 1;
        }
        else {
            NE1 = *stats.kspace_encode_step_1.rbegin() - *stats.kspace_encode_step_1.begin() + 1;
        }
    }
    
    //// NE2 ////
    uint16_t NE2;
    if ( ((encoding.trajectory == ISMRMRD::TrajectoryType::CARTESIAN)) || (encoding.trajectory == ISMRMRD::TrajectoryType::EPI) )
    {
        if (encoding.parallelImaging)
        {
            if (forref && (encoding.parallelImaging.get().calibrationMode.get() == "separate" || encoding.parallelImaging.get().calibrationMode.get() == "external"))
            {
                NE2 = encoding.encodingLimits.kspace_encoding_step_2->maximum - encoding.encodingLimits.kspace_encoding_step_2->minimum + 1;
            }
            else
            {
                NE2 = encoding.encodedSpace.matrixSize.z;
            }
        }
        else
        {
            if (encoding.encodingLimits.kspace_encoding_step_2.is_present())
            {
                NE2 = encoding.encodingLimits.kspace_encoding_step_2->maximum - encoding.encodingLimits.kspace_encoding_step_2->minimum + 1;
            }
            else
            {
                NE2 = encoding.encodedSpace.matrixSize.z;
            }
        }
    }
    else {
        if (encoding.encodingLimits.kspace_encoding_step_2.is_present())
        {
            NE2 = encoding.encodingLimits.kspace_encoding_step_2->maximum - encoding.encodingLimits.kspace_encoding_step_2->minimum + 1;
        }
        else
        {
            NE2 = *stats.kspace_encode_step_2.rbegin() - *stats.kspace_encode_step_2.begin() + 1;
        }
    }
    
    //// NCHA ////
    uint16_t NCHA = acqhdr.active_channels;
    
    //// NLOC ////
    uint16_t NLOC;
    if (split_slices_)
    {
        NLOC = 1;
    }
    else
    {
        
        if (encoding.encodingLimits.slice.is_present())
        {
            NLOC = encoding.encodingLimits.slice->maximum - encoding.encodingLimits.slice->minimum + 1;
        }
        else
        {
            NLOC = 1;
        }
        // if the AcquisitionAccumulateTriggerGadget sort by SLC, then the stats should be used to determine NLOC
        size_t NLOC_received = *stats.slice.rbegin() - *stats.slice.begin() + 1; // crash here if forref if true

        if (NLOC_received < NLOC)
        {
            NLOC = NLOC_received;
        }
    }
    //// NN ////
    uint16_t NN;
    switch (N_) {
    case PHASE:
      NN = *stats.phase.rbegin() - *stats.phase.begin() + 1;
      break;
    case CONTRAST:
      NN = *stats.contrast.rbegin() - *stats.contrast.begin() + 1;
      break;
    case REPETITION:
      NN = *stats.repetition.rbegin() - *stats.repetition.begin() + 1;
      break;
    case SET:
      NN = *stats.set.rbegin() - *stats.set.begin() + 1;
      break;
    case SEGMENT:
      NN = *stats.segment.rbegin() - *stats.segment.begin() + 1;
      break;
    case AVERAGE:
      NN = *stats.average.rbegin() - *stats.average.begin() + 1;
      break;
    case SLICE:
      NN =  *stats.slice.rbegin() - *stats.slice.begin() + 1;
      break;
    default:
      NN = 1;
    }
    //// NS ////
    uint16_t NS;
    switch (S_) {
    case PHASE:
      NS = *stats.phase.rbegin() - *stats.phase.begin() + 1;
      break;
    case CONTRAST:
      NS = *stats.contrast.rbegin() - *stats.contrast.begin() + 1;
      break;
    case REPETITION:
      NS = *stats.repetition.rbegin() - *stats.repetition.begin() + 1;
      break;
    case SET:
      NS = *stats.set.rbegin() - *stats.set.begin() + 1;
      break;
    case SEGMENT:
      NS = *stats.segment.rbegin() - *stats.segment.begin() + 1;
      break;
    case AVERAGE:
      NS = *stats.average.rbegin() - *stats.average.begin() + 1;
      break;
    case SLICE:
      NS =  *stats.slice.rbegin() - *stats.slice.begin() + 1;
      break;
    default:
      NS = 1;
    }
    
    return new uint16_t[7]{NE0,NE1,NE2,NCHA,NLOC,NN,NS};
}
  

  void MatlabBucketReconGadget::fillSamplingDescription(SamplingDescription & sampling, ISMRMRD::Encoding & encoding, IsmrmrdAcquisitionBucketStats & stats, ISMRMRD::AcquisitionHeader& acqhdr, bool forref)
  {
    // For cartesian trajectories, assume that any oversampling has been removed.
    if ( ((encoding.trajectory == ISMRMRD::TrajectoryType::CARTESIAN)) ) {
        sampling.encoded_FOV_[0] = encoding.reconSpace.fieldOfView_mm.x;
        sampling.encoded_matrix_[0] = encoding.reconSpace.matrixSize.x;
    } else {
        sampling.encoded_FOV_[0] = encoding.encodedSpace.fieldOfView_mm.x;
        sampling.encoded_matrix_[0] = encoding.encodedSpace.matrixSize.x;
    }

    sampling.encoded_FOV_[1] = encoding.encodedSpace.fieldOfView_mm.y;
    sampling.encoded_FOV_[2] = encoding.encodedSpace.fieldOfView_mm.z;

    sampling.encoded_matrix_[1] = encoding.encodedSpace.matrixSize.y;
    sampling.encoded_matrix_[2] = encoding.encodedSpace.matrixSize.z;

    sampling.recon_FOV_[0] = encoding.reconSpace.fieldOfView_mm.x;
    sampling.recon_FOV_[1] = encoding.reconSpace.fieldOfView_mm.y;
    sampling.recon_FOV_[2] = encoding.reconSpace.fieldOfView_mm.z;

    sampling.recon_matrix_[0] = encoding.reconSpace.matrixSize.x;
    sampling.recon_matrix_[1] = encoding.reconSpace.matrixSize.y;
    sampling.recon_matrix_[2] = encoding.reconSpace.matrixSize.z;

    // For cartesian trajectories, assume that any oversampling has been removed.
    if ( ((encoding.trajectory == ISMRMRD::TrajectoryType::CARTESIAN)) || (encoding.trajectory == ISMRMRD::TrajectoryType::EPI) ) {
        sampling.sampling_limits_[0].min_ = acqhdr.discard_pre;
        sampling.sampling_limits_[0].max_ = acqhdr.number_of_samples - acqhdr.discard_post - 1;
        sampling.sampling_limits_[0].center_ = acqhdr.number_of_samples / 2;
    } else {
        sampling.sampling_limits_[0].min_ = 0;
        sampling.sampling_limits_[0].max_ = encoding.encodedSpace.matrixSize.x - 1;
        sampling.sampling_limits_[0].center_ = encoding.encodedSpace.matrixSize.x / 2;
    }

    // if the scan is cartesian  
        if ( ( (encoding.trajectory == ISMRMRD::TrajectoryType::CARTESIAN) && (!forref || (forref && (encoding.parallelImaging.get().calibrationMode.get() == "embedded"))) )
        || ( (encoding.trajectory == ISMRMRD::TrajectoryType::EPI) && !forref) )
    {
        int16_t space_matrix_offset_E1 = 0;
        if (encoding.encodingLimits.kspace_encoding_step_1.is_present())
        {
            space_matrix_offset_E1 = (int16_t)encoding.encodedSpace.matrixSize.y / 2 - (int16_t)encoding.encodingLimits.kspace_encoding_step_1->center;
        }

        int16_t space_matrix_offset_E2 = 0;
        if (encoding.encodingLimits.kspace_encoding_step_2.is_present() && encoding.encodedSpace.matrixSize.z > 1)
        {
            space_matrix_offset_E2 = (int16_t)encoding.encodedSpace.matrixSize.z / 2 - (int16_t)encoding.encodingLimits.kspace_encoding_step_2->center;
        }

        // E1
        sampling.sampling_limits_[1].min_ = encoding.encodingLimits.kspace_encoding_step_1->minimum + space_matrix_offset_E1;
        sampling.sampling_limits_[1].max_ = encoding.encodingLimits.kspace_encoding_step_1->maximum + space_matrix_offset_E1;
        sampling.sampling_limits_[1].center_ = sampling.encoded_matrix_[1] / 2;

        GADGET_CHECK_THROW(sampling.sampling_limits_[1].min_ < encoding.encodedSpace.matrixSize.y);
        GADGET_CHECK_THROW(sampling.sampling_limits_[1].max_ >= sampling.sampling_limits_[1].min_);
        GADGET_CHECK_THROW(sampling.sampling_limits_[1].center_ >= sampling.sampling_limits_[1].min_);
        GADGET_CHECK_THROW(sampling.sampling_limits_[1].center_ <= sampling.sampling_limits_[1].max_);

        // E2
        sampling.sampling_limits_[2].min_ = encoding.encodingLimits.kspace_encoding_step_2->minimum + space_matrix_offset_E2;
        sampling.sampling_limits_[2].max_ = encoding.encodingLimits.kspace_encoding_step_2->maximum + space_matrix_offset_E2;
        sampling.sampling_limits_[2].center_ = sampling.encoded_matrix_[2] / 2;

        GADGET_CHECK_THROW(sampling.sampling_limits_[2].min_ < encoding.encodedSpace.matrixSize.y);
        GADGET_CHECK_THROW(sampling.sampling_limits_[2].max_ >= sampling.sampling_limits_[2].min_);
        GADGET_CHECK_THROW(sampling.sampling_limits_[2].center_ >= sampling.sampling_limits_[2].min_);
        GADGET_CHECK_THROW(sampling.sampling_limits_[2].center_ <= sampling.sampling_limits_[2].max_);
    }
    else
    {
        sampling.sampling_limits_[1].min_ = encoding.encodingLimits.kspace_encoding_step_1->minimum;
        sampling.sampling_limits_[1].max_ = encoding.encodingLimits.kspace_encoding_step_1->maximum;
        sampling.sampling_limits_[1].center_ = encoding.encodingLimits.kspace_encoding_step_1->center;

        sampling.sampling_limits_[2].min_ = encoding.encodingLimits.kspace_encoding_step_2->minimum;
        sampling.sampling_limits_[2].max_ = encoding.encodingLimits.kspace_encoding_step_2->maximum;
        sampling.sampling_limits_[2].center_ = encoding.encodingLimits.kspace_encoding_step_2->center;
    }

    if (verbose.value())
    {
        GDEBUG_STREAM("Encoding space : " << int(encoding.trajectory)
            << " - FOV : [ " << encoding.encodedSpace.fieldOfView_mm.x << " " << encoding.encodedSpace.fieldOfView_mm.y << " " << encoding.encodedSpace.fieldOfView_mm.z << " ] "
            << " - Matris size : [ " << encoding.encodedSpace.matrixSize.x << " " << encoding.encodedSpace.matrixSize.y << " " << encoding.encodedSpace.matrixSize.z << " ] ");

        GDEBUG_STREAM("Sampling limits : "
                << "- RO : [ " << sampling.sampling_limits_[0].min_ << " " << sampling.sampling_limits_[0].center_ << " " << sampling.sampling_limits_[0].max_
                << " ] - E1 : [ " << sampling.sampling_limits_[1].min_ << " " << sampling.sampling_limits_[1].center_ << " " << sampling.sampling_limits_[1].max_
                << " ] - E2 : [ " << sampling.sampling_limits_[2].min_ << " " << sampling.sampling_limits_[2].center_ << " " << sampling.sampling_limits_[2].max_ << " ]");
    }
  }

  

  GADGET_FACTORY_DECLARE(MatlabBucketReconGadget)

}

