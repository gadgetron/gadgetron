#include "GadgetIsmrmrdReadWrite.h"
#include "BucketToBufferGadget.h"
#include "Gadgetron.h"
#include "mri_core_data.h"

namespace Gadgetron{

BucketToBufferGadget::~BucketToBufferGadget()
{
    //The buckets array should be empty but just in case, let's make sure all the stuff is released.
}

int BucketToBufferGadget
::process_config(ACE_Message_Block* mb)
{

    std::string N_dimension = *this->get_string_value("N_dimension");
    std::string S_dimension = *this->get_string_value("S_dimension");
    
    if (N_dimension.size() == 0) {
        N_ = NONE;
    } else if (N_dimension.compare("average") == 0) {
        N_ = AVERAGE;
    } else if (N_dimension.compare("contrast") == 0) {
        N_ = CONTRAST;
    } else if (N_dimension.compare("phase") == 0) {
        N_ = PHASE;
    } else if (N_dimension.compare("repetition") == 0) {
        N_ = REPETITION;
    } else if (N_dimension.compare("set") == 0) {
        N_ = SET;
    } else if (N_dimension.compare("segment") == 0) {
        N_ = SEGMENT;
    } else {
        GADGET_DEBUG2("WARNING: Unknown N dimension (%s), N set to NONE", N_dimension.c_str());
        N_ = NONE;
    }
  
    GADGET_DEBUG2("N DIMENSION IS: %s (%d)\n", N_dimension.c_str(), N_);

    if (S_dimension.size() == 0) {
        S_ = NONE;
    } else if (S_dimension.compare("average") == 0) {
        S_ = AVERAGE;
    } else if (S_dimension.compare("contrast") == 0) {
        S_ = CONTRAST;
    } else if (S_dimension.compare("phase") == 0) {
        S_ = PHASE;
    } else if (S_dimension.compare("repetition") == 0) {
        S_ = REPETITION;
    } else if (S_dimension.compare("set") == 0) {
        S_ = SET;
    } else if (S_dimension.compare("segment") == 0) {
        S_ = SEGMENT;
    } else {
        GADGET_DEBUG2("WARNING: Unknown sort dimension (%s), sorting set to NONE\n", S_dimension.c_str());
        S_ = NONE;
    }
  
    GADGET_DEBUG2("S DIMENSION IS: %s (%d)\n", S_dimension.c_str(), S_);

    split_slices_  = this->get_bool_value("split_slices");
    GADGET_DEBUG2("SPLIT SLICES IS: %b\n", split_slices_);

    return GADGET_OK;
}

int BucketToBufferGadget
::process(GadgetContainerMessage<IsmrmrdAcquisitionBucket>* m1)
{

    size_t key;
    std::map<size_t, GadgetContainerMessage<IsmrmrdReconData>* > recon_data_buffers;
 
    //GADGET_DEBUG1("BucketToBufferGadget::process\n");
    //Iterate over the reference data of the bucket
    for(std::vector<IsmrmrdAcquisitionData>::iterator it = m1->getObjectPtr()->ref_.begin();
        it != m1->getObjectPtr()->ref_.end(); ++it)
    {
        // std::cout << it->head_->getObjectPtr()->idx.kspace_encode_step_1 << std::endl;
          
        //Generate the key to the corresponding ReconData buffer
        key = getKey(it->head_->getObjectPtr()->idx);
        //Look up the corresponding ReconData buffer
        if (recon_data_buffers.find(key) == recon_data_buffers.end()) {
            //ReconData buffer does not exist, create it
            recon_data_buffers[key] = new GadgetContainerMessage<IsmrmrdReconData>;
            // Set the size of the arrays
        }

        //Stuff it
        //recon_data_buffers[key]->getObjectPtr()->
    }
      
    ////Iterate over the imaging data of the bucket
    //for(std::vector<IsmrmrdAcquisitionData>::iterator it = m1->getObjectPtr()->data_.begin();
    //    it != m1->getObjectPtr()->data_.end(); ++it)
    //{
    //    // std::cout << it->head_->getObjectPtr()->idx.kspace_encode_step_1 << std::endl;
    //
    //    //Generate the key to the corresponding ReconData buffer
    //    key = getKey(it->head_->getObjectPtr()->idx);
    //    //Look up the corresponding ReconData buffer
    //    if (recon_data_buffers.find(key) == recon_data_buffers.end()) {
    //        //ReconData buffer does not exist, create it
    //        recon_data_buffers[key] = new GadgetContainerMessage<IsmrmrdReconData>;
    //        //Initialize the data array
    //        recon_data_buffers[key].data_.            
    //    }
    //    //Check and update the limits
    //    recon_data_buffers[key].data_[
    //}

    //Send all the ReconData messages
    GADGET_DEBUG2("End of bucket reached, sending out %d ReconData buffers\n", recon_data_buffers.size());
    for(std::map<size_t, GadgetContainerMessage<IsmrmrdReconData>* >::iterator it = recon_data_buffers.begin(); it != recon_data_buffers.end(); it++)
    {
        if (it->second) {
            if (this->next()->putq(it->second) == -1) {
                it->second->release();
                GADGET_DEBUG1("Failed to pass bucket down the chain\n");
                return GADGET_FAIL;
            }
        }
    }

    //Clear the recondata buffer map
    recon_data_buffers.clear();  // is this necessary?
      
    //We can release the incoming bucket now. This will release all of the data it contains.
    m1->release();

    return GADGET_OK;
}

int BucketToBufferGadget::close(unsigned long flags)
{
    
    int ret = Gadget::close(flags);
    GADGET_DEBUG1("BucketToBufferGadget::close\n");

    return ret;
}

size_t BucketToBufferGadget::getSlice(ISMRMRD::ISMRMRD_EncodingCounters idx)
{
    size_t index;
    
    if( split_slices_ ) {
        index = idx.slice;
    } else {
        index = 0;
    }
    
    return index;
}

size_t BucketToBufferGadget::getN(ISMRMRD::ISMRMRD_EncodingCounters idx)
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

size_t BucketToBufferGadget::getS(ISMRMRD::ISMRMRD_EncodingCounters idx)
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

size_t BucketToBufferGadget::getKey(ISMRMRD::ISMRMRD_EncodingCounters idx)
{
    //[RO, E1, E2, CHA, SLC, PHS, CON, REP, SET, SEG, AVE]
    //[SLC, PHS, CON, REP, SET, SEG, AVE]
    //        collapse acros two of them (N and S)
    
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

    if ((N_ == SEGMENT) || (S_ == SEGMENT)) {
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

GADGET_FACTORY_DECLARE(BucketToBufferGadget)

}
