#include "GadgetIsmrmrdReadWrite.h"
#include "BucketToBufferGadget.h"
#include "mri_core_data.h"
#include "hoNDArray_elemwise.h"
#include "hoNDArray_reductions.h"
namespace Gadgetron{

  BucketToBufferGadget::~BucketToBufferGadget()
  {
    //The buckets array should be empty but just in case, let's make sure all the stuff is released.
  }

  int BucketToBufferGadget
  ::process_config(ACE_Message_Block* mb)
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
    } else if (N_dimension.value().compare("slice") == 0){
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

    return GADGET_OK;
  }

  int BucketToBufferGadget
  ::process(GadgetContainerMessage<IsmrmrdAcquisitionBucket>* m1)
  {

    size_t key;
    std::map<size_t, GadgetContainerMessage<IsmrmrdReconData>* > recon_data_buffers;

    //GDEBUG("BucketToBufferGadget::process\n");

    //Some information about the bucket
    //GDEBUG_STREAM("The Reference part: " << m1->getObjectPtr()->refstats_.size() << std::endl);
    //GDEBUG_STREAM("   nslices: " << m1->getObjectPtr()->refstats_[0].slice.size() << std::endl);
    //for (int e=0; e<m1->getObjectPtr()->refstats_.size() ; e++) {
    //    for (std::set<uint16_t>::iterator it = m1->getObjectPtr()->refstats_[e].kspace_encode_step_1.begin();
    //         it != m1->getObjectPtr()->refstats_[e].kspace_encode_step_1.end(); ++it) {
    //        GDEBUG_STREAM("   K1: " <<  *it << std::endl);
    //    }
    //}
    //GDEBUG_STREAM("The data part: " << m1->getObjectPtr()->datastats_.size() << std::endl);
    //GDEBUG_STREAM("   nslices: " << m1->getObjectPtr()->datastats_[0].slice.size() << std::endl);
    //for (int e=0; e<m1->getObjectPtr()->datastats_.size() ; e++) {
    //    for (std::set<uint16_t>::iterator it = m1->getObjectPtr()->datastats_[e].kspace_encode_step_1.begin();
    //         it != m1->getObjectPtr()->datastats_[e].kspace_encode_step_1.end(); ++it) {
    //        GDEBUG_STREAM("   K1: " <<  *it << std::endl);
    //    }
    //}

    //Iterate over the reference data of the bucket
    size_t count = 0;
    for(std::vector<IsmrmrdAcquisitionData>::iterator it = m1->getObjectPtr()->ref_.begin();
        it != m1->getObjectPtr()->ref_.end(); ++it)
      {
        //Get a reference to the header for this acquisition
        ISMRMRD::AcquisitionHeader & acqhdr = *it->head_->getObjectPtr();

        //Generate the key to the corresponding ReconData buffer
        key = getKey(acqhdr.idx);

        //The storage is based on the encoding space
        uint16_t espace = acqhdr.encoding_space_ref;

        //GDEBUG_STREAM("espace: " << acqhdr.encoding_space_ref << std::endl);
        //GDEBUG_STREAM("slice: " << acqhdr.idx.slice << std::endl);
        //GDEBUG_STREAM("rep: " << acqhdr.idx.repetition << std::endl);
        //GDEBUG_STREAM("k1: " << acqhdr.idx.kspace_encode_step_1 << std::endl);
        //GDEBUG_STREAM("k2: " << acqhdr.idx.kspace_encode_step_2 << std::endl);
        //GDEBUG_STREAM("seg: " << acqhdr.idx.segment << std::endl);
        //GDEBUG_STREAM("key: " << key << std::endl);

        //Get some references to simplify the notation
        //the reconstruction bit corresponding to this ReconDataBuffer and encoding space
        IsmrmrdReconBit & rbit = getRBit(recon_data_buffers, key, espace);
        //and the corresponding data buffer for the reference data
        IsmrmrdDataBuffered & dataBuffer = rbit.ref_;
        //this encoding space's xml header info
        ISMRMRD::Encoding & encoding = hdr_.encoding[espace];
        //this bucket's reference stats
        IsmrmrdAcquisitionBucketStats & stats = m1->getObjectPtr()->refstats_[espace];

        //Fill the sampling description for this data buffer
        if(count == 0 )
            fillSamplingDescription(dataBuffer.sampling_, encoding, stats, acqhdr);

        //Make sure that the data storage for this data buffer has been allocated
        //TODO should this check the limits, or should that be done in the stuff function?
        allocateDataArrays(dataBuffer, acqhdr, encoding, stats, true);

        // Stuff the data, header and trajectory into this data buffer
        stuff(it, dataBuffer, encoding);

        count++;
      }


    //Iterate over the imaging data of the bucket
    // this is exactly the same code as for the reference data except for
    // the chunk of the data buffer.
    count = 0;
    for(std::vector<IsmrmrdAcquisitionData>::iterator it = m1->getObjectPtr()->data_.begin();
        it != m1->getObjectPtr()->data_.end(); ++it)
      {
        //Get a reference to the header for this acquisition
        ISMRMRD::AcquisitionHeader & acqhdr = *it->head_->getObjectPtr();

        //Generate the key to the corresponding ReconData buffer
        key = getKey(acqhdr.idx);

        //The storage is based on the encoding space
        uint16_t espace = acqhdr.encoding_space_ref;

        //GDEBUG_STREAM("espace: " << acqhdr.encoding_space_ref << std::endl);
        //GDEBUG_STREAM("slice: " << acqhdr.idx.slice << std::endl);
        //GDEBUG_STREAM("rep: " << acqhdr.idx.repetition << std::endl);
        //GDEBUG_STREAM("k1: " << acqhdr.idx.kspace_encode_step_1 << std::endl);
        //GDEBUG_STREAM("k2: " << acqhdr.idx.kspace_encode_step_2 << std::endl);
        //GDEBUG_STREAM("seg: " << acqhdr.idx.segment << std::endl);
        //GDEBUG_STREAM("key: " << key << std::endl);

        //Get some references to simplify the notation
        //the reconstruction bit corresponding to this ReconDataBuffer and encoding space
        IsmrmrdReconBit & rbit = getRBit(recon_data_buffers, key, espace);
        //and the corresponding data buffer for the imaging data
        IsmrmrdDataBuffered & dataBuffer = rbit.data_;
        //this encoding space's xml header info
        ISMRMRD::Encoding & encoding = hdr_.encoding[espace];
        //this bucket's imaging data stats
        IsmrmrdAcquisitionBucketStats & stats = m1->getObjectPtr()->datastats_[espace];

        //Fill the sampling description for this data buffer
        if (count == 0)
            fillSamplingDescription(dataBuffer.sampling_, encoding, stats, acqhdr);

        //Make sure that the data storage for this data buffer has been allocated
        //TODO should this check the limits, or should that be done in the stuff function?
        allocateDataArrays(dataBuffer, acqhdr, encoding, stats, false);

        // Stuff the data, header and trajectory into this data buffer
        stuff(it, dataBuffer, encoding);

        count++;
      }


    //Send all the ReconData messages
    GDEBUG("End of bucket reached, sending out %d ReconData buffers\n", recon_data_buffers.size());
    for(std::map<size_t, GadgetContainerMessage<IsmrmrdReconData>* >::iterator it = recon_data_buffers.begin(); it != recon_data_buffers.end(); it++)
      {
        //GDEBUG_STREAM("Sending: " << it->first << std::endl);
        if (it->second) {
            if (this->next()->putq(it->second) == -1) {
                it->second->release();
                throw std::runtime_error("Failed to pass bucket down the chain\n");
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
    GDEBUG("BucketToBufferGadget::close\n");

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

  IsmrmrdReconBit & BucketToBufferGadget::getRBit(std::map<size_t, GadgetContainerMessage<IsmrmrdReconData>* > & recon_data_buffers, size_t key, uint16_t espace)
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

  void BucketToBufferGadget::allocateDataArrays(IsmrmrdDataBuffered & dataBuffer, ISMRMRD::AcquisitionHeader & acqhdr, ISMRMRD::Encoding encoding, IsmrmrdAcquisitionBucketStats & stats, bool forref)
  {
    if (dataBuffer.data_.get_number_of_elements() == 0)
      {
        //Allocate the reference data array
        //7D,  fixed order [E0, E1, E2, CHA, N, S, LOC]
        //11D, fixed order [E0, E1, E2, CHA, SLC, PHS, CON, REP, SET, SEG, AVE]
        uint16_t NE0;
        if (encoding.trajectory.compare("cartesian") == 0) {
            // if seperate or external calibration mode, using the acq length for NE0
            if (forref && (encoding.parallelImaging.get().calibrationMode.get() == "separate"
                || encoding.parallelImaging.get().calibrationMode.get() == "external"))
            {
                NE0 = acqhdr.number_of_samples;
            }
            else
            {
                NE0 = encoding.reconSpace.matrixSize.x;
            }
        } else {
            NE0 = acqhdr.number_of_samples - acqhdr.discard_pre - acqhdr.discard_post;
        }

        uint16_t NE1;
        if (encoding.trajectory.compare("cartesian") == 0) {
            if (forref && (encoding.parallelImaging.get().calibrationMode.get() == "separate"
                || encoding.parallelImaging.get().calibrationMode.get() == "external") )
            {
                NE1 = *stats.kspace_encode_step_1.rbegin() - *stats.kspace_encode_step_1.begin() + 1;
            }
            else
            {
                NE1 = encoding.encodedSpace.matrixSize.y;
            }
        } else {
            if (encoding.encodingLimits.kspace_encoding_step_1.is_present()) {
                NE1 = encoding.encodingLimits.kspace_encoding_step_1->maximum - encoding.encodingLimits.kspace_encoding_step_1->minimum + 1;
            } else {
                NE1 = *stats.kspace_encode_step_1.rbegin() - *stats.kspace_encode_step_1.begin() + 1;
            }
        }

        uint16_t NE2;
        if (encoding.trajectory.compare("cartesian") == 0) {
            if (forref && (encoding.parallelImaging.get().calibrationMode.get() == "separate"
                || encoding.parallelImaging.get().calibrationMode.get() == "external"))
            {
                NE2 = encoding.encodingLimits.kspace_encoding_step_2->maximum - encoding.encodingLimits.kspace_encoding_step_2->minimum + 1;
            }
            else
            {
                NE2 = encoding.encodedSpace.matrixSize.z;
            }
        } else {
            if (encoding.encodingLimits.kspace_encoding_step_2.is_present()) {
                NE2 = encoding.encodingLimits.kspace_encoding_step_2->maximum - encoding.encodingLimits.kspace_encoding_step_2->minimum + 1;
            } else {
                NE2 = *stats.kspace_encode_step_2.rbegin() - *stats.kspace_encode_step_2.begin() + 1;
            }
        }

        uint16_t NCHA = acqhdr.active_channels;

        uint16_t NLOC;
        if (split_slices_) {
            NLOC = 1;
        } else {
            if (encoding.encodingLimits.slice.is_present()) {
                NLOC = encoding.encodingLimits.slice->maximum - encoding.encodingLimits.slice->minimum + 1;
            } else {
                NLOC = *stats.slice.rbegin() - *stats.slice.begin() + 1;
            }
        }

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

        //GDEBUG_STREAM("Data dimensions:" << std::endl);
        //GDEBUG_STREAM("   NE0:  " << NE0  << std::endl);
        //GDEBUG_STREAM("   NE1:  " << NE1  << std::endl);
        //GDEBUG_STREAM("   NE2:  " << NE2  << std::endl);
        //GDEBUG_STREAM("   NLOC: " << NLOC << std::endl);
        //GDEBUG_STREAM("   NCHA: " << NCHA << std::endl);
        //GDEBUG_STREAM("   NN:   " << NN   << std::endl);
        //GDEBUG_STREAM("   NS:   " << NS   << std::endl);

        //Allocate the array for the data
        dataBuffer.data_.create(NE0, NE1, NE2, NCHA, NN, NS, NLOC);
        clear(&dataBuffer.data_);

        //Allocate the array for the headers
        dataBuffer.headers_.create(NE1, NE2, NN, NS, NLOC);

        //Allocate the array for the trajectories
        uint16_t TRAJDIM = acqhdr.trajectory_dimensions;
        if (TRAJDIM > 0)
          {
            dataBuffer.trajectory_.create(TRAJDIM, NE0, NE1, NE2, NN, NS, NLOC);
            clear(&dataBuffer.trajectory_);
          }

        //boost::shared_ptr< std::vector<size_t> > dims =  dataBuffer.data_.get_dimensions();
        //GDEBUG_STREAM("NDArray dims: ");
        //for( std::vector<size_t>::const_iterator i = dims->begin(); i != dims->end(); ++i) {
        //    GDEBUG_STREAM(*i << ' ');
        //}
        //GDEBUG_STREAM(std::endl);
      }

  }

  void BucketToBufferGadget::fillSamplingDescription(SamplingDescription & sampling, ISMRMRD::Encoding & encoding, IsmrmrdAcquisitionBucketStats & stats, ISMRMRD::AcquisitionHeader& acqhdr)
  {
    // For cartesian trajectories, assume that any oversampling has been removed.
    if (encoding.trajectory.compare("cartesian") == 0) {
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
    if (encoding.trajectory.compare("cartesian") == 0) {

        size_t num = acqhdr.number_of_samples / 2;

        if (acqhdr.discard_pre == 0 && acqhdr.discard_post==0)
        {
            sampling.sampling_limits_[0].min_ = 0;
            sampling.sampling_limits_[0].max_ = acqhdr.number_of_samples - 1;
        }
        else if (acqhdr.discard_pre>0) // pre zeros
        {
            sampling.sampling_limits_[0].min_ = acqhdr.discard_pre;
            sampling.sampling_limits_[0].max_ = acqhdr.number_of_samples - 1;
        }
        else if (acqhdr.discard_post>0) // post zeros
        {
            sampling.sampling_limits_[0].min_ = 0;
            sampling.sampling_limits_[0].max_ = acqhdr.number_of_samples - acqhdr.discard_post - 1;
        }

        sampling.sampling_limits_[0].center_ = acqhdr.number_of_samples / 2;
    } else {
        sampling.sampling_limits_[0].min_ = 0;
        sampling.sampling_limits_[0].max_ = encoding.encodedSpace.matrixSize.x - 1;
        sampling.sampling_limits_[0].center_ = encoding.encodedSpace.matrixSize.x / 2;
    }

    sampling.sampling_limits_[1].min_ =
        encoding.encodingLimits.kspace_encoding_step_1->minimum;
    sampling.sampling_limits_[1].max_ =
        encoding.encodingLimits.kspace_encoding_step_1->maximum;
    sampling.sampling_limits_[1].center_ =
        encoding.encodingLimits.kspace_encoding_step_1->center;

    sampling.sampling_limits_[2].min_ =
        encoding.encodingLimits.kspace_encoding_step_2->minimum;
    sampling.sampling_limits_[2].max_ =
        encoding.encodingLimits.kspace_encoding_step_2->maximum;
    sampling.sampling_limits_[2].center_ =
        encoding.encodingLimits.kspace_encoding_step_2->center;
  }

  void BucketToBufferGadget::stuff(std::vector<IsmrmrdAcquisitionData>::iterator it, IsmrmrdDataBuffered & dataBuffer, ISMRMRD::Encoding encoding)
  {

    // The acquisition header and data
    ISMRMRD::AcquisitionHeader & acqhdr = *it->head_->getObjectPtr();
    hoNDArray< std::complex<float> > & acqdata = *it->data_->getObjectPtr();
    // we make one for the trajectory down below if we need it

    size_t slice_loc;
    if (split_slices_)
      {
        slice_loc = 0;
      }
    else
      {
        slice_loc = acqhdr.idx.slice;
      }

    //Stuff the data
    uint16_t npts_to_copy = acqhdr.number_of_samples - acqhdr.discard_pre - acqhdr.discard_post;
    long long offset;
    if (encoding.trajectory.compare("cartesian") == 0) {
        offset  = (long long) dataBuffer.sampling_.sampling_limits_[0].center_ - (long long) acqhdr.center_sample;
    } else {
        //TODO what about EPI with asymmetric readouts?
        //TODO any other sort of trajectory?
        offset = 0;
    }
    long long roffset = (long long) dataBuffer.data_.get_size(0) - npts_to_copy - offset;

    //GDEBUG_STREAM("Num_samp: "<< acqhdr.number_of_samples << ", pre: " << acqhdr.discard_pre << ", post" << acqhdr.discard_post << std::endl);
    //std::cout << "Sampling limits: "
    //    << "  min: " << dataBuffer.sampling_.sampling_limits_[0].min_
    //    << "  max: " << dataBuffer.sampling_.sampling_limits_[0].max_
    //    << "  center: " << dataBuffer.sampling_.sampling_limits_[0].center_
    //    << std::endl;
    //GDEBUG_STREAM("npts_to_copy = " << npts_to_copy  << std::endl);
    //GDEBUG_STREAM("offset = " << offset  << std::endl);
    //GDEBUG_STREAM("loffset = " << roffset << std::endl);

    if ((offset < 0) | (roffset < 0) )
      {
        throw std::runtime_error("Acquired reference data does not fit into the reference data buffer.\n");
      }

    std::complex<float> *dataptr;
    uint16_t NCHA = (uint16_t)dataBuffer.data_.get_size(3);
    uint16_t NN = (uint16_t)dataBuffer.data_.get_size(4);
    uint16_t NS = (uint16_t)dataBuffer.data_.get_size(5);

    uint16_t NUsed = (uint16_t)getN(acqhdr.idx);
    if (NUsed >= NN) NUsed = NN - 1;

    uint16_t SUsed = (uint16_t)getS(acqhdr.idx);
    if (SUsed >= NS) SUsed = NS - 1;

    for (uint16_t cha = 0; cha < NCHA; cha++)
      {
        dataptr = & dataBuffer.data_(
            offset, acqhdr.idx.kspace_encode_step_1, acqhdr.idx.kspace_encode_step_2, cha, NUsed, SUsed, slice_loc);


        memcpy(dataptr, &acqdata(acqhdr.discard_pre, cha), sizeof(std::complex<float>)*npts_to_copy);
      }

    //Stuff the header
    dataBuffer.headers_(acqhdr.idx.kspace_encode_step_1,
        acqhdr.idx.kspace_encode_step_2, NUsed, SUsed, slice_loc) = acqhdr;

    //Stuff the trajectory
    if (acqhdr.trajectory_dimensions > 0) {

        hoNDArray< float > & acqtraj = *it->traj_->getObjectPtr();  // TODO do we need to check this?

        float * trajptr;

        trajptr = &dataBuffer.trajectory_(0,
            offset, acqhdr.idx.kspace_encode_step_1, acqhdr.idx.kspace_encode_step_2, NUsed, SUsed, slice_loc);

        memcpy(trajptr, &acqtraj(0, acqhdr.discard_pre), sizeof(float)*npts_to_copy*acqhdr.trajectory_dimensions);

    }
  }

  GADGET_FACTORY_DECLARE(BucketToBufferGadget)

}

