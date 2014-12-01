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

    ignore_segment_  = this->get_bool_value("ignore_segment");
    GADGET_DEBUG2("IGNORE SEGMENT IS: %b\n", ignore_segment_);

    // keep a copy of the deserialized ismrmrd xml header for runtime
    ISMRMRD::deserialize(mb->rd_ptr(), hdr_);
    
    return GADGET_OK;
}

int BucketToBufferGadget
::process(GadgetContainerMessage<IsmrmrdAcquisitionBucket>* m1)
{

    size_t key;
    std::map<size_t, GadgetContainerMessage<IsmrmrdReconData>* > recon_data_buffers;
 
    //GADGET_DEBUG1("BucketToBufferGadget::process\n");

    //Some information about the bucket
    //std::cout << "The Reference part: " << m1->getObjectPtr()->refstats_.size() << std::endl;
    //std::cout << "   nslices: " << m1->getObjectPtr()->refstats_[0].slice.size() << std::endl;
    //for (int e=0; e<m1->getObjectPtr()->refstats_.size() ; e++) {
    //    for (std::set<uint16_t>::iterator it = m1->getObjectPtr()->refstats_[e].kspace_encode_step_1.begin();
    //         it != m1->getObjectPtr()->refstats_[e].kspace_encode_step_1.end(); ++it) {            
    //        std::cout << "   K1: " <<  *it << std::endl;
    //    }
    //}
    //std::cout << "The data part: " << m1->getObjectPtr()->datastats_.size() << std::endl;
    //std::cout << "   nslices: " << m1->getObjectPtr()->datastats_[0].slice.size() << std::endl;
    //for (int e=0; e<m1->getObjectPtr()->datastats_.size() ; e++) {
    //    for (std::set<uint16_t>::iterator it = m1->getObjectPtr()->datastats_[e].kspace_encode_step_1.begin();
    //         it != m1->getObjectPtr()->datastats_[e].kspace_encode_step_1.end(); ++it) {            
    //        std::cout << "   K1: " <<  *it << std::endl;
    //    }
    //}

    //Iterate over the reference data of the bucket
    for(std::vector<IsmrmrdAcquisitionData>::iterator it = m1->getObjectPtr()->ref_.begin();
        it != m1->getObjectPtr()->ref_.end(); ++it)
    {
        ISMRMRD::AcquisitionHeader & acqhdr = *it->head_->getObjectPtr();
        hoNDArray< std::complex<float> > & acqdata = *it->data_->getObjectPtr();
        
        uint16_t espace = acqhdr.encoding_space_ref;

        //Generate the key to the corresponding ReconData buffer
        key = getKey(acqhdr.idx);
        
        //std::cout << "espace: " << acqhdr.encoding_space_ref << std::endl;
        //std::cout << "slice: " << acqhdr.idx.slice << std::endl;
        //std::cout << "rep: " << acqhdr.idx.repetition << std::endl;
        //std::cout << "k1: " << acqhdr.idx.kspace_encode_step_1 << std::endl;
        //std::cout << "k2: " << acqhdr.idx.kspace_encode_step_2 << std::endl;
        //std::cout << "seg: " << acqhdr.idx.segment << std::endl;
        //std::cout << "key: " << key << std::endl;
        
        //Look up the corresponding ReconData buffer
        if (recon_data_buffers.find(key) == recon_data_buffers.end())
        {
            //ReconData buffer does not exist, create it
            recon_data_buffers[key] = new GadgetContainerMessage<IsmrmrdReconData>;
        }

        //Look up the DataBuffered entry corresponding to this encoding space
        // create if needed and set the fields of view and matrix size
        //std::cout << "RBIT " << recon_data_buffers[key]->getObjectPtr()->rbit_.size() << std::endl;
        if ( recon_data_buffers[key]->getObjectPtr()->rbit_.size() < (espace+1) )
        {
            recon_data_buffers[key]->getObjectPtr()->rbit_.resize(espace+1);
        }
        
        // For cartesian trajectories, assume that any oversampling has been removed.
        if (hdr_.encoding[espace].trajectory.compare("cartesian") == 0) {
            recon_data_buffers[key]->getObjectPtr()->rbit_[espace].data_.sampling_.encoded_FOV_[0] = hdr_.encoding[espace].reconSpace.fieldOfView_mm.x;
            recon_data_buffers[key]->getObjectPtr()->rbit_[espace].data_.sampling_.encoded_matrix_[0] = hdr_.encoding[espace].reconSpace.matrixSize.x;
        } else {
            recon_data_buffers[key]->getObjectPtr()->rbit_[espace].data_.sampling_.encoded_FOV_[0] = hdr_.encoding[espace].encodedSpace.fieldOfView_mm.x;
            recon_data_buffers[key]->getObjectPtr()->rbit_[espace].data_.sampling_.encoded_matrix_[0] = hdr_.encoding[espace].encodedSpace.matrixSize.x;
        }
        
        recon_data_buffers[key]->getObjectPtr()->rbit_[espace].ref_.sampling_.encoded_FOV_[1] = hdr_.encoding[espace].encodedSpace.fieldOfView_mm.y;
        recon_data_buffers[key]->getObjectPtr()->rbit_[espace].ref_.sampling_.encoded_FOV_[2] = hdr_.encoding[espace].encodedSpace.fieldOfView_mm.z;
        
        recon_data_buffers[key]->getObjectPtr()->rbit_[espace].ref_.sampling_.encoded_matrix_[1] = hdr_.encoding[espace].encodedSpace.matrixSize.y;
        recon_data_buffers[key]->getObjectPtr()->rbit_[espace].ref_.sampling_.encoded_matrix_[2] = hdr_.encoding[espace].encodedSpace.matrixSize.z;
        
        recon_data_buffers[key]->getObjectPtr()->rbit_[espace].ref_.sampling_.recon_FOV_[0] = hdr_.encoding[espace].reconSpace.fieldOfView_mm.x;
        recon_data_buffers[key]->getObjectPtr()->rbit_[espace].ref_.sampling_.recon_FOV_[1] = hdr_.encoding[espace].reconSpace.fieldOfView_mm.y;
        recon_data_buffers[key]->getObjectPtr()->rbit_[espace].ref_.sampling_.recon_FOV_[2] = hdr_.encoding[espace].reconSpace.fieldOfView_mm.z;
        
        recon_data_buffers[key]->getObjectPtr()->rbit_[espace].ref_.sampling_.recon_matrix_[0] = hdr_.encoding[espace].reconSpace.matrixSize.x;
        recon_data_buffers[key]->getObjectPtr()->rbit_[espace].ref_.sampling_.recon_matrix_[1] = hdr_.encoding[espace].reconSpace.matrixSize.y;
        recon_data_buffers[key]->getObjectPtr()->rbit_[espace].ref_.sampling_.recon_matrix_[2] = hdr_.encoding[espace].reconSpace.matrixSize.z;
        
        // For cartesian trajectories, assume that any oversampling has been removed.
        if (hdr_.encoding[espace].trajectory.compare("cartesian") == 0) {
            recon_data_buffers[key]->getObjectPtr()->rbit_[espace].data_.sampling_.sampling_limits_[0].min_ = 0;
            recon_data_buffers[key]->getObjectPtr()->rbit_[espace].data_.sampling_.sampling_limits_[0].max_ = hdr_.encoding[espace].reconSpace.matrixSize.x - 1;
            recon_data_buffers[key]->getObjectPtr()->rbit_[espace].data_.sampling_.sampling_limits_[0].center_ = hdr_.encoding[espace].reconSpace.matrixSize.x / 2;
        } else {
            recon_data_buffers[key]->getObjectPtr()->rbit_[espace].data_.sampling_.sampling_limits_[0].min_ = 0;
            recon_data_buffers[key]->getObjectPtr()->rbit_[espace].data_.sampling_.sampling_limits_[0].max_ = hdr_.encoding[espace].encodedSpace.matrixSize.x - 1;            
            recon_data_buffers[key]->getObjectPtr()->rbit_[espace].data_.sampling_.sampling_limits_[0].center_ = hdr_.encoding[espace].encodedSpace.matrixSize.x / 2;
        }
        
        recon_data_buffers[key]->getObjectPtr()->rbit_[espace].ref_.sampling_.sampling_limits_[1].min_ =
                hdr_.encoding[espace].encodingLimits.kspace_encoding_step_1->minimum;
        recon_data_buffers[key]->getObjectPtr()->rbit_[espace].ref_.sampling_.sampling_limits_[1].max_ =
                hdr_.encoding[espace].encodingLimits.kspace_encoding_step_1->maximum;
        recon_data_buffers[key]->getObjectPtr()->rbit_[espace].ref_.sampling_.sampling_limits_[1].center_ =
                hdr_.encoding[espace].encodingLimits.kspace_encoding_step_1->center;
        
        recon_data_buffers[key]->getObjectPtr()->rbit_[espace].ref_.sampling_.sampling_limits_[2].min_ =
                hdr_.encoding[espace].encodingLimits.kspace_encoding_step_2->minimum;
        recon_data_buffers[key]->getObjectPtr()->rbit_[espace].ref_.sampling_.sampling_limits_[2].max_ =
                hdr_.encoding[espace].encodingLimits.kspace_encoding_step_2->maximum;
        recon_data_buffers[key]->getObjectPtr()->rbit_[espace].ref_.sampling_.sampling_limits_[2].center_ =
                hdr_.encoding[espace].encodingLimits.kspace_encoding_step_2->center;

        if (recon_data_buffers[key]->getObjectPtr()->rbit_[espace].ref_.data_.get_number_of_elements() == 0)
        {
            //Allocate the reference data array
            //7D,  fixed order [RO, E1, E2, CHA, SLC, N, S]
            //11D, fixed order [RO, E1, E2, CHA, SLC, PHS, CON, REP, SET, SEG, AVE]
            uint16_t NRO;
            if (hdr_.encoding[espace].trajectory.compare("cartesian") == 0) {
                NRO = hdr_.encoding[espace].reconSpace.matrixSize.x;
            } else {
                NRO = acqhdr.number_of_samples - acqhdr.discard_pre - acqhdr.discard_post;
            }
            
            uint16_t NE1;
            if (hdr_.encoding[espace].trajectory.compare("cartesian") == 0) {
                NE1 = hdr_.encoding[espace].encodedSpace.matrixSize.y;
            } else {
                if (hdr_.encoding[espace].encodingLimits.kspace_encoding_step_1.is_present()) {
                NE1 = hdr_.encoding[espace].encodingLimits.kspace_encoding_step_1->maximum - hdr_.encoding[espace].encodingLimits.kspace_encoding_step_1->minimum + 1;
                } else {
                    NE1 = *m1->getObjectPtr()->refstats_[espace].kspace_encode_step_1.rbegin() - *m1->getObjectPtr()->refstats_[espace].kspace_encode_step_1.begin() + 1;
                }
            }
            
            uint16_t NE2;
            if (hdr_.encoding[espace].trajectory.compare("cartesian") == 0) {
                NE2 = hdr_.encoding[espace].encodedSpace.matrixSize.z;
            } else {
                if (hdr_.encoding[espace].encodingLimits.kspace_encoding_step_2.is_present()) {
                    NE2 = hdr_.encoding[espace].encodingLimits.kspace_encoding_step_2->maximum - hdr_.encoding[espace].encodingLimits.kspace_encoding_step_2->minimum + 1;
                } else {
                    NE2 = *m1->getObjectPtr()->refstats_[espace].kspace_encode_step_2.rbegin() - *m1->getObjectPtr()->refstats_[espace].kspace_encode_step_2.begin() + 1;
                }
            }
            
            uint16_t NCHA = acqhdr.active_channels;
            
            uint16_t NSLC;
            if (split_slices_) {
                NSLC = 1;
            } else {
                if (hdr_.encoding[espace].encodingLimits.slice.is_present()) {
                    NSLC = hdr_.encoding[espace].encodingLimits.slice->maximum - hdr_.encoding[espace].encodingLimits.slice->minimum + 1;
                } else {
                    NSLC = *m1->getObjectPtr()->refstats_[espace].slice.rbegin() - *m1->getObjectPtr()->refstats_[espace].slice.begin() + 1;
                }
            }
            
            uint16_t NN;
            switch (N_) {
                case PHASE:
                    NN = *m1->getObjectPtr()->refstats_[espace].phase.rbegin() - *m1->getObjectPtr()->refstats_[espace].phase.begin() + 1;
                    break;
                case CONTRAST:
                    NN = *m1->getObjectPtr()->refstats_[espace].contrast.rbegin() - *m1->getObjectPtr()->refstats_[espace].contrast.begin() + 1;
                    break;
                case REPETITION:
                    NN = *m1->getObjectPtr()->refstats_[espace].repetition.rbegin() - *m1->getObjectPtr()->refstats_[espace].repetition.begin() + 1;
                    break;
                case SET:
                    NN = *m1->getObjectPtr()->refstats_[espace].set.rbegin() - *m1->getObjectPtr()->refstats_[espace].set.begin() + 1;
                    break;
                case SEGMENT:
                    NN = *m1->getObjectPtr()->refstats_[espace].segment.rbegin() - *m1->getObjectPtr()->refstats_[espace].segment.begin() + 1;
                    break;
                case AVERAGE:
                    NN = *m1->getObjectPtr()->refstats_[espace].average.rbegin() - *m1->getObjectPtr()->refstats_[espace].average.begin() + 1;
                    break;
                default:
                    NN = 1;
            }
            
            uint16_t NS;
            switch (S_) {
                case PHASE:
                    NS = *m1->getObjectPtr()->refstats_[espace].phase.rbegin() - *m1->getObjectPtr()->refstats_[espace].phase.begin() + 1;
                    break;
                case CONTRAST:
                    NS = *m1->getObjectPtr()->refstats_[espace].contrast.rbegin() - *m1->getObjectPtr()->refstats_[espace].contrast.begin() + 1;
                    break;
                case REPETITION:
                    NS = *m1->getObjectPtr()->refstats_[espace].repetition.rbegin() - *m1->getObjectPtr()->refstats_[espace].repetition.begin() + 1;
                    break;
                case SET:
                    NS = *m1->getObjectPtr()->refstats_[espace].set.rbegin() - *m1->getObjectPtr()->refstats_[espace].set.begin() + 1;
                    break;
                case SEGMENT:
                    NS = *m1->getObjectPtr()->refstats_[espace].segment.rbegin() - *m1->getObjectPtr()->refstats_[espace].segment.begin() + 1;
                    break;
                case AVERAGE:
                    NS = *m1->getObjectPtr()->refstats_[espace].average.rbegin() - *m1->getObjectPtr()->refstats_[espace].average.begin() + 1;
                    break;
                default:
                    NS = 1;
            }

            //std::cout << "Reference dimensions:" << std::endl;
            //std::cout << "   NRO:  " << NRO  << std::endl;
            //std::cout << "   NE1:  " << NE1  << std::endl;
            //std::cout << "   NE2:  " << NE2  << std::endl;
            //std::cout << "   NSLC: " << NSLC << std::endl;
            //std::cout << "   NCHA: " << NCHA << std::endl;
            //std::cout << "   NN:   " << NN   << std::endl;
            //std::cout << "   NS:   " << NS   << std::endl;

            //Allocate the array for the data
            recon_data_buffers[key]->getObjectPtr()->rbit_[espace].ref_.data_.create(NRO, NE1, NE2, NCHA, NSLC, NN, NS);
            recon_data_buffers[key]->getObjectPtr()->rbit_[espace].ref_.data_.fill(std::complex<float>(0.0f,0.0f));

            //Allocate the array for the headers
            recon_data_buffers[key]->getObjectPtr()->rbit_[espace].ref_.headers_.create(NE1, NE2, NSLC, NN, NS);
            
            //Allocate the array for the trajectories
            uint16_t TRAJDIM = acqhdr.trajectory_dimensions;
            if (TRAJDIM > 0)
            {
                recon_data_buffers[key]->getObjectPtr()->rbit_[espace].ref_.trajectory_.create(NRO, TRAJDIM, NE1, NE2, NSLC, NN, NS);
                recon_data_buffers[key]->getObjectPtr()->rbit_[espace].ref_.trajectory_.fill(0.0f);
            }
            
            //boost::shared_ptr< std::vector<size_t> > dims =  recon_data_buffers[key]->getObjectPtr()->rbit_[espace].ref_.data_.get_dimensions();
            //std::cout << "Ref NDArray dims: ";
            //for( std::vector<size_t>::const_iterator i = dims->begin(); i != dims->end(); ++i) {
            //    std::cout << *i << ' ';
            //}
            //std::cout << std::endl;
        }
        
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
        if (hdr_.encoding[espace].trajectory.compare("cartesian") == 0) {
            offset  = (long long) recon_data_buffers[key]->getObjectPtr()->rbit_[espace].ref_.sampling_.sampling_limits_[0].center_ - (long long) acqhdr.center_sample;
        } else {
            //TODO what about EPI with asymmetric readouts?
            //TODO any other sort of trajectory?
            offset = 0;
        }
        long long roffset = (long long) recon_data_buffers[key]->getObjectPtr()->rbit_[espace].ref_.data_.get_size(0) - npts_to_copy - offset;

        //std::cout << "REFERENCE" << std::endl;
        //std::cout << "Num_samp: "<< acqhdr.number_of_samples << ", pre: " << acqhdr.discard_pre << ", post" << acqhdr.discard_post << std::endl;        
        //std::cout << "Sampling limits: " 
        //          << "  min: " << recon_data_buffers[key]->getObjectPtr()->rbit_[espace].ref_.sampling_.sampling_limits_[0].min_
        //          << "  max: " << recon_data_buffers[key]->getObjectPtr()->rbit_[espace].ref_.sampling_.sampling_limits_[0].max_
        //          << "  center: " << recon_data_buffers[key]->getObjectPtr()->rbit_[espace].ref_.sampling_.sampling_limits_[0].center_
        //          << std::endl;
        //std::cout << "npts_to_copy = " << npts_to_copy  << std::endl;
        //std::cout << "offset = " << offset  << std::endl;
        //std::cout << "loffset = " << roffset << std::endl;
        
        if ((offset < 0) | (roffset < 0) )
        {
            throw std::runtime_error("Acquired reference data does not fit into the reference data buffer.\n");
        }
        
        std::complex<float> *dataptr;
        uint16_t NCHA = recon_data_buffers[key]->getObjectPtr()->rbit_[espace].ref_.data_.get_size(3);
        for (uint16_t cha = 0; cha < NCHA; cha++)
        {
            dataptr = & recon_data_buffers[key]->getObjectPtr()->rbit_[espace].ref_.data_(
                offset, acqhdr.idx.kspace_encode_step_1, acqhdr.idx.kspace_encode_step_2, cha, slice_loc, getN(acqhdr.idx),  getS(acqhdr.idx));
            memcpy(dataptr, &acqdata(acqhdr.discard_pre, cha), sizeof(std::complex<float>)*npts_to_copy);
        }

        //Stuff the header
        recon_data_buffers[key]->getObjectPtr()->rbit_[espace].ref_.headers_(acqhdr.idx.kspace_encode_step_1,
                                 acqhdr.idx.kspace_encode_step_2, slice_loc, getN(acqhdr.idx),  getS(acqhdr.idx)) = acqhdr;

        //Stuff the trajectory
        if (acqhdr.trajectory_dimensions > 0) {
            float * trajptr;
            hoNDArray< float > & acqtraj = *it->traj_->getObjectPtr();  // TODO do we need to check this?            
            for (uint16_t tdim = 0; tdim < acqhdr.trajectory_dimensions; tdim++) {
                trajptr = &recon_data_buffers[key]->getObjectPtr()->rbit_[espace].ref_.trajectory_(
                    offset, tdim, acqhdr.idx.kspace_encode_step_1, acqhdr.idx.kspace_encode_step_2, slice_loc, getN(acqhdr.idx),  getS(acqhdr.idx));
                memcpy(trajptr, & acqtraj(acqhdr.discard_pre, tdim), sizeof(float)*npts_to_copy);
            }
        }
        
    }


    
    //Iterate over the imaging data of the bucket
    for(std::vector<IsmrmrdAcquisitionData>::iterator it = m1->getObjectPtr()->data_.begin();
        it != m1->getObjectPtr()->data_.end(); ++it)
    {
        ISMRMRD::AcquisitionHeader & acqhdr = *it->head_->getObjectPtr();
        hoNDArray< std::complex<float> > & acqdata = *it->data_->getObjectPtr();
        
        uint16_t espace = acqhdr.encoding_space_ref;

        //Generate the key to the corresponding ReconData buffer
        key = getKey(acqhdr.idx);
        
        //std::cout << "espace: " << acqhdr.encoding_space_ref << std::endl;
        //std::cout << "slice: " << acqhdr.idx.slice << std::endl;
        //std::cout << "rep: " << acqhdr.idx.repetition << std::endl;
        //std::cout << "k1: " << acqhdr.idx.kspace_encode_step_1 << std::endl;
        //std::cout << "k2: " << acqhdr.idx.kspace_encode_step_2 << std::endl;
        //std::cout << "seg: " << acqhdr.idx.segment << std::endl;
        //std::cout << "key: " << key << std::endl;
        
        //Look up the corresponding ReconData buffer
        if (recon_data_buffers.find(key) == recon_data_buffers.end())
        {
            //ReconData buffer does not exist, create it
            recon_data_buffers[key] = new GadgetContainerMessage<IsmrmrdReconData>;
        }

        //Look up the DataBuffered entry corresponding to this encoding space
        // create if needed and set the fields of view and matrix size
        //std::cout << "RBIT " << recon_data_buffers[key]->getObjectPtr()->rbit_.size() << std::endl;
        if ( recon_data_buffers[key]->getObjectPtr()->rbit_.size() < (espace+1) )
        {
            recon_data_buffers[key]->getObjectPtr()->rbit_.resize(espace+1);
        }
        
        // For cartesian trajectories, assume that any oversampling has been removed.
        if (hdr_.encoding[espace].trajectory.compare("cartesian") == 0) {
            recon_data_buffers[key]->getObjectPtr()->rbit_[espace].data_.sampling_.encoded_FOV_[0] = hdr_.encoding[espace].reconSpace.fieldOfView_mm.x;
            recon_data_buffers[key]->getObjectPtr()->rbit_[espace].data_.sampling_.encoded_matrix_[0] = hdr_.encoding[espace].reconSpace.matrixSize.x;
        } else {
            recon_data_buffers[key]->getObjectPtr()->rbit_[espace].data_.sampling_.encoded_FOV_[0] = hdr_.encoding[espace].encodedSpace.fieldOfView_mm.x;
            recon_data_buffers[key]->getObjectPtr()->rbit_[espace].data_.sampling_.encoded_matrix_[0] = hdr_.encoding[espace].encodedSpace.matrixSize.x;
        }

        recon_data_buffers[key]->getObjectPtr()->rbit_[espace].data_.sampling_.encoded_FOV_[1] = hdr_.encoding[espace].encodedSpace.fieldOfView_mm.y;
        recon_data_buffers[key]->getObjectPtr()->rbit_[espace].data_.sampling_.encoded_FOV_[2] = hdr_.encoding[espace].encodedSpace.fieldOfView_mm.z;
        
        recon_data_buffers[key]->getObjectPtr()->rbit_[espace].data_.sampling_.encoded_matrix_[1] = hdr_.encoding[espace].encodedSpace.matrixSize.y;
        recon_data_buffers[key]->getObjectPtr()->rbit_[espace].data_.sampling_.encoded_matrix_[2] = hdr_.encoding[espace].encodedSpace.matrixSize.z;
        
        recon_data_buffers[key]->getObjectPtr()->rbit_[espace].data_.sampling_.recon_FOV_[0] = hdr_.encoding[espace].reconSpace.fieldOfView_mm.x;
        recon_data_buffers[key]->getObjectPtr()->rbit_[espace].data_.sampling_.recon_FOV_[1] = hdr_.encoding[espace].reconSpace.fieldOfView_mm.y;
        recon_data_buffers[key]->getObjectPtr()->rbit_[espace].data_.sampling_.recon_FOV_[2] = hdr_.encoding[espace].reconSpace.fieldOfView_mm.z;
        
        recon_data_buffers[key]->getObjectPtr()->rbit_[espace].data_.sampling_.recon_matrix_[0] = hdr_.encoding[espace].reconSpace.matrixSize.x;
        recon_data_buffers[key]->getObjectPtr()->rbit_[espace].data_.sampling_.recon_matrix_[1] = hdr_.encoding[espace].reconSpace.matrixSize.y;
        recon_data_buffers[key]->getObjectPtr()->rbit_[espace].data_.sampling_.recon_matrix_[2] = hdr_.encoding[espace].reconSpace.matrixSize.z;
        
        // For cartesian trajectories, assume that any oversampling has been removed.
        if (hdr_.encoding[espace].trajectory.compare("cartesian") == 0) {
            recon_data_buffers[key]->getObjectPtr()->rbit_[espace].data_.sampling_.sampling_limits_[0].min_ = 0;
            recon_data_buffers[key]->getObjectPtr()->rbit_[espace].data_.sampling_.sampling_limits_[0].max_ = hdr_.encoding[espace].reconSpace.matrixSize.x - 1;
            recon_data_buffers[key]->getObjectPtr()->rbit_[espace].data_.sampling_.sampling_limits_[0].center_ = hdr_.encoding[espace].reconSpace.matrixSize.x / 2;
        } else {
            recon_data_buffers[key]->getObjectPtr()->rbit_[espace].data_.sampling_.sampling_limits_[0].min_ = 0;
            recon_data_buffers[key]->getObjectPtr()->rbit_[espace].data_.sampling_.sampling_limits_[0].max_ = hdr_.encoding[espace].encodedSpace.matrixSize.x - 1;            
            recon_data_buffers[key]->getObjectPtr()->rbit_[espace].data_.sampling_.sampling_limits_[0].center_ = hdr_.encoding[espace].encodedSpace.matrixSize.x / 2;
        }
        
        recon_data_buffers[key]->getObjectPtr()->rbit_[espace].data_.sampling_.sampling_limits_[1].min_ =
                hdr_.encoding[espace].encodingLimits.kspace_encoding_step_1->minimum;
        recon_data_buffers[key]->getObjectPtr()->rbit_[espace].data_.sampling_.sampling_limits_[1].max_ =
                hdr_.encoding[espace].encodingLimits.kspace_encoding_step_1->maximum;
        recon_data_buffers[key]->getObjectPtr()->rbit_[espace].data_.sampling_.sampling_limits_[1].center_ =
                hdr_.encoding[espace].encodingLimits.kspace_encoding_step_1->center;
        
        recon_data_buffers[key]->getObjectPtr()->rbit_[espace].data_.sampling_.sampling_limits_[2].min_ =
                hdr_.encoding[espace].encodingLimits.kspace_encoding_step_2->minimum;
        recon_data_buffers[key]->getObjectPtr()->rbit_[espace].data_.sampling_.sampling_limits_[2].max_ =
                hdr_.encoding[espace].encodingLimits.kspace_encoding_step_2->maximum;
        recon_data_buffers[key]->getObjectPtr()->rbit_[espace].data_.sampling_.sampling_limits_[2].center_ =
                hdr_.encoding[espace].encodingLimits.kspace_encoding_step_2->center;

        if (recon_data_buffers[key]->getObjectPtr()->rbit_[espace].data_.data_.get_number_of_elements() == 0)
        {
            //Allocate the reference data array
            //7D,  fixed order [RO, E1, E2, CHA, SLC, N, S]
            //11D, fixed order [RO, E1, E2, CHA, SLC, PHS, CON, REP, SET, SEG, AVE]
            uint16_t NRO;
            if (hdr_.encoding[espace].trajectory.compare("cartesian") == 0) {
                NRO = hdr_.encoding[espace].reconSpace.matrixSize.x;
            } else {
                NRO = acqhdr.number_of_samples - acqhdr.discard_pre - acqhdr.discard_post;
            }
            
            uint16_t NE1;
            if (hdr_.encoding[espace].trajectory.compare("cartesian") == 0) {
                NE1 = hdr_.encoding[espace].encodedSpace.matrixSize.y;
            } else {
                if (hdr_.encoding[espace].encodingLimits.kspace_encoding_step_1.is_present()) {
                    NE1 = hdr_.encoding[espace].encodingLimits.kspace_encoding_step_1->maximum - hdr_.encoding[espace].encodingLimits.kspace_encoding_step_1->minimum + 1;
                } else {
                    NE1 = *m1->getObjectPtr()->datastats_[espace].kspace_encode_step_1.rbegin() - *m1->getObjectPtr()->datastats_[espace].kspace_encode_step_1.begin() + 1;
                }
            }
            
            uint16_t NE2;
            if (hdr_.encoding[espace].trajectory.compare("cartesian") == 0) {
                NE2 = hdr_.encoding[espace].encodedSpace.matrixSize.z;
            } else {
                if (hdr_.encoding[espace].encodingLimits.kspace_encoding_step_2.is_present()) {
                    NE2 = hdr_.encoding[espace].encodingLimits.kspace_encoding_step_2->maximum - hdr_.encoding[espace].encodingLimits.kspace_encoding_step_2->minimum + 1;
                } else {
                    NE2 = *m1->getObjectPtr()->datastats_[espace].kspace_encode_step_2.rbegin() - *m1->getObjectPtr()->datastats_[espace].kspace_encode_step_2.begin() + 1;
                }
            }
            
            uint16_t NCHA = acqhdr.active_channels;
            
            uint16_t NSLC;
            if (split_slices_) {
                NSLC = 1;
            } else {
                if (hdr_.encoding[espace].encodingLimits.slice.is_present()) {
                    NSLC = hdr_.encoding[espace].encodingLimits.slice->maximum - hdr_.encoding[espace].encodingLimits.slice->minimum + 1;
                } else {
                    NSLC = *m1->getObjectPtr()->datastats_[espace].slice.rbegin() - *m1->getObjectPtr()->datastats_[espace].slice.begin() + 1;
                }
            }
            
            uint16_t NN;
            switch (N_) {
                case PHASE:
                    NN = *m1->getObjectPtr()->datastats_[espace].phase.rbegin() - *m1->getObjectPtr()->datastats_[espace].phase.begin() + 1;
                    break;
                case CONTRAST:
                    NN = *m1->getObjectPtr()->datastats_[espace].contrast.rbegin() - *m1->getObjectPtr()->datastats_[espace].contrast.begin() + 1;
                    break;
                case REPETITION:
                    NN = *m1->getObjectPtr()->datastats_[espace].repetition.rbegin() - *m1->getObjectPtr()->datastats_[espace].repetition.begin() + 1;
                    break;
                case SET:
                    NN = *m1->getObjectPtr()->datastats_[espace].set.rbegin() - *m1->getObjectPtr()->datastats_[espace].set.begin() + 1;
                    break;
                case SEGMENT:
                    NN = *m1->getObjectPtr()->datastats_[espace].segment.rbegin() - *m1->getObjectPtr()->datastats_[espace].segment.begin() + 1;
                    break;
                case AVERAGE:
                    NN = *m1->getObjectPtr()->datastats_[espace].average.rbegin() - *m1->getObjectPtr()->datastats_[espace].average.begin() + 1;
                    break;
                default:
                    NN = 1;
            }
            
            uint16_t NS;
            switch (S_) {
                case PHASE:
                    NS = *m1->getObjectPtr()->datastats_[espace].phase.rbegin() - *m1->getObjectPtr()->datastats_[espace].phase.begin() + 1;
                    break;
                case CONTRAST:
                    NS = *m1->getObjectPtr()->datastats_[espace].contrast.rbegin() - *m1->getObjectPtr()->datastats_[espace].contrast.begin() + 1;
                    break;
                case REPETITION:
                    NS = *m1->getObjectPtr()->datastats_[espace].repetition.rbegin() - *m1->getObjectPtr()->datastats_[espace].repetition.begin() + 1;
                    break;
                case SET:
                    NS = *m1->getObjectPtr()->datastats_[espace].set.rbegin() - *m1->getObjectPtr()->datastats_[espace].set.begin() + 1;
                    break;
                case SEGMENT:
                    NS = *m1->getObjectPtr()->datastats_[espace].segment.rbegin() - *m1->getObjectPtr()->datastats_[espace].segment.begin() + 1;
                    break;
                case AVERAGE:
                    NS = *m1->getObjectPtr()->datastats_[espace].average.rbegin() - *m1->getObjectPtr()->datastats_[espace].average.begin() + 1;
                    break;
                default:
                    NS = 1;
            }

            //std::cout << "Data dimensions:" << std::endl;
            //std::cout << "   NRO:  " << NRO  << std::endl;
            //std::cout << "   NE1:  " << NE1  << std::endl;
            //std::cout << "   NE2:  " << NE2  << std::endl;
            //std::cout << "   NSLC: " << NSLC << std::endl;
            //std::cout << "   NCHA: " << NCHA << std::endl;
            //std::cout << "   NN:   " << NN   << std::endl;
            //std::cout << "   NS:   " << NS   << std::endl;

            //Allocate the array for the data
            recon_data_buffers[key]->getObjectPtr()->rbit_[espace].data_.data_.create(NRO, NE1, NE2, NCHA, NSLC, NN, NS);
            recon_data_buffers[key]->getObjectPtr()->rbit_[espace].data_.data_.fill(std::complex<float>(0.0f,0.0f));

            //Allocate the array for the headers
            recon_data_buffers[key]->getObjectPtr()->rbit_[espace].data_.headers_.create(NE1, NE2, NSLC, NN, NS);
            
            //Allocate the array for the trajectories
            uint16_t TRAJDIM = acqhdr.trajectory_dimensions;
            if (TRAJDIM > 0)
            {
                recon_data_buffers[key]->getObjectPtr()->rbit_[espace].data_.trajectory_.create(NRO, TRAJDIM, NE1, NE2, NSLC, NN, NS);
                recon_data_buffers[key]->getObjectPtr()->rbit_[espace].data_.trajectory_.fill(0.0f);
            }
        }
        
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
        if (hdr_.encoding[espace].trajectory.compare("cartesian") == 0) {
            offset  = (long long) recon_data_buffers[key]->getObjectPtr()->rbit_[espace].data_.sampling_.sampling_limits_[0].center_ - (long long) acqhdr.center_sample;
        } else {
            //TODO what about EPI with asymmetric readouts?
            //TODO any other sort of trajectory?
            offset = 0;
        }
        long long roffset = (long long) recon_data_buffers[key]->getObjectPtr()->rbit_[espace].data_.data_.get_size(0) - npts_to_copy - offset;

        //std::cout << "DATA" << std::endl;
        //std::cout << "Num_samp: "<< acqhdr.number_of_samples << ", pre: " << acqhdr.discard_pre << ", post" << acqhdr.discard_post << std::endl;        
        //std::cout << "Sampling limits: " 
        //          << "  min: " << recon_data_buffers[key]->getObjectPtr()->rbit_[espace].data_.sampling_.sampling_limits_[0].min_
        //          << "  max: " << recon_data_buffers[key]->getObjectPtr()->rbit_[espace].data_.sampling_.sampling_limits_[0].max_
        //          << "  center: " << recon_data_buffers[key]->getObjectPtr()->rbit_[espace].data_.sampling_.sampling_limits_[0].center_
        //          << std::endl;
        //std::cout << "npts_to_copy = " << npts_to_copy  << std::endl;
        //std::cout << "offset = " << offset  << std::endl;
        //std::cout << "loffset = " << roffset << std::endl;
                
        if ((offset < 0) | (roffset < 0) )
        {
            throw std::runtime_error("Acquired reference data does not fit into the reference data buffer.\n");
        }
        
        std::complex<float> *dataptr;
        uint16_t NCHA = recon_data_buffers[key]->getObjectPtr()->rbit_[espace].data_.data_.get_size(3);
        for (uint16_t cha = 0; cha < NCHA; cha++)
        {
            dataptr = & recon_data_buffers[key]->getObjectPtr()->rbit_[espace].data_.data_(
                offset, acqhdr.idx.kspace_encode_step_1, acqhdr.idx.kspace_encode_step_2, cha, slice_loc, getN(acqhdr.idx),  getS(acqhdr.idx));
            memcpy(dataptr, &acqdata(acqhdr.discard_pre, cha), sizeof(std::complex<float>)*npts_to_copy);
        }

        //Stuff the header
        recon_data_buffers[key]->getObjectPtr()->rbit_[espace].data_.headers_(acqhdr.idx.kspace_encode_step_1,
                                 acqhdr.idx.kspace_encode_step_2, slice_loc, getN(acqhdr.idx),  getS(acqhdr.idx)) = acqhdr;

        //Stuff the trajectory
        if (acqhdr.trajectory_dimensions > 0) {
            float * trajptr;
            hoNDArray< float > & acqtraj = *it->traj_->getObjectPtr();            
            for (uint16_t tdim = 0; tdim < acqhdr.trajectory_dimensions; tdim++) {
                trajptr = &recon_data_buffers[key]->getObjectPtr()->rbit_[espace].data_.trajectory_(
                    offset, tdim, acqhdr.idx.kspace_encode_step_1, acqhdr.idx.kspace_encode_step_2, slice_loc, getN(acqhdr.idx),  getS(acqhdr.idx));
                memcpy(trajptr, & acqtraj(acqhdr.discard_pre, tdim), sizeof(float)*npts_to_copy);
            }
        }
    }

    
    //Send all the ReconData messages
    //GADGET_DEBUG2("End of bucket reached, sending out %d ReconData buffers\n", recon_data_buffers.size());
    for(std::map<size_t, GadgetContainerMessage<IsmrmrdReconData>* >::iterator it = recon_data_buffers.begin(); it != recon_data_buffers.end(); it++)
    {
        //std::cout << "Sending: " << it->first << std::endl;
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

GADGET_FACTORY_DECLARE(BucketToBufferGadget)

}
