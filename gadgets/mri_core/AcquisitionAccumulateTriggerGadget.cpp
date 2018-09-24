#include "GadgetIsmrmrdReadWrite.h"
#include "AcquisitionAccumulateTriggerGadget.h"
#include "mri_core_data.h"
#include "log.h"

namespace Gadgetron {

    AcquisitionAccumulateTriggerGadget::~AcquisitionAccumulateTriggerGadget()
    {
        //The buckets array should be empty but just in case, let's make sure all the stuff is released.
        for (map_type_::iterator it = buckets_.begin(); it != buckets_.end(); it++) {
            if (it->second) {
                it->second->release();
            }
        }
    }

    int AcquisitionAccumulateTriggerGadget
        ::process_config(ACE_Message_Block* mb)
    {

        std::string trigger_dimension_local = trigger_dimension.value();
        std::string sorting_dimension_local = sorting_dimension.value();

        if (trigger_dimension_local.size() == 0) {
            trigger_ = NONE;
        }
        else if (trigger_dimension_local.compare("kspace_encode_step_1") == 0) {
            trigger_ = KSPACE_ENCODE_STEP_1;
        }
        else if (trigger_dimension_local.compare("kspace_encode_step_2") == 0) {
            trigger_ = KSPACE_ENCODE_STEP_2;
        }
        else if (trigger_dimension_local.compare("average") == 0) {
            trigger_ = AVERAGE;
        }
        else if (trigger_dimension_local.compare("slice") == 0) {
            trigger_ = SLICE;
        }
        else if (trigger_dimension_local.compare("contrast") == 0) {
            trigger_ = CONTRAST;
        }
        else if (trigger_dimension_local.compare("phase") == 0) {
            trigger_ = PHASE;
        }
        else if (trigger_dimension_local.compare("repetition") == 0) {
            trigger_ = REPETITION;
        }
        else if (trigger_dimension_local.compare("set") == 0) {
            trigger_ = SET;
        }
        else if (trigger_dimension_local.compare("segment") == 0) {
            trigger_ = SEGMENT;
        }
        else if (trigger_dimension_local.compare("user_0") == 0) {
            trigger_ = USER_0;
        }
        else if (trigger_dimension_local.compare("user_1") == 0) {
            trigger_ = USER_1;
        }
        else if (trigger_dimension_local.compare("user_2") == 0) {
            trigger_ = USER_2;
        }
        else if (trigger_dimension_local.compare("user_3") == 0) {
            trigger_ = USER_3;
        }
        else if (trigger_dimension_local.compare("user_4") == 0) {
            trigger_ = USER_4;
        }
        else if (trigger_dimension_local.compare("user_5") == 0) {
            trigger_ = USER_5;
        }
        else if (trigger_dimension_local.compare("user_6") == 0) {
            trigger_ = USER_6;
        }
        else if (trigger_dimension_local.compare("user_7") == 0) {
            trigger_ = USER_7;
        }
        else {
            GDEBUG("WARNING: Unknown trigger dimension (%s), trigger condition set to NONE (end of scan)", trigger_dimension_local.c_str());
            trigger_ = NONE;
        }

        GDEBUG("TRIGGER DIMENSION IS: %s (%d)\n", trigger_dimension_local.c_str(), trigger_);

        if (sorting_dimension_local.size() == 0) {
            sort_ = NONE;
        }
        else if (sorting_dimension_local.compare("kspace_encode_step_1") == 0) {
            sort_ = KSPACE_ENCODE_STEP_1;
        }
        else if (sorting_dimension_local.compare("kspace_encode_step_2") == 0) {
            sort_ = KSPACE_ENCODE_STEP_2;
        }
        else if (sorting_dimension_local.compare("average") == 0) {
            sort_ = AVERAGE;
        }
        else if (sorting_dimension_local.compare("slice") == 0) {
            sort_ = SLICE;
        }
        else if (sorting_dimension_local.compare("contrast") == 0) {
            sort_ = CONTRAST;
        }
        else if (sorting_dimension_local.compare("phase") == 0) {
            sort_ = PHASE;
        }
        else if (sorting_dimension_local.compare("repetition") == 0) {
            sort_ = REPETITION;
        }
        else if (sorting_dimension_local.compare("set") == 0) {
            sort_ = SET;
        }
        else if (sorting_dimension_local.compare("segment") == 0) {
            sort_ = SEGMENT;
        }
        else if (sorting_dimension_local.compare("user_0") == 0) {
            sort_ = USER_0;
        }
        else if (sorting_dimension_local.compare("user_1") == 0) {
            sort_ = USER_1;
        }
        else if (sorting_dimension_local.compare("user_2") == 0) {
            sort_ = USER_2;
        }
        else if (sorting_dimension_local.compare("user_3") == 0) {
            sort_ = USER_3;
        }
        else if (sorting_dimension_local.compare("user_4") == 0) {
            sort_ = USER_4;
        }
        else if (sorting_dimension_local.compare("user_5") == 0) {
            sort_ = USER_5;
        }
        else if (sorting_dimension_local.compare("user_6") == 0) {
            sort_ = USER_6;
        }
        else if (sorting_dimension_local.compare("user_7") == 0) {
            sort_ = USER_7;
        }
        else {
            GDEBUG("WARNING: Unknown sort dimension (%s), sorting set to NONE\n", sorting_dimension_local.c_str());
            sort_ = NONE;
        }

        GDEBUG("SORTING DIMENSION IS: %s (%d)\n", sorting_dimension_local.c_str(), sort_);

        trigger_events_ = 0;

        return GADGET_OK;
    }

    int AcquisitionAccumulateTriggerGadget
        ::process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1)
    {

        //Ignore noise scans
        if (m1->getObjectPtr()->isFlagSet(ISMRMRD::ISMRMRD_ACQ_IS_NOISE_MEASUREMENT)) {
            m1->release();
            return GADGET_OK;
        }

        GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2 = AsContainerMessage< hoNDArray< std::complex<float> > >(m1->cont());
        if (!m2)
        {
            GDEBUG("Error casting acquisition data package");
            return GADGET_FAIL;
        }

        //It is enough to put the first one, since they are linked
        unsigned short sorting_index = 0;
        switch (sort_) {
        case KSPACE_ENCODE_STEP_1:
            sorting_index = m1->getObjectPtr()->idx.kspace_encode_step_1;
            break;
        case KSPACE_ENCODE_STEP_2:
            sorting_index = m1->getObjectPtr()->idx.kspace_encode_step_2;
            break;
        case AVERAGE:
            sorting_index = m1->getObjectPtr()->idx.average;
            break;
        case SLICE:
            sorting_index = m1->getObjectPtr()->idx.slice;
            break;
        case CONTRAST:
            sorting_index = m1->getObjectPtr()->idx.contrast;
            break;
        case PHASE:
            sorting_index = m1->getObjectPtr()->idx.phase;
            break;
        case REPETITION:
            sorting_index = m1->getObjectPtr()->idx.repetition;
            break;
        case SET:
            sorting_index = m1->getObjectPtr()->idx.set;
            break;
        case SEGMENT:
            sorting_index = m1->getObjectPtr()->idx.segment;
            break;
        case USER_0:
            sorting_index = m1->getObjectPtr()->idx.user[0];
            break;
        case USER_1:
            sorting_index = m1->getObjectPtr()->idx.user[1];
            break;
        case USER_2:
            sorting_index = m1->getObjectPtr()->idx.user[2];
            break;
        case USER_3:
            sorting_index = m1->getObjectPtr()->idx.user[3];
            break;
        case USER_4:
            sorting_index = m1->getObjectPtr()->idx.user[4];
            break;
        case USER_5:
            sorting_index = m1->getObjectPtr()->idx.user[5];
            break;
        case USER_6:
            sorting_index = m1->getObjectPtr()->idx.user[6];
            break;
        case USER_7:
            sorting_index = m1->getObjectPtr()->idx.user[7];
            break;
        case NONE:
            sorting_index = 0;
            break;
        default:
            GDEBUG("Unknown sorting condition %d\n", sort_);
            m1->release();
            return GADGET_FAIL;
        }

        //Create the data structure that will go in the bucket
        IsmrmrdAcquisitionData d(m1, m2, AsContainerMessage< hoNDArray<float> >(m2->cont()));

        //Now let's figure out if a trigger condition has occurred.
        if (prev_.head_) { //Make sure this is not the first acquisition we are receiving
            switch (trigger_) {
            case KSPACE_ENCODE_STEP_1:
                if (prev_.head_->getObjectPtr()->idx.kspace_encode_step_1 !=
                    d.head_->getObjectPtr()->idx.kspace_encode_step_1) {
                    trigger();
                }
                break;
            case KSPACE_ENCODE_STEP_2:
                if (prev_.head_->getObjectPtr()->idx.kspace_encode_step_2 !=
                    d.head_->getObjectPtr()->idx.kspace_encode_step_2) {
                    trigger();
                }
                break;
            case AVERAGE:
                if (prev_.head_->getObjectPtr()->idx.average !=
                    d.head_->getObjectPtr()->idx.average) {
                    trigger();
                }
                break;
            case SLICE:
                if (prev_.head_->getObjectPtr()->idx.slice !=
                    d.head_->getObjectPtr()->idx.slice) {
                    trigger();
                }
                break;
            case CONTRAST:
                if (prev_.head_->getObjectPtr()->idx.contrast !=
                    d.head_->getObjectPtr()->idx.contrast) {
                    trigger();
                }
                break;
            case PHASE:
                if (prev_.head_->getObjectPtr()->idx.phase !=
                    d.head_->getObjectPtr()->idx.phase) {
                    trigger();
                }
                break;
            case REPETITION:
                if (prev_.head_->getObjectPtr()->idx.repetition !=
                    d.head_->getObjectPtr()->idx.repetition) {
                    trigger();
                }
                break;
            case SET:
                if (prev_.head_->getObjectPtr()->idx.set !=
                    d.head_->getObjectPtr()->idx.set) {
                    trigger();
                }
                break;
            case SEGMENT:
                if (prev_.head_->getObjectPtr()->idx.segment !=
                    d.head_->getObjectPtr()->idx.segment) {
                    trigger();
                }
                break;
            case USER_0:
                if (prev_.head_->getObjectPtr()->idx.user[0] !=
                    d.head_->getObjectPtr()->idx.user[0]) {
                    trigger();
                }
                break;
            case USER_1:
                if (prev_.head_->getObjectPtr()->idx.user[1] !=
                    d.head_->getObjectPtr()->idx.user[1]) {
                    trigger();
                }
                break;
            case USER_2:
                if (prev_.head_->getObjectPtr()->idx.user[2] !=
                    d.head_->getObjectPtr()->idx.user[2]) {
                    trigger();
                }
                break;
            case USER_3:
                if (prev_.head_->getObjectPtr()->idx.user[3] !=
                    d.head_->getObjectPtr()->idx.user[3]) {
                    trigger();
                }
                break;
            case USER_4:
                if (prev_.head_->getObjectPtr()->idx.user[4] !=
                    d.head_->getObjectPtr()->idx.user[4]) {
                    trigger();
                }
                break;
            case USER_5:
                if (prev_.head_->getObjectPtr()->idx.user[5] !=
                    d.head_->getObjectPtr()->idx.user[5]) {
                    trigger();
                }
                break;
            case USER_6:
                if (prev_.head_->getObjectPtr()->idx.user[6] !=
                    d.head_->getObjectPtr()->idx.user[6]) {
                    trigger();
                }
                break;
            case USER_7:
                if (prev_.head_->getObjectPtr()->idx.user[7] !=
                    d.head_->getObjectPtr()->idx.user[7]) {
                    trigger();
                }
                break;
            case NONE:
                break;
            default:
                GDEBUG("Unknown trigger condition %d\n", trigger_);
                return GADGET_FAIL;
            }
        }

        //Now we can update the previous data item that we store for 
        //purposes of determining if trigger condition has occurred. 
        prev_ = d;

        //Find the bucket the data should go in
        map_type_::iterator it = buckets_.find(sorting_index);
        if (it == buckets_.end()) {
            //Bucket does not exist, create it
            buckets_[sorting_index] = new GadgetContainerMessage<IsmrmrdAcquisitionBucket>;
        }
        IsmrmrdAcquisitionBucket* bucket = buckets_[sorting_index]->getObjectPtr();

        uint16_t espace = m1->getObjectPtr()->encoding_space_ref;

        if (!(
            ISMRMRD::FlagBit(ISMRMRD::ISMRMRD_ACQ_IS_PARALLEL_CALIBRATION).isSet(m1->getObjectPtr()->flags) ||
            ISMRMRD::FlagBit(ISMRMRD::ISMRMRD_ACQ_IS_PHASECORR_DATA).isSet(m1->getObjectPtr()->flags)
            ))
        {
            bucket->data_.push_back(d);
            if (bucket->datastats_.size() < (espace + 1)) {
                bucket->datastats_.resize(espace + 1);
            }
            bucket->datastats_[espace].kspace_encode_step_1.insert(m1->getObjectPtr()->idx.kspace_encode_step_1);
            bucket->datastats_[espace].kspace_encode_step_2.insert(m1->getObjectPtr()->idx.kspace_encode_step_2);
            bucket->datastats_[espace].slice.insert(m1->getObjectPtr()->idx.slice);
            bucket->datastats_[espace].phase.insert(m1->getObjectPtr()->idx.phase);
            bucket->datastats_[espace].contrast.insert(m1->getObjectPtr()->idx.contrast);
            bucket->datastats_[espace].set.insert(m1->getObjectPtr()->idx.set);
            bucket->datastats_[espace].segment.insert(m1->getObjectPtr()->idx.segment);
            bucket->datastats_[espace].average.insert(m1->getObjectPtr()->idx.average);
            bucket->datastats_[espace].repetition.insert(m1->getObjectPtr()->idx.repetition);
        }

        if (ISMRMRD::FlagBit(ISMRMRD::ISMRMRD_ACQ_IS_PARALLEL_CALIBRATION).isSet(m1->getObjectPtr()->flags) ||
            ISMRMRD::FlagBit(ISMRMRD::ISMRMRD_ACQ_IS_PARALLEL_CALIBRATION_AND_IMAGING).isSet(m1->getObjectPtr()->flags))
        {
            bucket->ref_.push_back(d);
            if (bucket->refstats_.size() < (espace + 1)) {
                bucket->refstats_.resize(espace + 1);
            }
            bucket->refstats_[espace].kspace_encode_step_1.insert(m1->getObjectPtr()->idx.kspace_encode_step_1);
            bucket->refstats_[espace].kspace_encode_step_2.insert(m1->getObjectPtr()->idx.kspace_encode_step_2);
            bucket->refstats_[espace].slice.insert(m1->getObjectPtr()->idx.slice);
            bucket->refstats_[espace].phase.insert(m1->getObjectPtr()->idx.phase);
            bucket->refstats_[espace].contrast.insert(m1->getObjectPtr()->idx.contrast);
            bucket->refstats_[espace].set.insert(m1->getObjectPtr()->idx.set);
            bucket->refstats_[espace].segment.insert(m1->getObjectPtr()->idx.segment);
            bucket->refstats_[espace].average.insert(m1->getObjectPtr()->idx.average);
            bucket->refstats_[espace].repetition.insert(m1->getObjectPtr()->idx.repetition);
        }

        //We can release the data now. It is reference counted and counter have been incremented through operations above. 
        m1->release();

        //TODO: 
        // At this point it would make sense to check the data flags for trigger conditions. 

        return GADGET_OK;
    }

    int AcquisitionAccumulateTriggerGadget::process(GadgetContainerMessage<ISMRMRD::ISMRMRD_WaveformHeader>* m1)
    {
        GadgetContainerMessage< hoNDArray< uint32_t > >* m2 = AsContainerMessage< hoNDArray< uint32_t > >(m1->cont());
        if (!m2)
        {
            GDEBUG("Error casting waveform data package");
            return GADGET_FAIL;
        }

        ISMRMRD::Waveform ismrmrd_wav;
        ismrmrd_wav.head = *m1->getObjectPtr();
        ismrmrd_wav.data = m2->getObjectPtr()->begin();

        // buffer the waveform data
        wav_buf_.push_back(ismrmrd_wav);

        ismrmrd_wav.data = NULL;
        m1->release();

        return GADGET_OK;
    }

    int AcquisitionAccumulateTriggerGadget::trigger()
    {
        //We will keep track of the triggers we encounter
        trigger_events_++;

        GDEBUG("Trigger (%d) occurred, sending out %d buckets\n", trigger_events_, buckets_.size());
        //Pass all buckets down the chain
        for (const auto& it : buckets_)
        {
            if (it.second)
            {
                // attach the waveform
                if (!this->wav_buf_.empty())
                {
                    it.second->getObjectPtr()->waveform_ = this->wav_buf_;
                    this->wav_buf_.clear();
                }

                if (this->next()->putq(it.second) == -1)
                {
                    it.second->release();
                    GDEBUG("Failed to pass bucket down the chain\n");
                    return GADGET_FAIL;
                }
            }
        }

        buckets_.clear();
        prev_ = IsmrmrdAcquisitionData(); //Reset previous so that we don't end up triggering again
        return GADGET_OK;
    }

    int AcquisitionAccumulateTriggerGadget::close(unsigned long flags)
    {

        int ret = Gadget::close(flags);

        if (flags != 0) {
            GDEBUG("AcquisitionAccumulateTriggerGadget::close\n");
            trigger();
        }
        return ret;
    }

    GADGET_FACTORY_DECLARE(AcquisitionAccumulateTriggerGadget)
}
