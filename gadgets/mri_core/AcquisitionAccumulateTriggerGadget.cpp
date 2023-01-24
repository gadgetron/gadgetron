#include "AcquisitionAccumulateTriggerGadget.h"
#include "log.h"
#include "mri_core_data.h"
#include <boost/algorithm/string.hpp>

namespace Gadgetron {
    using TriggerDimension = AcquisitionAccumulateTriggerGadget::TriggerDimension;
    namespace {
        bool is_noise(Core::Acquisition& acq) {
            return std::get<ISMRMRD::AcquisitionHeader>(acq).isFlagSet(ISMRMRD::ISMRMRD_ACQ_IS_NOISE_MEASUREMENT);
        }

        unsigned short get_index(const ISMRMRD::AcquisitionHeader& header, TriggerDimension index) {
            switch (index) {
            case TriggerDimension::kspace_encode_step_1: return header.idx.kspace_encode_step_1;
            case TriggerDimension::kspace_encode_step_2: return header.idx.kspace_encode_step_2;
            case TriggerDimension::average: return header.idx.average;
            case TriggerDimension::slice: return header.idx.slice;
            case TriggerDimension::contrast: return header.idx.contrast;
            case TriggerDimension::phase: return header.idx.phase;
            case TriggerDimension::repetition: return header.idx.repetition;
            case TriggerDimension::set: return header.idx.set;
            case TriggerDimension::segment: return header.idx.segment;
            case TriggerDimension::user_0: return header.idx.user[0];
            case TriggerDimension::user_1: return header.idx.user[1];
            case TriggerDimension::user_2: return header.idx.user[2];
            case TriggerDimension::user_3: return header.idx.user[3];
            case TriggerDimension::user_4: return header.idx.user[4];
            case TriggerDimension::user_5: return header.idx.user[5];
            case TriggerDimension::user_6: return header.idx.user[6];
            case TriggerDimension::user_7: return header.idx.user[7];
            case TriggerDimension::parallel_calibration: return header.isFlagSet(ISMRMRD::ISMRMRD_ACQ_IS_PARALLEL_CALIBRATION);
            case TriggerDimension::n_acquisitions: return 0;
            case TriggerDimension::none: return 0;
            }
            throw std::runtime_error("Illegal enum");
        }





        struct EqualityTrigger {
            explicit EqualityTrigger(TriggerDimension trig) : trigger{ trig } {}
            const TriggerDimension trigger;
            Core::optional<unsigned short> previous_trigger;
            bool trigger_before(const ISMRMRD::AcquisitionHeader& head) {
                auto acq_index   = get_index(head, trigger);
                auto result      = (previous_trigger != acq_index && previous_trigger != Core::none);
                previous_trigger = acq_index;
                return result;
            }

            static bool trigger_after(const ISMRMRD::AcquisitionHeader& head) {
                return false;
            }
        };
        struct FlagRemovedTrigger {
            explicit FlagRemovedTrigger(TriggerDimension trig) : trigger{trig} {}
            const TriggerDimension trigger;
            bool previous_trigger = false;
            bool trigger_before(const ISMRMRD::AcquisitionHeader& head) {
                bool flag_active = (bool) get_index(head, trigger);
                bool result = (previous_trigger && !flag_active);
                // GDEBUG("cnt: %d, Previous Trigger: %d, Trigger: %d, result: %d\n", head.scan_counter, previous_trigger, flag_active, result);
                previous_trigger = flag_active;
                return result;
            }

            static bool trigger_after(const ISMRMRD::AcquisitionHeader& head) {
                return false;
            }
        };
        struct NumAcquisitionsTrigger {
            explicit NumAcquisitionsTrigger(size_t target_acquisitions_first, size_t target_acquisitions) : target_acquisitions_first{target_acquisitions_first},target_acquisitions{ target_acquisitions } {}
            const size_t target_acquisitions_first;
            const size_t target_acquisitions;
            size_t num_acquisitions = 0;
            bool first = true;

            static bool trigger_before(const ISMRMRD::AcquisitionHeader& head) {
                return false;
            }
            bool trigger_after(const ISMRMRD::AcquisitionHeader& head) {
                // Handle possible n_acq trigger _after_ pushing data - all others come before
                size_t current_target_acquisitions = first ? target_acquisitions_first : target_acquisitions;
                bool result = ++num_acquisitions >= current_target_acquisitions;
                if (result)
                {
                    first = false;
                    num_acquisitions = 0;
                }
                return result;
            }
        };
        struct NoneTrigger {
            static bool trigger_before(const ISMRMRD::AcquisitionHeader& head) {
                return false;
            }
            static bool trigger_after(const ISMRMRD::AcquisitionHeader& head) {
                return false;
            }
        };

        using Trigger = Core::variant<EqualityTrigger, FlagRemovedTrigger, NumAcquisitionsTrigger, NoneTrigger>;

        Trigger get_trigger(const AcquisitionAccumulateTriggerGadget& gadget) {
            switch (gadget.trigger_dimension) {

            case TriggerDimension::kspace_encode_step_1:
            case TriggerDimension::kspace_encode_step_2:
            case TriggerDimension::average:
            case TriggerDimension::slice:
            case TriggerDimension::contrast:
            case TriggerDimension::phase:
            case TriggerDimension::repetition:
            case TriggerDimension::set:
            case TriggerDimension::segment:
            case TriggerDimension::user_0:
            case TriggerDimension::user_1:
            case TriggerDimension::user_2:
            case TriggerDimension::user_3:
            case TriggerDimension::user_4:
            case TriggerDimension::user_5:
            case TriggerDimension::user_6:
            case TriggerDimension::user_7: return EqualityTrigger(gadget.trigger_dimension);
            case TriggerDimension::parallel_calibration: return FlagRemovedTrigger(gadget.trigger_dimension);
            case TriggerDimension::n_acquisitions: return NumAcquisitionsTrigger(gadget.n_acquisitions_before_trigger,gadget.n_acquisitions_before_ongoing_trigger);
            case TriggerDimension::none: return NoneTrigger();
            default: throw std::runtime_error("ENUM TriggerDimension is in an invalid state.");
            }
        }

        bool trigger_before(Trigger& trigger, const ISMRMRD::AcquisitionHeader& head) {
            return Core::visit([&](auto& var) { return var.trigger_before(head); }, trigger);
        }
        bool trigger_after(Trigger& trigger, const ISMRMRD::AcquisitionHeader& acq) {
            return Core::visit([&](auto& var) { return var.trigger_after(acq); }, trigger);
        }


    }

    void AcquisitionAccumulateTriggerGadget::send_data(Core::OutputChannel& out, std::map<unsigned short, AcquisitionBucket>& buckets,
                                                       std::vector<Core::Waveform>& waveforms) {
        trigger_events++;
        GDEBUG("Trigger (%d) occurred, sending out %d buckets\n", trigger_events, buckets.size());
        buckets.begin()->second.waveform_ = std::move(waveforms);
        // Pass all buckets down the chain
        for (auto& bucket : buckets)
            out.push(std::move(bucket.second));

        buckets.clear();
    }
    void AcquisitionAccumulateTriggerGadget ::process(
        Core::InputChannel<Core::variant<Core::Acquisition, Core::Waveform>>& in, Core::OutputChannel& out) {

        auto waveforms = std::vector<Core::Waveform>{};
        auto buckets   = std::map<unsigned short, AcquisitionBucket>{};
        auto trigger   = get_trigger(*this);

        for (auto message : in) {
            if (Core::holds_alternative<Core::Waveform>(message)) {
                waveforms.emplace_back(std::move(Core::get<Core::Waveform>(message)));
                continue;
            }

            auto& acq = Core::get<Core::Acquisition>(message);
            if (is_noise(acq))
                continue;
            auto head = std::get<ISMRMRD::AcquisitionHeader>(acq);

            if (trigger_before(trigger, head))
                send_data(out, buckets, waveforms);
            // It is enough to put the first one, since they are linked
            unsigned short sorting_index = get_index(head, sorting_dimension);

            AcquisitionBucket& bucket = buckets[sorting_index];
            bucket.add_acquisition(std::move(acq));

            if (trigger_after(trigger, head))
                send_data(out, buckets, waveforms);
        }
        send_data(out,buckets,waveforms);
    }
    GADGETRON_GADGET_EXPORT(AcquisitionAccumulateTriggerGadget);

    namespace {
        const std::map<std::string, TriggerDimension> triggerdimension_from_name = {

            { "kspace_encode_step_1", TriggerDimension::kspace_encode_step_1 },
            { "kspace_encode_step_2", TriggerDimension::kspace_encode_step_2 },
            { "average", TriggerDimension::average }, { "slice", TriggerDimension::slice },
            { "contrast", TriggerDimension::contrast }, { "phase", TriggerDimension::phase },
            { "repetition", TriggerDimension::repetition }, { "set", TriggerDimension::set },
            { "segment", TriggerDimension::segment }, { "user_0", TriggerDimension::user_0 },
            { "user_1", TriggerDimension::user_1 }, { "user_2", TriggerDimension::user_2 },
            { "user_3", TriggerDimension::user_3 }, { "user_4", TriggerDimension::user_4 },
            { "user_5", TriggerDimension::user_5 }, { "user_6", TriggerDimension::user_6 },
            { "user_7", TriggerDimension::user_7 }, { "n_acquisitions", TriggerDimension::n_acquisitions },
            { "parallel_calibration", TriggerDimension::parallel_calibration },
            { "none", TriggerDimension::none }, { "", TriggerDimension::none }
        };
    }

    void from_string(const std::string& str, TriggerDimension& trigger) {
        auto lower = str;
        boost::to_lower(lower);
        trigger = triggerdimension_from_name.at(lower);
    }

}
