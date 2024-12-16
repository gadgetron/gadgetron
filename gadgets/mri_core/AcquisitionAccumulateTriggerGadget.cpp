#include "AcquisitionAccumulateTriggerGadget.h"
#include "log.h"
#include "mri_core_utility.h"
#include <boost/algorithm/string.hpp>

namespace Gadgetron {
    using TriggerDimension = AcquisitionAccumulateTriggerGadget::TriggerDimension;
    namespace {
        bool is_noise(mrd::Acquisition& acq) {
            return acq.head.flags.HasFlags(mrd::AcquisitionFlags::kIsNoiseMeasurement);
        }

        unsigned int get_index(const mrd::AcquisitionHeader& head, TriggerDimension index) {
            switch (index) {
            case TriggerDimension::kspace_encode_step_1: return head.idx.kspace_encode_step_1.value_or(0);
            case TriggerDimension::kspace_encode_step_2: return head.idx.kspace_encode_step_2.value_or(0);
            case TriggerDimension::average: return head.idx.average.value_or(0);
            case TriggerDimension::slice: return head.idx.slice.value_or(0);
            case TriggerDimension::contrast: return head.idx.contrast.value_or(0);
            case TriggerDimension::phase: return head.idx.phase.value_or(0);
            case TriggerDimension::repetition: return head.idx.repetition.value_or(0);
            case TriggerDimension::set: return head.idx.set.value_or(0);
            case TriggerDimension::segment: return head.idx.segment.value_or(0);
            case TriggerDimension::user_0: return head.idx.user[0];
            case TriggerDimension::user_1: return head.idx.user[1];
            case TriggerDimension::user_2: return head.idx.user[2];
            case TriggerDimension::user_3: return head.idx.user[3];
            case TriggerDimension::user_4: return head.idx.user[4];
            case TriggerDimension::user_5: return head.idx.user[5];
            case TriggerDimension::user_6: return head.idx.user[6];
            case TriggerDimension::user_7: return head.idx.user[7];
            case TriggerDimension::n_acquisitions: return 0;
            case TriggerDimension::none: return 0;
            }
            throw std::runtime_error("Illegal enum");
        }

        struct EqualityTrigger {
            explicit EqualityTrigger(TriggerDimension trig) : trigger{ trig } {}
            const TriggerDimension trigger;
            std::optional<unsigned int> previous_trigger;
            bool trigger_before(const mrd::Acquisition& acq) {
                auto acq_index   = get_index(acq.head, trigger);
                auto result      = (previous_trigger != acq_index && previous_trigger != std::nullopt);
                previous_trigger = acq_index;
                return result;
            }

            static bool trigger_after(const mrd::Acquisition& acq) {
                return false;
            }
        };

        struct NumAcquisitionsTrigger {
            explicit NumAcquisitionsTrigger(size_t target_acquisitions_first, size_t target_acquisitions) : target_acquisitions_first{target_acquisitions_first},target_acquisitions{ target_acquisitions } {}
            const size_t target_acquisitions_first;
            const size_t target_acquisitions;
            size_t num_acquisitions = 0;
            bool first = true;

            static bool trigger_before(const mrd::Acquisition& acq) {
                return false;
            }
            bool trigger_after(const mrd::Acquisition& acq) {
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
            static bool trigger_before(const mrd::Acquisition& acq) {
                return false;
            }
            static bool trigger_after(const mrd::Acquisition& acq) {
                return false;
            }
        };

        using Trigger = std::variant<EqualityTrigger, NumAcquisitionsTrigger, NoneTrigger>;

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
            case TriggerDimension::n_acquisitions: return NumAcquisitionsTrigger(gadget.n_acquisitions_before_trigger,gadget.n_acquisitions_before_ongoing_trigger);
            case TriggerDimension::none: return NoneTrigger();
            default: throw std::runtime_error("ENUM TriggerDimension is in an invalid state.");
            }
        }

        bool trigger_before(Trigger& trigger, const mrd::Acquisition& acq) {
            return std::visit([&](auto& var) { return var.trigger_before(acq); }, trigger);
        }
        bool trigger_after(Trigger& trigger, const mrd::Acquisition& acq) {
            return std::visit([&](auto& var) { return var.trigger_after(acq); }, trigger);
        }
    }

    void AcquisitionAccumulateTriggerGadget::send_data(Core::OutputChannel& out, std::map<unsigned int, mrd::AcquisitionBucket>& buckets,
                                                       std::vector<mrd::WaveformUint32>& waveforms)
    {
        trigger_events++;
        GDEBUG_STREAM("Trigger " << trigger_events << " occurred, sending out " << buckets.size() << " buckets, " << waveforms.size() << " waveforms ... ");
        if(!waveforms.empty()) {
            buckets.begin()->second.waveforms = std::move(waveforms);
        }
        // Pass all buckets down the chain
        for (auto& bucket : buckets) {
            out.push(std::move(bucket.second));
        }

        buckets.clear();
    }
    void AcquisitionAccumulateTriggerGadget::process(Core::InputChannel<std::variant<mrd::Acquisition, mrd::WaveformUint32>>& in, Core::OutputChannel& out)
    {
        auto waveforms = std::vector<mrd::WaveformUint32>{};
        auto buckets   = std::map<unsigned int, mrd::AcquisitionBucket>{};
        auto trigger   = get_trigger(*this);

        size_t count = 0;
        for (auto message : in) {
            if (std::holds_alternative<mrd::WaveformUint32>(message)) {
                waveforms.emplace_back(std::move(std::get<mrd::WaveformUint32>(message)));
                continue;
            }

            auto& acq = std::get<mrd::Acquisition>(message);
            if (is_noise(acq)) {
                continue;
            }

            if (trigger_before(trigger, acq)) {
                send_data(out, buckets, waveforms);
            }
            // It is enough to put the first one, since they are linked
            auto sorting_index = get_index(acq.head, sorting_dimension);

            mrd::AcquisitionBucket& bucket = buckets[sorting_index];
            Gadgetron::add_acquisition_to_bucket(bucket, std::move(acq));

            if (trigger_after(trigger, acq)) {
                send_data(out, buckets, waveforms);
            }
            count++;
        }
        GDEBUG_STREAM("AcquisitionAccumulateTriggerGadget processed " << count << " Acquisitions total");
        send_data(out, buckets, waveforms);
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
            { "none", TriggerDimension::none }, { "", TriggerDimension::none }
        };
    }

    void from_string(const std::string& str, TriggerDimension& trigger) {
        auto lower = str;
        boost::to_lower(lower);
        trigger = triggerdimension_from_name.at(lower);
    }

} // namespace Gadgetron