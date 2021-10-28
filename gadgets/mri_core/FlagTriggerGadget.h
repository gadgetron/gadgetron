//
// Created by david on 24/08/2020.
//
#pragma once
#include "Node.h"

namespace Gadgetron {

class FlagTriggerGadget : public Core::ChannelGadget<Core::Acquisition> {
  public:
    enum class TriggerFlags : uint64_t {
        first_in_encode_step1 = 0,
        last_in_encode_step1 = 1,
        first_in_encode_step2 = 2,
        last_in_encode_step2 = 3,
        first_in_average = 4,
        last_in_average = 5,
        first_in_slice = 6,
        last_in_slice = 7,
        first_in_contrast = 8,
        last_in_contrast = 9,
        first_in_phase = 10,
        last_in_phase = 11,
        first_in_repetition = 12,
        last_in_repetition = 13,
        first_in_set = 14,
        last_in_set = 15,
        first_in_segment = 16,
        last_in_segment = 17,
        is_noise_measurement = 18,
        is_parallel_calibration = 19,
        is_parallel_calibration_and_imaging = 20,
        is_reverse = 21,
        is_navigation_data = 22,
        is_phasecorr_data = 23,
        last_in_measurement = 24,
        is_hpfeedback_data = 25,
        is_dummyscan_data = 26,
        is_rtfeedback_data = 27,
        is_surfacecoilcorrectionscan_data = 29,

        compression1 = 52,
        compression2 = 53,
        compression3 = 54,
        compression4 = 55,
        user1 = 56,
        user2 = 57,
        user3 = 58,
        user4 = 59,
        user5 = 60,
        user6 = 61,
        user7 = 62,
        user8 = 63
    };

    FlagTriggerGadget(const Core::Context& context, const Core::GadgetProperties& props);
    ~FlagTriggerGadget() override = default;

    void process(Core::InputChannel<Core::Acquisition>& input,
                 Core::OutputChannel& output) override;

    NODE_PROPERTY(trigger_flags, std::string, "Trigger flags (separated by comma)", "");

    static std::function<bool(const Core::Acquisition& acq)> create_trigger_filter(const std::string& trigger_string);

  private:
    std::function<bool(const Core::Acquisition&)> predicate;
};


} // namespace Gadgetron
