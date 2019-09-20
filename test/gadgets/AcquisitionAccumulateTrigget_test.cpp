//
// Created by dchansen on 9/20/19.
//
#include "../../gadgets/mri_core/AcquisitionAccumulateTriggerGadget.h"
#include "setup_gadget.h"
#include <future>
#include <gtest/gtest.h>
using namespace Gadgetron;
using namespace Gadgetron::Test;
using namespace std::string_literals;
using namespace std::chrono_literals;
//class AcquisitionAccumulateTriggerTest : public ::testing::Test {
//public:
//    AcquisitionAccumulateTriggerTest()
//        : channels{ setup_gadget<AcquisitionAccumulateTriggerGadget>({ { "trigger_dimension"s, "slice"s } }) } {}
//    GadgetChannels<AcquisitionAccumulateTriggerGadget> channels;
//};

TEST(AcquisitionAccumulateTriggerTest, slice_trigger) {

    try {
        auto channels = setup_gadget<AcquisitionAccumulateTriggerGadget>({ { "trigger_dimension"s, "slice"s } });

        for (size_t i = 0; i < 11; i++) {
            auto acq                      = generate_acquisition(192, 16);
            auto& head                    = std::get<ISMRMRD::AcquisitionHeader>(acq);
            head.idx.kspace_encode_step_1 = i;
            channels.input.push(acq);
        }

        auto acq   = generate_acquisition(192, 162);
        auto& head = std::get<ISMRMRD::AcquisitionHeader>(acq);
        head.idx.slice++;
        channels.input.push(acq);

        auto message_future = std::async([&]() { return channels.output.pop(); });

        auto ec = message_future.wait_for(1000ms);
        ASSERT_EQ(ec, std::future_status::ready);
        auto message = message_future.get();

        ASSERT_TRUE(Core::convertible_to<AcquisitionBucket>(message));
        auto bucket = Core::force_unpack<AcquisitionBucket>(std::move(message));
        ASSERT_EQ(bucket.data_.size(), 11);
    } catch (const Core::ChannelClosed&){}

}