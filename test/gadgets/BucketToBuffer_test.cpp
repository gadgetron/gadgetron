//
// Created by dchansen on 9/20/19.
//
#include "../../gadgets/mri_core/BucketToBufferGadget.h"
#include "setup_gadget.h"
#include <future>
#include <gtest/gtest.h>
using namespace Gadgetron;
using namespace Gadgetron::Test;
using namespace std::string_literals;
using namespace std::chrono_literals;
// class AcquisitionAccumulateTriggerTest : public ::testing::Test {
// public:
//    AcquisitionAccumulateTriggerTest()
//        : channels{ setup_gadget<AcquisitionAccumulateTriggerGadget>({ { "trigger_dimension"s, "slice"s } }) } {}
//    GadgetChannels<AcquisitionAccumulateTriggerGadget> channels;
//};

AcquisitionBucket generate_bucket(size_t slice, size_t set) {
    auto bucket = AcquisitionBucket{};
    for (size_t i = 0; i < 11; i++) {
        auto acq                      = generate_acquisition(192, 16);
        auto& head                    = std::get<ISMRMRD::AcquisitionHeader>(acq);
        head.idx.kspace_encode_step_1 = i;
        head.idx.slice = slice;
        head.idx.set = set;
        bucket.data_.push_back(std::move(acq));
    }
    return bucket;



}

TEST(BucketToBufferGadgetTest, simpletest) {

    try {
        auto channels
            = setup_gadget<BucketToBufferGadget>({ { "N_dimension"s, "slice"s }, { "S_dimensions"s, "set" } });

        for (size_t i = 0; i < 3; i++) {
            channels.input.push(generate_bucket(i,0));

        }

        auto message_future = std::async([&]() { return channels.output.pop(); });

        auto ec = message_future.wait_for(1000ms);
        ASSERT_EQ(ec, std::future_status::ready);
        auto message = message_future.get();

        ASSERT_TRUE(Core::convertible_to<IsmrmrdReconData>(message));
        auto buffer = Core::force_unpack<IsmrmrdReconData>(std::move(message));
        ASSERT_EQ(buffer.rbit_.size(), 3);
    } catch (const Core::ChannelClosed&) {
    }
}
