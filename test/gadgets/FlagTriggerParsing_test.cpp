//
// Created by dchansen on 9/9/20.
//
#include <gtest/gtest.h>
#include "../../gadgets/mri_core/FlagTriggerGadget.h"

using namespace Gadgetron;

TEST(FlagTrigger,simple){

    auto func = Gadgetron::FlagTriggerGadget::create_trigger_filter("last_in_slice");

    auto acquisition = Core::Acquisition();
    auto& head = std::get<ISMRMRD::AcquisitionHeader>(acquisition);


    ASSERT_FALSE(func(acquisition));
    head.setFlag(ISMRMRD::ISMRMRD_ACQ_LAST_IN_SLICE);
    //bool tmp = func(acquisition);
    ASSERT_TRUE(func(acquisition));

    head.clearAllFlags();
    head.setFlag(ISMRMRD::ISMRMRD_ACQ_FIRST_IN_ENCODE_STEP1);

    ASSERT_FALSE(func(acquisition));

    head.setFlag(ISMRMRD::ISMRMRD_ACQ_LAST_IN_SLICE);

    ASSERT_TRUE(func(acquisition));
}

TEST(FlagTrigger,or_test){

    auto func = Gadgetron::FlagTriggerGadget::create_trigger_filter("last_in_slice  ||   last_in_repetition");

    auto acquisition = Core::Acquisition();
    auto& head = std::get<ISMRMRD::AcquisitionHeader>(acquisition);


    ASSERT_FALSE(func(acquisition));
    head.setFlag(ISMRMRD::ISMRMRD_ACQ_LAST_IN_SLICE);
    //bool tmp = func(acquisition);
    ASSERT_TRUE(func(acquisition));

    head.clearAllFlags();
    head.setFlag(ISMRMRD::ISMRMRD_ACQ_LAST_IN_REPETITION);

    ASSERT_TRUE(func(acquisition));

    head.setFlag(ISMRMRD::ISMRMRD_ACQ_LAST_IN_SLICE);

    ASSERT_TRUE(func(acquisition));
}

TEST(FlagTrigger,and_test){

    auto func = Gadgetron::FlagTriggerGadget::create_trigger_filter("last_in_slice  &&   last_in_repetition");

    auto acquisition = Core::Acquisition();
    auto& head = std::get<ISMRMRD::AcquisitionHeader>(acquisition);


    ASSERT_FALSE(func(acquisition));
    head.setFlag(ISMRMRD::ISMRMRD_ACQ_LAST_IN_SLICE);
    //bool tmp = func(acquisition);
    ASSERT_FALSE(func(acquisition));

    head.clearAllFlags();
    head.setFlag(ISMRMRD::ISMRMRD_ACQ_LAST_IN_REPETITION);

    ASSERT_FALSE(func(acquisition));

    head.setFlag(ISMRMRD::ISMRMRD_ACQ_LAST_IN_SLICE);

    ASSERT_TRUE(func(acquisition));
}

TEST(FlagTrigger,not_test){

auto func = Gadgetron::FlagTriggerGadget::create_trigger_filter("last_in_slice   &&  !last_in_repetition");

auto acquisition = Core::Acquisition();
auto& head = std::get<ISMRMRD::AcquisitionHeader>(acquisition);


ASSERT_FALSE(func(acquisition));
head.setFlag(ISMRMRD::ISMRMRD_ACQ_LAST_IN_SLICE);
//bool tmp = func(acquisition);
ASSERT_TRUE(func(acquisition));

head.clearAllFlags();
head.setFlag(ISMRMRD::ISMRMRD_ACQ_LAST_IN_REPETITION);

ASSERT_FALSE(func(acquisition));

head.setFlag(ISMRMRD::ISMRMRD_ACQ_LAST_IN_SLICE);

ASSERT_FALSE(func(acquisition));
}


TEST(FlagTrigger,presedence_test){

    auto func = Gadgetron::FlagTriggerGadget::create_trigger_filter("last_in_slice || last_in_repetition && first_in_slice");

    auto acquisition = Core::Acquisition();
    auto& head = std::get<ISMRMRD::AcquisitionHeader>(acquisition);


    ASSERT_FALSE(func(acquisition));
    head.setFlag(ISMRMRD::ISMRMRD_ACQ_LAST_IN_SLICE);
//bool tmp = func(acquisition);
    ASSERT_TRUE(func(acquisition));

    head.clearAllFlags();
    head.setFlag(ISMRMRD::ISMRMRD_ACQ_LAST_IN_REPETITION);

    ASSERT_FALSE(func(acquisition));

    head.setFlag(ISMRMRD::ISMRMRD_ACQ_FIRST_IN_SLICE);

    ASSERT_TRUE(func(acquisition));

    head.setFlag(ISMRMRD::ISMRMRD_ACQ_LAST_IN_SLICE);

    ASSERT_TRUE(func(acquisition));
}

TEST(FlagTrigger,parenthesis_test){

    auto func = Gadgetron::FlagTriggerGadget::create_trigger_filter("( last_in_slice || last_in_repetition ) && first_in_slice");

    auto acquisition = Core::Acquisition();
    auto& head = std::get<ISMRMRD::AcquisitionHeader>(acquisition);


    ASSERT_FALSE(func(acquisition));
    head.setFlag(ISMRMRD::ISMRMRD_ACQ_LAST_IN_SLICE);
    ASSERT_FALSE(func(acquisition));
    head.setFlag(ISMRMRD::ISMRMRD_ACQ_FIRST_IN_SLICE);
    ASSERT_TRUE(func(acquisition));

}


