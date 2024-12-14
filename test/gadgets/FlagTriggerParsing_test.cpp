//
// Created by dchansen on 9/9/20.
//
#include <gtest/gtest.h>
#include "../../gadgets/mri_core/FlagTriggerGadget.h"

using namespace Gadgetron;

TEST(FlagTrigger,simple){

    auto func = Gadgetron::FlagTriggerGadget::create_trigger_filter("last_in_slice");

    auto acquisition = mrd::Acquisition();
    auto& head = acquisition.head;

    ASSERT_FALSE(func(acquisition));
    head.flags.SetFlags(mrd::AcquisitionFlags::kLastInSlice);
    //bool tmp = func(acquisition);
    ASSERT_TRUE(func(acquisition));

    head.flags.Clear();
    head.flags.SetFlags(mrd::AcquisitionFlags::kFirstInEncodeStep1);
    ASSERT_FALSE(func(acquisition));

    head.flags.SetFlags(mrd::AcquisitionFlags::kLastInSlice);
    ASSERT_TRUE(func(acquisition));
}

TEST(FlagTrigger,or_test){

    auto func = Gadgetron::FlagTriggerGadget::create_trigger_filter("last_in_slice  ||   last_in_repetition");

    auto acquisition = mrd::Acquisition();
    auto& head = acquisition.head;

    ASSERT_FALSE(func(acquisition));
    head.flags.SetFlags(mrd::AcquisitionFlags::kLastInSlice);
    //bool tmp = func(acquisition);
    ASSERT_TRUE(func(acquisition));

    head.flags.Clear();
    head.flags.SetFlags(mrd::AcquisitionFlags::kLastInRepetition);
    ASSERT_TRUE(func(acquisition));

    head.flags.SetFlags(mrd::AcquisitionFlags::kLastInSlice);
    ASSERT_TRUE(func(acquisition));
}

TEST(FlagTrigger,and_test){

    auto func = Gadgetron::FlagTriggerGadget::create_trigger_filter("last_in_slice  &&   last_in_repetition");

    auto acquisition = mrd::Acquisition();
    auto& head = acquisition.head;

    ASSERT_FALSE(func(acquisition));
    head.flags.SetFlags(mrd::AcquisitionFlags::kLastInSlice);
    //bool tmp = func(acquisition);
    ASSERT_FALSE(func(acquisition));

    head.flags.Clear();
    head.flags.SetFlags(mrd::AcquisitionFlags::kLastInRepetition);
    ASSERT_FALSE(func(acquisition));

    head.flags.SetFlags(mrd::AcquisitionFlags::kLastInSlice);
    ASSERT_TRUE(func(acquisition));
}

TEST(FlagTrigger,not_test){
    auto func = Gadgetron::FlagTriggerGadget::create_trigger_filter("last_in_slice   &&  !last_in_repetition");

    auto acquisition = mrd::Acquisition();
    auto& head = acquisition.head;

    ASSERT_FALSE(func(acquisition));
    head.flags.SetFlags(mrd::AcquisitionFlags::kLastInSlice);
    //bool tmp = func(acquisition);
    ASSERT_TRUE(func(acquisition));

    head.flags.Clear();
    head.flags.SetFlags(mrd::AcquisitionFlags::kLastInRepetition);
    ASSERT_FALSE(func(acquisition));

    head.flags.SetFlags(mrd::AcquisitionFlags::kLastInSlice);
    ASSERT_FALSE(func(acquisition));
}


TEST(FlagTrigger,presedence_test){

    auto func = Gadgetron::FlagTriggerGadget::create_trigger_filter("last_in_slice || last_in_repetition && first_in_slice");

    auto acquisition = mrd::Acquisition();
    auto& head = acquisition.head;

    ASSERT_FALSE(func(acquisition));
    head.flags.SetFlags(mrd::AcquisitionFlags::kLastInSlice);
    //bool tmp = func(acquisition);
    ASSERT_TRUE(func(acquisition));

    head.flags.Clear();
    head.flags.SetFlags(mrd::AcquisitionFlags::kLastInRepetition);
    ASSERT_FALSE(func(acquisition));

    head.flags.SetFlags(mrd::AcquisitionFlags::kFirstInSlice);
    ASSERT_TRUE(func(acquisition));

    head.flags.SetFlags(mrd::AcquisitionFlags::kLastInSlice);
    ASSERT_TRUE(func(acquisition));
}

TEST(FlagTrigger,parenthesis_test){

    auto func = Gadgetron::FlagTriggerGadget::create_trigger_filter("( last_in_slice || last_in_repetition ) && first_in_slice");

    auto acquisition = mrd::Acquisition();
    auto& head = acquisition.head;

    ASSERT_FALSE(func(acquisition));
    head.flags.SetFlags(mrd::AcquisitionFlags::kLastInSlice);
    ASSERT_FALSE(func(acquisition));
    head.flags.SetFlags(mrd::AcquisitionFlags::kFirstInSlice);
    ASSERT_TRUE(func(acquisition));
}
