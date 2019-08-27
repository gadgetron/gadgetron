//
// Created by dchansen on 4/24/19.
//
#include "ChannelAlgorithms.h"

#include <gtest/gtest.h>
#include <string>

using testing::Types;
using namespace Gadgetron::Core;
using namespace Gadgetron::Core::Algorithm;
using namespace std::string_literals;

class ChannelAlgortihmsTest : public ::testing::Test {
protected:
    void SetUp() override {

        {
            auto channels = make_channel<MessageChannel>();
            in            = std::make_shared<GenericInputChannel>(std::move(channels.input));
            out           = std::make_shared<OutputChannel>(std::move(channels.output));
            tin           = std::make_shared<InputChannel<std::string>>(*in, *out);
        }

        out->push("Penguins"s);
        out->push("Are"s);
        out->push("Awesome"s);
        out->push("Cats"s);
        out->push("Are"s);
        out->push("Weird"s);
    }

    std::shared_ptr<InputChannel<std::string>> tin;
    std::shared_ptr<GenericInputChannel> in;
    std::shared_ptr<OutputChannel> out;
};

TEST_F(ChannelAlgortihmsTest, takeWhile) {

    auto take_while_channel = take_while(*this->tin, [](std::string& message) {
        return message != "Awesome";
    });

    int i = 0;
    for (auto message : take_while_channel) {
        i++;
    }

    ASSERT_EQ(i, 3);
}
TEST_F(ChannelAlgortihmsTest, takeUntil) {

    auto take_while_channel = take_until(*this->tin, [](std::string& message) {
        return message == "Awesome";
    });

    int i = 0;
    for (auto message : take_while_channel) {
        i++;
    }

    ASSERT_EQ(i, 3);
}
