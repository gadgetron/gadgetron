
#include <gtest/gtest.h>
#include "Message.h"
#include "Channel.h"

TEST(TypeTests,multitype){
    using namespace Gadgetron::Core;

    MessageChannel channel;
    InputChannel<Message>& inputChannel = channel;
    OutputChannel& outputChannel = channel;

    outputChannel.push(std::make_unique<std::string>("test"),std::make_unique<int>(4));


    auto message = inputChannel.pop();

    bool convertible = convertible_to<std::string,int>(*message);
    EXPECT_TRUE(convertible);

}


TEST(TypeTests,singletype){
    using namespace Gadgetron::Core;

    MessageChannel channel;
    InputChannel<Message>& inputChannel = channel;
    OutputChannel& outputChannel = channel;

    outputChannel.push(std::make_unique<std::string>("test"));


    auto message = inputChannel.pop();

    bool convertible = convertible_to<std::string>(*message);
    EXPECT_TRUE(convertible);

}
TEST(TypeTests,optionaltype){
    using namespace Gadgetron::Core;

    MessageChannel channel;
    InputChannel<Message>& inputChannel = channel;
    OutputChannel& outputChannel = channel;

    outputChannel.push(std::make_unique<std::string>("test"), std::make_unique<int>(1));
    {
        auto message = inputChannel.pop();
        bool convertible = convertible_to<std::string, boost::optional<int>>(*message);
        EXPECT_TRUE(convertible);
    }

    {
        outputChannel.push(std::make_unique<std::string>("test"), std::make_unique<std::string>("test"));

        auto message = inputChannel.pop();
        bool convertible = convertible_to<std::string, boost::optional<int>, std::string>(*message);
        EXPECT_TRUE(convertible);
    }

}

TEST(TypeTests,optionaltype2){
    using namespace Gadgetron::Core;

    MessageChannel channel;
    InputChannel<Message>& inputChannel = channel;
    OutputChannel& outputChannel = channel;


    {
        outputChannel.push(std::make_unique<std::string>("test"));

        auto message = inputChannel.pop();
        bool convertible = convertible_to<std::string, boost::optional<int>>(*message);
        EXPECT_TRUE(convertible);
    }
}
