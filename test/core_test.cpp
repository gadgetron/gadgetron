
#include <gtest/gtest.h>
#include "Message.h"
#include "Channel.h"

TEST(TypeTests,multitype){
    using namespace Gadgetron::Core;

    MessageChannel channel;
    InputChannel& inputChannel = channel;
    OutputChannel& outputChannel = channel;

    outputChannel.push(std::make_unique<std::string>("test"),std::make_unique<int>(4));


    auto message = inputChannel.pop();

    bool convertible = convertible_to<std::string,int>(*message);
    EXPECT_TRUE(convertible);

}


TEST(TypeTests,singletype){
    using namespace Gadgetron::Core;

    MessageChannel channel;
    InputChannel& inputChannel = channel;
    OutputChannel& outputChannel = channel;

    outputChannel.push(std::make_unique<std::string>("test"));


    auto message = inputChannel.pop();

    bool convertible = convertible_to<std::string>(*message);
    EXPECT_TRUE(convertible);

}
TEST(TypeTests,optionaltype){
    using namespace Gadgetron::Core;

    MessageChannel channel;
    InputChannel& inputChannel = channel;
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
    InputChannel& inputChannel = channel;
    OutputChannel& outputChannel = channel;


    {
        outputChannel.push(std::make_unique<std::string>("test"));

        auto message = inputChannel.pop();
        bool convertible = convertible_to<std::string, boost::optional<int>>(*message);
        EXPECT_TRUE(convertible);
    }
}


TEST(TypeTests,converttype){

    using namespace Gadgetron::Core;


    MessageChannel channel;
    InputChannel& inputChannel = channel;
    OutputChannel& outputChannel = channel;


    outputChannel.push(std::make_unique<std::string>("hello"),std::make_unique<std::string>("world"));

    auto message = inputChannel.pop();

    auto pack = unpack<std::string,std::string>(std::move(message));

    EXPECT_EQ(*std::get<0>(pack),"hello");
}
