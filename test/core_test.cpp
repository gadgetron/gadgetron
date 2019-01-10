
#include <gtest/gtest.h>
#include "Message.h"
#include "Channel.h"
#include "Types.h"

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
        bool convertible = convertible_to<std::string, optional<int>>(*message);
        EXPECT_TRUE(convertible);
    }

    {
        outputChannel.push(std::make_unique<std::string>("test"), std::make_unique<std::string>("test"));

        auto message = inputChannel.pop();
        bool convertible = convertible_to<std::string, optional<int>, std::string>(*message);
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
        bool convertible = convertible_to<std::string, optional<int>>(*message);
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


TEST(TypeTests,optionaltype3){

    using namespace Gadgetron::Core;


    MessageChannel channel;
    InputChannel& inputChannel = channel;
    OutputChannel& outputChannel = channel;


    outputChannel.push(std::make_unique<std::string>("hello"));

    auto message = inputChannel.pop();

    bool convertible = convertible_to<optional<std::string>,std::string>(*message);
    EXPECT_TRUE(convertible);
}


TEST(TypeTests,varianttype){
    using namespace Gadgetron::Core;


    MessageChannel channel;
    InputChannel& inputChannel = channel;
    OutputChannel& outputChannel = channel;

    outputChannel.push(std::make_unique<std::string>("hello"));

    {
        auto message = inputChannel.pop();

        bool convertible = convertible_to<variant<std::string, int>>(*message);
        EXPECT_TRUE(convertible);
    }

    {
        outputChannel.push(std::make_unique<std::string>("hello"));

        auto message = inputChannel.pop();

        bool convertible = convertible_to<variant<int,std::string>>(*message);
        EXPECT_TRUE(convertible);
    }
}
TEST(TypeTests,varianttype2){
    using namespace Gadgetron::Core;


    MessageChannel channel;
    InputChannel& inputChannel = channel;
    OutputChannel& outputChannel = channel;

    outputChannel.push(std::make_unique<std::string>("hello"));

    {
        auto message = inputChannel.pop();

        auto variation = unpack<variant<std::string, int>>(std::move(message));
        EXPECT_FALSE(variation->empty());
        std::cout << variation->type().name() << std::endl;
        EXPECT_EQ(variation->which(),1);
    }

}

TEST(TypeTests,tupletype){
    using namespace Gadgetron::Core;


    MessageChannel channel;
    InputChannel& inputChannel = channel;
    OutputChannel& outputChannel = channel;

    outputChannel.push(std::make_unique<std::string>("hello"),std::make_unique<float>(1.0f),std::make_unique<int>(42));

    auto message = inputChannel.pop();

    bool convertible = convertible_to<tuple<std::string,float, int>>(*message);
    EXPECT_TRUE(convertible);
}


TEST(TypeTests,tuplevarianttype){
    using namespace Gadgetron::Core;


    MessageChannel channel;
    InputChannel& inputChannel = channel;
    OutputChannel& outputChannel = channel;

    outputChannel.push(std::make_unique<std::string>("hello"),std::make_unique<float>(1.0f),std::make_unique<int>(42));

    auto message = inputChannel.pop();

    bool convertible = convertible_to<variant<tuple<std::string,float, int>,float>>(*message);
    EXPECT_TRUE(convertible);
}
