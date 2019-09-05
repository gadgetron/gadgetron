
#include <gtest/gtest.h>
#include "Message.h"
#include "Channel.h"
#include "Types.h"

TEST(TypeTests, multitype) {
    using namespace Gadgetron::Core;

    auto channel = make_channel<MessageChannel>();
    GenericInputChannel inputChannel = std::move(channel.input);
    OutputChannel outputChannel = std::move(channel.output);

    outputChannel.push(std::string("test"), int(4));


    auto message = inputChannel.pop();

    bool convertible = convertible_to<std::string, int>(message);
    EXPECT_TRUE(convertible);

}


TEST(TypeTests, singletype) {
    using namespace Gadgetron::Core;
    auto channel = make_channel<MessageChannel>();
    GenericInputChannel inputChannel = std::move(channel.input);
    OutputChannel outputChannel = std::move(channel.output);


    outputChannel.push(std::string("test"));


    auto message = inputChannel.pop();

    bool convertible = convertible_to<std::string>(message);
    EXPECT_TRUE(convertible);

}

TEST(TypeTests, optionaltype) {
    using namespace Gadgetron::Core;

    auto channel = make_channel<MessageChannel>();
    GenericInputChannel inputChannel = std::move(channel.input);
    OutputChannel outputChannel = std::move(channel.output);

    outputChannel.push(std::string("test"), int(1));
    {
        auto message = inputChannel.pop();
        bool convertible = convertible_to<std::string, optional<int>>(message);
        EXPECT_TRUE(convertible);
    }

    {
        outputChannel.push(std::string("test"), std::string("test"));

        auto message = inputChannel.pop();
        bool convertible = convertible_to<std::string, optional<int>, std::string>(message);
        EXPECT_TRUE(convertible);
    }

}

TEST(TypeTests, optionaltype2) {
    using namespace Gadgetron::Core;

    auto channel = make_channel<MessageChannel>();
    GenericInputChannel inputChannel = std::move(channel.input);
    OutputChannel outputChannel = std::move(channel.output);


    {
        outputChannel.push(std::string("test"));

        auto message = inputChannel.pop();
        bool convertible = convertible_to<std::string, optional<int>>(message);
        EXPECT_TRUE(convertible);
    }
}


TEST(TypeTests, converttype) {

    using namespace Gadgetron::Core;


    auto channel = make_channel<MessageChannel>();
    GenericInputChannel inputChannel = std::move(channel.input);
    OutputChannel outputChannel = std::move(channel.output);


    outputChannel.push(std::string("hello"), std::string("world"));

    auto message = inputChannel.pop();

    auto pack = unpack<std::string, std::string>(std::move(message));

    EXPECT_EQ(std::get<0>(*pack), "hello");
}


TEST(TypeTests, optionaltype3) {

    using namespace Gadgetron::Core;

    auto channel = make_channel<MessageChannel>();
    GenericInputChannel inputChannel = std::move(channel.input);
    OutputChannel outputChannel = std::move(channel.output);


    outputChannel.push(std::string("hello"));

    auto message = inputChannel.pop();

    bool convertible = convertible_to<optional<std::string>, std::string>(message);
    EXPECT_TRUE(convertible);
}


TEST(TypeTests, varianttype) {
    using namespace Gadgetron::Core;
    auto channel = make_channel<MessageChannel>();
    GenericInputChannel inputChannel = std::move(channel.input);
    OutputChannel outputChannel = std::move(channel.output);

    outputChannel.push(std::string("hello"));

    {
        auto message = inputChannel.pop();

        bool convertible = convertible_to<variant<std::string, int>>(message);
        EXPECT_TRUE(convertible);
    }

    {
        outputChannel.push(std::string("hello"));

        auto message = inputChannel.pop();

        bool convertible = convertible_to<variant<int, std::string>>(message);
        EXPECT_TRUE(convertible);
    }
}

TEST(TypeTests, varianttype2) {
    using namespace Gadgetron::Core;

    auto channel = make_channel<MessageChannel>();
    GenericInputChannel inputChannel = std::move(channel.input);
    OutputChannel outputChannel = std::move(channel.output);

    outputChannel.push(std::string("hello"));

    {
        auto message = inputChannel.pop();

        auto variation = force_unpack<variant<std::string, int>>(std::move(message));
        EXPECT_EQ(variation.index(), 0);
    }

}

TEST(TypeTests, tupletype) {
    using namespace Gadgetron::Core;


    auto channel = make_channel<MessageChannel>();
    GenericInputChannel inputChannel = std::move(channel.input);
    OutputChannel outputChannel = std::move(channel.output);

    outputChannel.push(std::string("hello"), 1.0f, int(42));

    auto message = inputChannel.pop();

    bool convertible = convertible_to<tuple<std::string, float, int>>(message);
    EXPECT_TRUE(convertible);
}


TEST(TypeTests, tuplevarianttype) {
    using namespace Gadgetron::Core;

    auto channel = make_channel<MessageChannel>();
    GenericInputChannel inputChannel = std::move(channel.input);
    OutputChannel outputChannel = std::move(channel.output);

    outputChannel.push(std::string("hello"), float(1.0f), int(42));

    auto message = inputChannel.pop();

    bool convertible = convertible_to<variant<tuple<std::string, float, int>, float>>(message);
    EXPECT_TRUE(convertible);
}

TEST(TypeTests, tupletype2) {
    using namespace Gadgetron::Core;

    auto channel = make_channel<MessageChannel>();
    GenericInputChannel inputChannel = std::move(channel.input);
    OutputChannel outputChannel = std::move(channel.output);

    outputChannel.push(std::make_tuple(std::string("hello"), 1.0f, int(42)));

    auto message = inputChannel.pop();

    bool convertible = convertible_to<tuple<std::string, float, int>>(message);
    EXPECT_TRUE(convertible);
}


TEST(TypeTests, tupletype3) {
    using namespace Gadgetron::Core;
    auto channel = make_channel<MessageChannel>();
    GenericInputChannel inputChannel = std::move(channel.input);
    OutputChannel outputChannel = std::move(channel.output);


    outputChannel.push(std::make_tuple(std::string("hello"), 1.0f, int(42)));

    auto message = inputChannel.pop();

    auto converted = force_unpack<tuple<std::string, float, int>>(std::move(message));

    EXPECT_EQ(std::get<1>(converted), 1.0f);
}


