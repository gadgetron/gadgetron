#include <chrono>
#include <random>
#include <thread>

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <date/date.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <range/v3/range.hpp>

#include "storage.h"

using namespace Gadgetron::Storage;
using namespace Gadgetron;

using namespace boost::filesystem;
using namespace boost::program_options;

class StorageTest : public ::testing::Test {
  protected:
    static void SetUpTestSuite() {

        temp_dir = boost::filesystem::temp_directory_path() / "gadgetron_session_test";
        boost::filesystem::remove_all(temp_dir);
        boost::filesystem::create_directory(temp_dir);

        options_description desc("Storage options");
        // clang-format off
        desc.add_options()
        ("storage_port", value<unsigned short>()->default_value(26589))
        ("database_dir,D", value<path>()->default_value(temp_dir))
        ("storage_dir,S", value<path>()->default_value(temp_dir));
        // clang-format on

        variables_map args;
        store(parse_environment(desc, "GADGETRON_STORAGE_TEST_"), args);
        notify(args);

        auto [address, process] = Server::ensure_storage_server(args);

        ISMRMRD::IsmrmrdHeader header;
        header.subjectInformation = ISMRMRD::SubjectInformation{{}, {}, std::string("Penny the Pirate"), {}, {}};
        header.studyInformation = ISMRMRD::StudyInformation{{}, {}, std::string("YAAARH")};

        storage_address = address;
        server = std::move(*process);
        storage = Gadgetron::Server::setup_storage_spaces(address, header);
    }

    static void TearDownTestSuite() {
        boost::filesystem::remove_all(temp_dir);
        server.terminate();
        server.wait();
    }

    inline static std::string storage_address;
    inline static StorageSpaces storage;
    inline static boost::process::child server;
    inline static boost::filesystem::path temp_dir;
};

std::vector<char> generate_random_vector(size_t size) {
    std::vector<char> data(size, 0);
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<> distrib(0, 255);
    for (size_t n = 0; n < data.size(); n++) {
        data[n] = (char)distrib(gen);
    }
    return data;
}

TEST_F(StorageTest, storage_client_should_return_inserted_data_and_all_meta_tags) {
    StorageClient storage_client(storage_address);
    auto data = generate_random_vector(256);
    auto datastream = std::stringstream(std::string(data.begin(), data.end()));
    std::string session_id = "mysession3";
    auto tags = StorageItemTags::Builder("mypatient")
                    .with_device("mydevice")
                    .with_session(session_id)
                    .with_name("myvariable1")
                    .with_custom_tag("tagone", "first")
                    .with_custom_tag("tagtwo", "second")
                    .with_custom_tag("multival", "value1")
                    .with_custom_tag("multival", "value2")
                    .build();

    auto resp = storage_client.store_item(tags, datastream);
    auto out = storage_client.get_latest_item(tags);
    std::vector<char> out_bytes;
    std::copy(std::istreambuf_iterator<char>(*out), std::istreambuf_iterator<char>(), std::back_inserter(out_bytes));
    EXPECT_EQ(data, out_bytes);
    auto outByUrl = storage_client.get_item_by_url(resp.data);
    std::vector<char> outBytesByUrl;
    std::copy(std::istreambuf_iterator<char>(*outByUrl), std::istreambuf_iterator<char>(),
              std::back_inserter(outBytesByUrl));
    EXPECT_EQ(data, outBytesByUrl);

    // Get the list of items
    auto list = storage_client.list_items(
        StorageItemTags::Builder("mypatient").with_session(session_id).with_name("myvariable1").build());
    EXPECT_TRUE(list.complete);
    EXPECT_EQ(list.items.size(), 1);

    // check tags
    EXPECT_EQ(list.items[0].tags.subject, "mypatient");
    EXPECT_EQ(list.items[0].tags.device, "mydevice");
    EXPECT_EQ(list.items[0].tags.session, session_id);
    EXPECT_EQ(list.items[0].tags.name, "myvariable1");
    EXPECT_EQ(list.items[0].tags.custom_tags.find("tagone")->second, "first");
    EXPECT_EQ(list.items[0].tags.custom_tags.find("tagtwo")->second, "second");
    auto multival_iterators = list.items[0].tags.custom_tags.equal_range("multival");
    std::vector<std::string> multival_values;
    std::transform(multival_iterators.first, multival_iterators.second, std::back_inserter(multival_values),
                   [](auto p) { return p.second; });
    ASSERT_THAT(multival_values, testing::ElementsAre("value1", "value2"));
}

TEST_F(StorageTest, storage_client_supports_paging) {
    StorageClient storage_client(storage_address);
    auto data = generate_random_vector(256);
    auto datastream = std::stringstream(std::string(data.begin(), data.end()), std::ios::binary);
    std::string session_id = "mysession4";
    auto tags = StorageItemTags::Builder("mypatient")
                    .with_device("mydevice")
                    .with_session(session_id)
                    .with_name("myvariable1")
                    .build();

    const int items_to_insert = 20;
    for (int i = 0; i < items_to_insert; i++) {
        storage_client.store_item(tags, datastream);
    }

    auto list = storage_client.list_items(tags, 5);
    EXPECT_FALSE(list.complete);
    EXPECT_EQ(list.items.size(), 5);

    for (int i = 0; i < 3; i++) {
        list = storage_client.get_next_page_of_items(list);
        EXPECT_EQ(list.items.size(), 5);
        EXPECT_EQ(list.complete, i == 2);
    }

    list = storage_client.list_items(tags, items_to_insert);
    EXPECT_TRUE(list.complete);
    EXPECT_EQ(list.items.size(), items_to_insert);
}

TEST_F(StorageTest, storage_client_should_store_items_and_return_list) {
    StorageClient storage_client(storage_address);
    auto data = generate_random_vector(256);
    auto datastream = std::stringstream(std::string(data.begin(), data.end()), std::ios::binary);
    std::string session_id = "mysession5";

    auto tags1 = StorageItemTags::Builder("mypatient")
                     .with_device("mydevice")
                     .with_session(session_id)
                     .with_name("myvariable1")
                     .build();

    auto tags2 = StorageItemTags::Builder("mypatient")
                     .with_device("mydevice")
                     .with_session(session_id)
                     .with_name("myvariable2")
                     .build();

    auto tags_session = StorageItemTags::Builder("mypatient").with_device("mydevice").with_session(session_id).build();

    auto resp = storage_client.store_item(tags1, datastream);

    datastream.clear();
    datastream.seekg(0);
    resp = storage_client.store_item(tags1, datastream);

    datastream.clear();
    datastream.seekg(0);
    resp = storage_client.store_item(tags1, datastream);

    datastream.clear();
    datastream.seekg(0);
    resp = storage_client.store_item(tags2, datastream);

    datastream.clear();
    datastream.seekg(0);
    resp = storage_client.store_item(tags2, datastream);

    auto list = storage_client.list_items(tags_session);
    EXPECT_EQ(list.items.size(), 5);
    list = storage_client.list_items(tags1);
    EXPECT_EQ(list.items.size(), 3);
    list = storage_client.list_items(tags2);
    EXPECT_EQ(list.items.size(), 2);
}

TEST_F(StorageTest, storage_client_supports_time_to_live) {
    StorageClient storage_client(storage_address);
    auto data = generate_random_vector(256);
    auto datastream = std::stringstream(std::string(data.begin(), data.end()), std::ios::binary);
    std::string session_id = "mysession6";

    auto tags = StorageItemTags::Builder("mypatient")
                    .with_device("mydevice")
                    .with_session(session_id)
                    .with_name("myname")
                    .build();

    auto resp = storage_client.store_item(tags, datastream, std::chrono::seconds(10));
    ASSERT_TRUE(resp.expires.has_value());

    resp = storage_client.store_item(tags, datastream);
    ASSERT_FALSE(resp.expires.has_value());
}

TEST_F(StorageTest, storage_client_should_return_empty_ptr_when_not_found) {
    StorageClient storage_client(storage_address);
    std::string session_id = "mysession7";
    auto tags = StorageItemTags::Builder("mypatient")
                    .with_device("mydevice")
                    .with_session(session_id)
                    .with_name("myvariable1")
                    .build();

    ASSERT_FALSE(storage_client.get_latest_item(tags));
    ASSERT_FALSE(storage_client.get_item_by_url(storage_address + "/foo"));;
}

TEST_F(StorageTest, storage_client_exception_contains_server_error_message) {
    StorageClient storage_client(storage_address);
    std::cout << storage_address << std::endl;
    auto tags = StorageItemTags::Builder("")
                    .build();

    auto datastream = std::stringstream("abc");
    EXPECT_THROW(
        try{
            storage_client.store_item(tags, datastream);
        } catch (std::runtime_error const& e) {
            EXPECT_THAT(e.what(), testing::HasSubstr("InvalidSubject"));
            throw;
        },
        std::runtime_error
    );
}

TEST_F(StorageTest, storage_client_handles_trailing_slash_in_address) {
    StorageClient storage_client(storage_address + "/");
    ASSERT_FALSE(storage_client.health_check().has_value());
}

TEST_F(StorageTest, storage_client_healthcheck_returns_error_when_address_invalid) {
    StorageClient storage_client(storage_address + "/foobar");
    ASSERT_TRUE(storage_client.health_check().has_value());
}

TEST_F(StorageTest, basic_storage) {
    hoNDArray<float> x(10);
    std::fill(x.begin(), x.end(), 23);
    storage.session->store("stuff", x);

    auto item = storage.session->get_latest<hoNDArray<float>>("stuff");
    ASSERT_TRUE(item.has_value());
    ASSERT_EQ(x, *item);

    hoNDArray<float> y(2, 2);
    std::fill(x.begin(), x.end(), 23);
    storage.session->store("stuff", y);

    item = storage.session->get_latest<hoNDArray<float>>("stuff");
    ASSERT_TRUE(item.has_value());
    ASSERT_EQ(y, *item);
}

TEST_F(StorageTest, larger_storage) {
    hoNDArray<float> x(1024 * 255);
    std::fill(x.begin(), x.end(), 23);
    storage.session->store("larger_storage_test", x);

    auto item = storage.session->get_latest<hoNDArray<float>>("larger_storage_test");
    ASSERT_EQ(x, *item);
}

TEST_F(StorageTest, image_test) {
    Core::Image<float> image;

    auto& [header, data, meta] = image;

    data = hoNDArray<float>(2, 2, 1, 4);
    std::fill(data.begin(), data.end(), 3);
    storage.session->store("image", image);

    auto item = storage.session->get_latest<Core::Image<float>>("image");

    auto [stored_header, stored_data, stored_meta] = *item;
    ASSERT_EQ(data, stored_data);
}