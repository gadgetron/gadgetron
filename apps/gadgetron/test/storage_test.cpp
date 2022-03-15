
#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include "storage.h"

#include <chrono>
#include <random>
#include <thread>

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <range/v3/range.hpp>

using namespace Gadgetron::Storage;
using namespace Gadgetron;

using namespace boost::filesystem;
using namespace boost::program_options;

using testing::ElementsAre;

class StorageTest : public ::testing::Test {
  protected:
    static void SetUpTestSuite() {

        temp_dir = boost::filesystem::temp_directory_path() / "gadgetron_session_test";
        boost::filesystem::remove_all(temp_dir);
        boost::filesystem::create_directory(temp_dir);

        options_description desc("Storage options");
        desc.add_options()("storage_port,s", value<unsigned short>()->default_value(26589),
                           "Port on which to run the storage server.")(
            "database_dir,D", value<path>()->default_value(temp_dir),
            "Directory in which to store the storage server database.")(
            "storage_dir,S", value<path>()->default_value(temp_dir), "Directory in which to store data blobs.");

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
    std::transform(multival_iterators.first, multival_iterators.second, std::back_inserter(multival_values), [](auto p) { return p.second; });
    ASSERT_THAT(multival_values, ElementsAre("value1", "value2"));
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

    for (int i=0; i<3; i++) {
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

TEST_F(StorageTest, basic_storage) {
    hoNDArray<float> x(10);
    std::fill(x.begin(), x.end(), 23);
    storage.session.store("stuff", x);

    auto storage_list = storage.session.fetch<hoNDArray<float>>("stuff");
    ASSERT_EQ(storage_list.size(), 1);
    auto fetched = storage_list[0];
    ASSERT_EQ(x, fetched);

    hoNDArray<float> y(2, 2);
    std::fill(x.begin(), x.end(), 23);
    storage.session.store("stuff", y);

    storage_list = storage.session.fetch<hoNDArray<float>>("stuff");
    ASSERT_EQ(storage_list.size(), 2);
    fetched = storage_list[0];
    ASSERT_EQ(fetched, y);
    fetched = storage_list[1];
    ASSERT_EQ(fetched, x);
}

TEST_F(StorageTest, larger_storage) {
    hoNDArray<float> x(1024 * 255);
    std::fill(x.begin(), x.end(), 23);
    storage.session.store("larger_storage_test", x);

    auto storage_list = storage.session.fetch<hoNDArray<float>>("larger_storage_test");
    ASSERT_EQ(storage_list.size(), 1);
    auto fetched = storage_list[0];
    ASSERT_EQ(x, fetched);
}

TEST_F(StorageTest, range_test) {
    hoNDArray<float> y(2, 2);
    std::fill(y.begin(), y.end(), 23);
    storage.session.store("range_test", y);
    storage.session.store("range_test", y);

    auto storage_list = storage.session.fetch<hoNDArray<float>>("range_test");
    auto storage_vector = storage_list | ranges::to<std::vector>();
    auto storage_vector2 = storage.session.fetch<hoNDArray<float>>("range_test") | ranges::to<std::vector>();
    ASSERT_EQ(storage_list.size(), 2);
    ASSERT_EQ(storage_list.size(), std::distance(storage_list.begin(), storage_list.end()));
    ASSERT_EQ(storage_list.size(), storage_vector.size());
    ASSERT_EQ(storage_list.size(), storage_vector2.size());
}

TEST_F(StorageTest, image_test) {
    Core::Image<float> image;

    auto& [header, data, meta] = image;

    data = hoNDArray<float>(2, 2, 1, 4);
    std::fill(data.begin(), data.end(), 3);
    storage.session.store("image", image);

    auto storage_list = storage.session.fetch<Core::Image<float>>("image");

    auto [stored_header, stored_data, stored_meta] = storage_list[0];
    ASSERT_EQ(data, stored_data);
}