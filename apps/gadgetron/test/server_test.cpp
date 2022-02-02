
#include "gtest/gtest.h"

#include "storage.h"

#include <chrono>
#include <thread>

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <range/v3/range.hpp>

using namespace Gadgetron::Storage;
using namespace Gadgetron;

using namespace boost::filesystem;
using namespace boost::program_options;

class ServerTest : public ::testing::Test {
protected:
    static void SetUpTestSuite() {

        temp_dir = boost::filesystem::temp_directory_path() / "gadgetron_session_test";
        boost::filesystem::remove_all(temp_dir);
        boost::filesystem::create_directory(temp_dir);

        options_description desc("Storage options");
        desc.add_options()
                ("storage_port,s",
                 value<unsigned short>()->default_value(26589),
                 "Port on which to run the storage server.")
                ("database_dir,D",
                 value<path>()->default_value(temp_dir),
                 "Directory in which to store the storage server database.")
                ("storage_dir,S",
                 value<path>()->default_value(temp_dir),
                 "Directory in which to store data blobs.");

        variables_map args;
        store(parse_environment(desc, "GADGETRON_STORAGE_TEST_"), args);
        notify(args);

        auto [address, process] = Server::ensure_storage_server(args);

        // We need to make sure the storage server is up and ready before proceeding.
        std::this_thread::sleep_for(std::chrono::milliseconds(1000));

        ISMRMRD::IsmrmrdHeader header;
        header.subjectInformation = ISMRMRD::SubjectInformation{{},{},std::string("Penny the Pirate"),{},{}};
        header.studyInformation = ISMRMRD::StudyInformation{{},{},std::string("YAAARH")};

        server = std::move(*process);
        storage = Gadgetron::Server::setup_storage_spaces(address, header);
    }

    static void TearDownTestSuite() {
        boost::filesystem::remove_all(temp_dir);
        server.terminate();
        server.wait();
    }

    inline static StorageSpaces storage;
    inline static boost::process::child server;
    inline static boost::filesystem::path temp_dir;
};

TEST_F(ServerTest, do_nothing) {
    // This should pass.
}

TEST_F(ServerTest, basic_storage){
    hoNDArray<float> x(10);
    std::fill(x.begin(),x.end(),23);
    storage.session.store("stuff",x);

    auto storage_list = storage.session.fetch<hoNDArray<float>>("stuff");
    ASSERT_EQ(storage_list.size(),1);
    auto fetched = storage_list[0];
    ASSERT_EQ(x,fetched);

    hoNDArray<float> y(2,2);
    std::fill(x.begin(),x.end(),23);
    storage.session.store("stuff",y);

    storage_list = storage.session.fetch<hoNDArray<float>>("stuff");
    ASSERT_EQ(storage_list.size(),2);
    fetched = storage_list[0];
    ASSERT_EQ(fetched,y);
    fetched = storage_list[1];
    ASSERT_EQ(fetched,x);
}

TEST_F(ServerTest, larger_storage){
    hoNDArray<float> x(1024*255);
    std::fill(x.begin(),x.end(),23);
    storage.session.store("larger_storage_test",x);

    auto storage_list = storage.session.fetch<hoNDArray<float>>("larger_storage_test");
    ASSERT_EQ(storage_list.size(),1);
    auto fetched = storage_list[0];
    ASSERT_EQ(x,fetched);
}

TEST_F(ServerTest, range_test){
    hoNDArray<float> y(2,2);
    std::fill(y.begin(),y.end(),23);
    storage.session.store("range_test",y);
    storage.session.store("range_test",y);

    auto storage_list = storage.session.fetch<hoNDArray<float>>("range_test");
    auto storage_vector = storage_list | ranges::to<std::vector>();
    auto storage_vector2 = storage.session.fetch<hoNDArray<float>>("range_test") | ranges::to<std::vector>();
    ASSERT_EQ(storage_list.size(), 2);
    ASSERT_EQ(storage_list.size(), std::distance(storage_list.begin(),storage_list.end()));
    ASSERT_EQ(storage_list.size(), storage_vector.size());
    ASSERT_EQ(storage_list.size(), storage_vector2.size());
}

TEST_F(ServerTest, image_test){
    Core::Image<float> image;

    auto& [header,data,meta] = image;

    data = hoNDArray<float>(2,2,1,4);
    std::fill(data.begin(),data.end(),3);
    storage.session.store("image",image);

    auto storage_list = storage.session.fetch<Core::Image<float>>("image");

    auto [stored_header,stored_data,stored_meta] = storage_list[0];
    ASSERT_EQ(data,stored_data);
}