
#include "gtest/gtest.h"

#include "StorageServer.h"
#include "storage.h"

#include <range/v3/range.hpp>

using namespace Gadgetron::Storage;
using namespace Gadgetron;


class ServerTest : public ::testing::Test {
protected:
    void SetUp() override {
        temp_dir = boost::filesystem::temp_directory_path() / "gadgetron_session_test";
        boost::filesystem::remove_all(temp_dir);
        boost::filesystem::create_directory(temp_dir);
        server = std::make_unique<StorageServer>(0,temp_dir/"database", temp_dir);
        ISMRMRD::IsmrmrdHeader header;
        header.subjectInformation = ISMRMRD::SubjectInformation{{},{},std::string("Penny the Pirate"),{},{}};
        header.studyInformation = ISMRMRD::StudyInformation{{},{},std::string("YAAARH")};
        storage = Gadgetron::Server::setup_storage_spaces("https://localhost:" + std::to_string(server->port()), header);
    }

    void TearDown() override {
        server = nullptr;
        boost::filesystem::remove_all(temp_dir);
    }

    StorageSpaces storage;
    std::unique_ptr<StorageServer> server;
    boost::filesystem::path temp_dir;

};

TEST_F(ServerTest,basic_storage){
    hoNDArray<float> x(10);
    std::fill(x.begin(),x.end(),23);
    this->storage.session.store("stuff",x);

    auto storage_list = this->storage.session.fetch<hoNDArray<float>>("stuff");

    ASSERT_EQ(storage_list.size(),1);

    auto fetched = storage_list[0];
    ASSERT_EQ(x,fetched);
    hoNDArray<float> y(2,2);
    std::fill(x.begin(),x.end(),23);
    this->storage.session.store("stuff",y);

    storage_list = this->storage.session.fetch<hoNDArray<float>>("stuff");

    ASSERT_EQ(storage_list.size(),2);
    fetched = storage_list[0];
    ASSERT_EQ(fetched,y);
    fetched = storage_list[1];
    ASSERT_EQ(fetched,x);
}

TEST_F(ServerTest,larger_storage){
    hoNDArray<float> x(1024*255);
    std::fill(x.begin(),x.end(),23);
    this->storage.session.store("stuff",x);

    auto storage_list = this->storage.session.fetch<hoNDArray<float>>("stuff");

    ASSERT_EQ(storage_list.size(),1);

    auto fetched = storage_list[0];
    ASSERT_EQ(x,fetched);


}
TEST_F(ServerTest,range_test){


    hoNDArray<float> y(2,2);
    std::fill(y.begin(),y.end(),23);
    this->storage.session.store("stuff",y);
    this->storage.session.store("stuff",y);

    auto storage_list = this->storage.session.fetch<hoNDArray<float>>("stuff");

    auto storage_vector = storage_list | ranges::to<std::vector>();


    auto storage_vector2 = this->storage.session.fetch<hoNDArray<float>>("stuff") | ranges::to<std::vector>();
    ASSERT_EQ(storage_list.size(),2);
    ASSERT_EQ(storage_list.size(), std::distance(storage_list.begin(),storage_list.end()));
    ASSERT_EQ(storage_list.size(),storage_vector.size());
    ASSERT_EQ(storage_list.size(),storage_vector2.size());


}

TEST_F(ServerTest,image_test){
    Core::Image<float> image;

    auto& [header,data,meta] = image;

    data = hoNDArray<float>(2,2,1,4);
    std::fill(data.begin(),data.end(),3);
    this->storage.session.store("image",image);

    auto storage_list = this->storage.session.fetch<Core::Image<float>>("image");

    auto [stored_header,stored_data,stored_meta] = storage_list[0];
    ASSERT_EQ(data,stored_data);
}