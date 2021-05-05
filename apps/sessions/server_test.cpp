//
// Created by dch on 5/5/21.
//

#include "gtest/gtest.h"
#include "SessionServer.h"
#include "RESTStorageClient.h"

using namespace Gadgetron::Sessions;
using namespace Gadgetron;

class ServerTest : public ::testing::Test {
protected:
    void SetUp() override {
        temp_dir = std::filesystem::temp_directory_path() / "gadgetron_session_test";
        std::filesystem::remove_all(temp_dir);
        std::filesystem::create_directory(temp_dir);
        server = std::make_unique<SessionServer>(0,temp_dir/"database", temp_dir);
        ISMRMRD::IsmrmrdHeader header;
        header.subjectInformation = ISMRMRD::SubjectInformation{{},{},std::string("Penny the Pirate"),{},{}};
        storage = Storage::setup_storage("localhost",std::to_string(server->port()),header);

    }

    void TearDown() override {
        server = nullptr;
        std::filesystem::remove_all(temp_dir);
    }

    Gadgetron::Core::Storage storage;
    std::unique_ptr<SessionServer> server;
    std::filesystem::path temp_dir;

};

TEST_F(ServerTest,basic_storage){
    hoNDArray<float> x(10);
    std::fill(x.begin(),x.end(),23);
    this->storage.session.store("stuff",x);

    auto storage_list = this->storage.session.fetch<hoNDArray<float>>("stuff");

    ASSERT_EQ(storage_list.size(),1);

    auto fetched = storage_list[0];
    ASSERT_EQ(x,fetched);


}

