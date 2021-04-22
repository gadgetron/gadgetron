//
// Created by dch on 4/21/21.
//

#include "gtest/gtest.h"
#include "DB.h"

#include <filesystem>

using namespace Gadgetron;
using namespace Gadgetron::Sessions::DB;
class SessionsTest : public ::testing::Test {
protected:
    void SetUp() override {
        temp_dir = std::filesystem::temp_directory_path() / "gadgetron_session_test";
        std::filesystem::create_directory(temp_dir);
        db = std::make_unique<DB>(temp_dir);

    }

    void TearDown() override {
        db = nullptr;
        std::filesystem::remove_all(temp_dir);

    }

    std::unique_ptr<DB> db;
    std::filesystem::path temp_dir;

};

TEST(temporary,basics){
    rocksdb::DB* db;
    rocksdb::DBOptions options;
    options.create_if_missing= true;
    options.create_missing_column_families = true;
    std::vector<rocksdb::ColumnFamilyDescriptor> cfd = { {"default",rocksdb::ColumnFamilyOptions()}, {"test",rocksdb::ColumnFamilyOptions()}};
    std::vector<rocksdb::ColumnFamilyHandle*> handles;
    auto status = rocksdb::DB::Open(options,"/tmp/testdb",cfd, &handles, &db);
    EXPECT_TRUE(status.ok());

    EXPECT_EQ(handles[1]->GetName(), "test");

    status = db->Put(rocksdb::WriteOptions{},handles[1],"penguins","I can be anything");
    EXPECT_TRUE(status.ok());
    std::string return_value;
    status = db->Get(rocksdb::ReadOptions{},handles[1],"penguins",&return_value);
    EXPECT_TRUE(status.ok());
    EXPECT_EQ(return_value,"I can be anything");

    delete db;


}


TEST_F(SessionsTest,pendingwrites){

    using clock = boost::posix_time::second_clock;
    auto meta = BlobMeta{"penguin",clock::universal_time(),clock::universal_time()};

    db->pending_writes.set("monkey",{clock::universal_time(),meta});

    auto retrieved = db->pending_writes["monkey"];

    EXPECT_TRUE(retrieved);
    EXPECT_EQ(meta.blob_id,retrieved->meta.blob_id);
    EXPECT_EQ(meta.creation_time,retrieved->meta.creation_time);
    EXPECT_EQ(meta.deletion_time,retrieved->meta.deletion_time);

}

TEST_F(SessionsTest,pendingwrites_delete){

    using clock = boost::posix_time::second_clock;
    auto meta = BlobMeta{"penguin",clock::universal_time(),clock::universal_time()};

    auto retrieved0 = db->pending_writes["monkey"];
    EXPECT_FALSE(retrieved0);

    db->pending_writes.set("monkey",{clock::universal_time(),meta});

    auto retrieved = db->pending_writes["monkey"];

    EXPECT_TRUE(retrieved);

    db->pending_writes.delete_key("monkey");

    auto retrieved2  = db->pending_writes["monkey"];

    EXPECT_FALSE(retrieved2);
}

TEST_F(SessionsTest,range_test){

    auto& info = db->db_info;

    info.set("/sessions/penguin/0",{1});
    info.set("/sessions/penguin/00",{2});
    info.set("/sessions/apenguin/0",{0});

    std::vector<json> values;
    for (auto [key,value] : info){
        values.push_back(value);
    }

    EXPECT_EQ(values[0],json{0});
    EXPECT_EQ(values[1],json{1});
    EXPECT_EQ(values[2],json{2});

}
TEST_F(SessionsTest,blobs_append) {

    using clock = boost::posix_time::second_clock;
    auto meta = BlobMeta{"penguin", clock::universal_time(), clock::universal_time()};
    auto meta2 = BlobMeta{"penguin2", clock::universal_time(), clock::universal_time()};

    db->blobs.push_back("blob", meta);
    db->blobs.push_back("blob", meta2);

    auto metas = db->blobs["blob"];
    ASSERT_EQ(metas.size(),2);
    EXPECT_EQ(metas[0],meta);
    EXPECT_EQ(metas[1],meta2);
}



