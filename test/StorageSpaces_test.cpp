#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <ismrmrd/xml.h>

#include <stdlib.h>

#include "IsmrmrdContextVariables.h"
#include "storage.h"

using namespace Gadgetron;
using namespace Gadgetron::Storage;
using namespace testing;

class MockStorageClient : public StorageClient {
  public:
    using StorageClient::StorageClient;

    MOCK_METHOD(StorageItemList, list_items, (StorageItemTags const& tags, size_t limit), (override));

    MOCK_METHOD(StorageItemList, get_next_page_of_items, (StorageItemList const& page), (override));

    MOCK_METHOD(std::shared_ptr<std::istream>, get_latest_item, (StorageItemTags const& tags), (override));

    MOCK_METHOD(std::shared_ptr<std::istream>, get_item_by_url, (std::string const& url), (override));

    MOCK_METHOD(StorageItem, store_item,
                (StorageItemTags const& tags, std::istream& data, std::optional<std::chrono::seconds> time_to_live),
                (override));

    MOCK_METHOD(std::optional<std::string>, health_check, (), (override));
};

bool tags_equal(StorageItemTags a, StorageItemTags b) {
    return a.subject == b.subject && a.device == b.device && a.session == b.session && a.name == b.name &&
           std::equal(a.custom_tags.begin(), a.custom_tags.end(), b.custom_tags.begin(), b.custom_tags.end());
}

MATCHER_P(TagsEq, cmp, "") { return tags_equal(arg, cmp); }

TEST(StorageSpacesTest, session_space_all_specified) {
    IsmrmrdContextVariables vars("mysubject", "mydevice", "mysession", "mymeasurement");

    auto expected_tags = StorageItemTags::Builder(vars.subject_id())
                             .with_device(vars.device_id())
                             .with_session(vars.session_id())
                             .with_name("myname")
                             .build();

    auto client = std::make_shared<MockStorageClient>("address");
    EXPECT_CALL(*client, store_item(TagsEq(expected_tags), _, Eq(std::chrono::seconds(3600))));
    EXPECT_CALL(*client, get_latest_item(TagsEq(expected_tags)));

    SessionSpace space(client, vars, std::chrono::hours(1));
    space.store("myname", "mydata");
    space.get_latest<std::vector<char>>("myname");
}

TEST(StorageSpacesTest, measurement_space_all_specified) {
    IsmrmrdContextVariables vars("mysubject", "mydevice", "mysession", "mymeasurement");

    auto expected_store_tags = StorageItemTags::Builder(vars.subject_id())
                                   .with_device(vars.device_id())
                                   .with_session(vars.session_id())
                                   .with_name("myname")
                                   .with_custom_tag("measurement", "mymeasurement")
                                   .build();

    auto expected_read_tags = StorageItemTags::Builder(vars.subject_id())
                                  .with_device(vars.device_id())
                                  .with_session(vars.session_id())
                                  .with_name("myname")
                                  .with_custom_tag("measurement", "mydependency")
                                  .build();

    auto client = std::make_shared<MockStorageClient>("address");
    EXPECT_CALL(*client, store_item(TagsEq(expected_store_tags), _, Eq(std::chrono::seconds(3600))));
    EXPECT_CALL(*client, get_latest_item(TagsEq(expected_read_tags)));

    MeasurementSpace space(client, vars, std::chrono::hours(1));
    space.store("myname", "mydata");
    space.get_latest<std::vector<char>>("mydependency", "myname");
}

TEST(StorageSpacesTest, scanner_space_all_specified) {
    IsmrmrdContextVariables vars("mysubject", "mydevice", "mysession", "mymeasurement");

    auto expected_tags = StorageItemTags::Builder("$null").with_device(vars.device_id()).with_name("myname").build();

    auto client = std::make_shared<MockStorageClient>("address");
    EXPECT_CALL(*client, store_item(TagsEq(expected_tags), _, Eq(std::chrono::seconds(3600))));
    EXPECT_CALL(*client, get_latest_item(TagsEq(expected_tags)));

    ScannerSpace space(client, vars, std::chrono::hours(1));
    space.store("myname", "mydata");
    space.get_latest<std::vector<char>>("myname");
}

TEST(StorageSpacesTest, session_space_fields_missing) {
    IsmrmrdContextVariables vars("", "", "", "");
    auto client = std::make_shared<MockStorageClient>("address");

    ScannerSpace space(client, vars, std::chrono::hours(1));

    EXPECT_THROW(space.store("myname", "mydata"), IncompleteStorageContextException);
    EXPECT_FALSE(space.get_latest<std::vector<char>>("myname").has_value());
}

TEST(StorageSpacesTest, measurement_space_subject_missing) {
    IsmrmrdContextVariables vars("", "mydevice", "mysession", "mymeasurement");
    auto client = std::make_shared<MockStorageClient>("address");

    MeasurementSpace space(client, vars, std::chrono::hours(1));

    EXPECT_THROW(space.store("myname", "mydata"), IncompleteStorageContextException);
    EXPECT_FALSE(space.get_latest<std::vector<char>>("mydependency", "myname").has_value());
}

TEST(StorageSpacesTest, measurement_space_measurement_missing) {
    IsmrmrdContextVariables vars("measurement", "mydevice", "mysession", "");
    auto client = std::make_shared<MockStorageClient>("address");

    auto expected_read_tags = StorageItemTags::Builder(vars.subject_id())
                                  .with_device(vars.device_id())
                                  .with_session(vars.session_id())
                                  .with_name("myname")
                                  .with_custom_tag("measurement", "mydependency")
                                  .build();

    MeasurementSpace space(client, vars, std::chrono::hours(1));

    EXPECT_THROW(space.store("myname", "mydata"), IncompleteStorageContextException);

    EXPECT_CALL(*client, get_latest_item(TagsEq(expected_read_tags)));
    space.get_latest<std::vector<char>>("mydependency", "myname");
}

TEST(StorageSpacesTest, scanner_space_device_missing) {
    IsmrmrdContextVariables vars("mysubject", "", "mysession", "mymeasurement");
    auto client = std::make_shared<MockStorageClient>("address");

    ScannerSpace space(client, vars, std::chrono::hours(1));

    EXPECT_THROW(space.store("myname", "mydata"), IncompleteStorageContextException);
    EXPECT_FALSE(space.get_latest<std::vector<char>>("myname").has_value());
}
