
#include <gtest/gtest.h>
#include "hoNDImage.h"
#include "hoMRImage.h"
#include "hoNDArray_elemwise.h"
#include "io/primitives.h"
#include "cmr_file_and_directory_handling.h"

TEST(TypeTests, hoNDImageTest) {

    Gadgetron::hoNDImage<double, 2> im2d(12, 3);

    im2d.set_pixel_size(0, 1.2);
    im2d.set_pixel_size(1, 3.2);

    im2d.set_origin(0, 0.23);
    im2d.set_origin(1, 12.45);

    im2d.set_axis(0, 0, 0.8);
    im2d.set_axis(0, 1, 0.2);

    im2d.set_axis(1, 0, 0.1);
    im2d.set_axis(1, 1, 0.9);

    Gadgetron::fill(im2d, 3.42);

    std::ostringstream ostr;
    Gadgetron::Core::IO::write(ostr, im2d);

    std::string content = ostr.str();

    std::istringstream istr;
    istr.str(content);

    Gadgetron::hoNDImage<double, 2> im2d_readed;
    Gadgetron::Core::IO::read(istr, im2d_readed);

    EXPECT_EQ(im2d.get_size(0), im2d_readed.get_size(0));
    EXPECT_EQ(im2d.get_size(1), im2d_readed.get_size(1));

    EXPECT_DOUBLE_EQ(im2d.get_pixel_size(0), im2d_readed.get_pixel_size(0));
    EXPECT_DOUBLE_EQ(im2d.get_pixel_size(1), im2d_readed.get_pixel_size(1));

    EXPECT_DOUBLE_EQ(im2d.get_origin(0), im2d_readed.get_origin(0));
    EXPECT_DOUBLE_EQ(im2d.get_origin(1), im2d_readed.get_origin(1));

    EXPECT_DOUBLE_EQ(im2d.get_axis(0, 0), im2d_readed.get_axis(0, 0));
    EXPECT_DOUBLE_EQ(im2d.get_axis(0, 1), im2d_readed.get_axis(0, 1));

    EXPECT_DOUBLE_EQ(im2d.get_axis(1, 0), im2d_readed.get_axis(1, 0));
    EXPECT_DOUBLE_EQ(im2d.get_axis(1, 1), im2d_readed.get_axis(1, 1));

    EXPECT_DOUBLE_EQ(im2d_readed(0, 0), 3.42);
    EXPECT_DOUBLE_EQ(im2d_readed(0, 1), 3.42);
    EXPECT_DOUBLE_EQ(im2d_readed(8, 2), 3.42);
    EXPECT_DOUBLE_EQ(im2d_readed(11, 2), 3.42);
}

TEST(TypeTests, hoMRImageTest) {

    Gadgetron::hoMRImage<double, 2> im2d(12, 3);

    im2d.set_pixel_size(0, 1.2);
    im2d.set_pixel_size(1, 3.2);

    im2d.set_origin(0, 0.23);
    im2d.set_origin(1, 12.45);

    im2d.set_axis(0, 0, 0.8);
    im2d.set_axis(0, 1, 0.2);

    im2d.set_axis(1, 0, 0.1);
    im2d.set_axis(1, 1, 0.9);

    memset(&im2d.header_, 0, sizeof(ISMRMRD::ISMRMRD_ImageHeader));

    im2d.header_.flags = 11;
    im2d.header_.patient_table_position[0] = 1.1;
    im2d.header_.patient_table_position[1] = 2.1;
    im2d.header_.patient_table_position[2] = 3.1;

    Gadgetron::fill(im2d, 3.42);

    im2d.attrib_.set("Test1", 1.0);
    im2d.attrib_.append("Test1", 2.0);

    im2d.attrib_.set("Test2", (long)2);
    im2d.attrib_.set("Test3", "This is a test");

    std::ostringstream ostr;
    Gadgetron::Core::IO::write(ostr, im2d);
    std::string content = ostr.str();

    std::stringstream str;
    ISMRMRD::serialize( const_cast<ISMRMRD::MetaContainer&>(im2d.attrib_), str);
    std::string attrib_content = str.str();

    GDEBUG_STREAM("im2d.attrib_ is " << attrib_content)

    std::istringstream istr;
    istr.str(content);

    Gadgetron::hoMRImage<double, 2> im2d_readed;
    Gadgetron::Core::IO::read(istr, im2d_readed);

    std::stringstream str_readed;
    ISMRMRD::serialize( const_cast<ISMRMRD::MetaContainer&>(im2d_readed.attrib_), str_readed);
    std::string attrib_content_readed = str_readed.str();
    GDEBUG_STREAM("im2d_readed.attrib_ is " << attrib_content_readed)

    EXPECT_EQ(im2d.get_size(0), im2d_readed.get_size(0));
    EXPECT_EQ(im2d.get_size(1), im2d_readed.get_size(1));

    EXPECT_DOUBLE_EQ(im2d.get_pixel_size(0), im2d_readed.get_pixel_size(0));
    EXPECT_DOUBLE_EQ(im2d.get_pixel_size(1), im2d_readed.get_pixel_size(1));

    EXPECT_DOUBLE_EQ(im2d.get_origin(0), im2d_readed.get_origin(0));
    EXPECT_DOUBLE_EQ(im2d.get_origin(1), im2d_readed.get_origin(1));

    EXPECT_DOUBLE_EQ(im2d.get_axis(0, 0), im2d_readed.get_axis(0, 0));
    EXPECT_DOUBLE_EQ(im2d.get_axis(0, 1), im2d_readed.get_axis(0, 1));

    EXPECT_DOUBLE_EQ(im2d.get_axis(1, 0), im2d_readed.get_axis(1, 0));
    EXPECT_DOUBLE_EQ(im2d.get_axis(1, 1), im2d_readed.get_axis(1, 1));

    EXPECT_DOUBLE_EQ(im2d_readed(0, 0), 3.42);
    EXPECT_DOUBLE_EQ(im2d_readed(0, 1), 3.42);
    EXPECT_DOUBLE_EQ(im2d_readed(8, 2), 3.42);
    EXPECT_DOUBLE_EQ(im2d_readed(11, 2), 3.42);

    EXPECT_EQ(im2d.header_.flags, im2d_readed.header_.flags);
    EXPECT_EQ(im2d.header_.patient_table_position[0], im2d_readed.header_.patient_table_position[0]);
    EXPECT_EQ(im2d.header_.patient_table_position[1], im2d_readed.header_.patient_table_position[1]);
    EXPECT_EQ(im2d.header_.patient_table_position[2], im2d_readed.header_.patient_table_position[2]);

    EXPECT_STREQ(attrib_content.c_str(), attrib_content_readed.c_str());
}

TEST(TypeTests, hoNDArrayObjectTest) {

    typedef Gadgetron::hoMRImage<float, 3> ImageType;

    ImageType im(12, 3, 4);

    im.set_pixel_size(0, 1.2);
    im.set_pixel_size(1, 3.2);
    im.set_pixel_size(2, 2.2);

    im.set_origin(0, 0.23);
    im.set_origin(1, 12.45);
    im.set_origin(2, 3.45);

    im.set_axis(0, 0, 0.8);
    im.set_axis(0, 1, 0.2);
    im.set_axis(0, 2, 0.1);

    im.set_axis(1, 0, 0.1);
    im.set_axis(1, 1, 0.9);
    im.set_axis(1, 2, 0.0);

    im.set_axis(2, 0, 0.85);
    im.set_axis(2, 1, 0.01);
    im.set_axis(2, 2, 0.03);

    memset(&im.header_, 0, sizeof(ISMRMRD::ISMRMRD_ImageHeader));

    im.header_.flags = 11;
    im.header_.patient_table_position[0] = 1.1;
    im.header_.patient_table_position[1] = 2.1;
    im.header_.patient_table_position[2] = 3.1;

    Gadgetron::fill(im, 3.42f);

    im.attrib_.set("hoNDObjectArray_Test1", 1.0);
    im.attrib_.append("hoNDObjectArray_Test1", 2.0);

    im.attrib_.set("hoNDObjectArray_Test2", (long)2);
    im.attrib_.set("hoNDObjectArray_Test3", "This is a test");

    Gadgetron::hoNDArray<ImageType> array;
    array.create(3, 4);
    for (auto i=0; i<12; i++)
        array[i] = im;

    std::ostringstream ostr;
    Gadgetron::Core::IO::write(ostr, array);
    std::string content = ostr.str();

    std::stringstream str;
    ISMRMRD::serialize( const_cast<ISMRMRD::MetaContainer&>(im.attrib_), str);
    std::string attrib_content = str.str();

    GDEBUG_STREAM("im.attrib_ is \n" << attrib_content)

    std::istringstream istr;
    istr.str(content);

    Gadgetron::hoNDArray<ImageType> array_readed;
    Gadgetron::Core::IO::read(istr, array_readed);

    for (auto i=0; i<12; i++) {
        std::stringstream str_readed;
        ISMRMRD::serialize( const_cast<ISMRMRD::MetaContainer&>(array_readed[i].attrib_), str_readed);
        std::string attrib_content_readed = str_readed.str();
        if(i==0) {
            GDEBUG_STREAM("array_readed[0].attrib_ is \n" << attrib_content_readed)
        }

        ImageType im_readed = array_readed[i];

        EXPECT_EQ(im.get_size(0), im_readed.get_size(0));
        EXPECT_EQ(im.get_size(1), im_readed.get_size(1));
        EXPECT_EQ(im.get_size(2), im_readed.get_size(2));

        EXPECT_DOUBLE_EQ(im.get_pixel_size(0), im_readed.get_pixel_size(0));
        EXPECT_DOUBLE_EQ(im.get_pixel_size(1), im_readed.get_pixel_size(1));
        EXPECT_DOUBLE_EQ(im.get_pixel_size(2), im_readed.get_pixel_size(2));

        EXPECT_DOUBLE_EQ(im.get_origin(0), im_readed.get_origin(0));
        EXPECT_DOUBLE_EQ(im.get_origin(1), im_readed.get_origin(1));
        EXPECT_DOUBLE_EQ(im.get_origin(2), im_readed.get_origin(2));

        EXPECT_DOUBLE_EQ(im.get_axis(0, 0), im_readed.get_axis(0, 0));
        EXPECT_DOUBLE_EQ(im.get_axis(0, 1), im_readed.get_axis(0, 1));
        EXPECT_DOUBLE_EQ(im.get_axis(0, 2), im_readed.get_axis(0, 2));

        EXPECT_DOUBLE_EQ(im.get_axis(1, 0), im_readed.get_axis(1, 0));
        EXPECT_DOUBLE_EQ(im.get_axis(1, 1), im_readed.get_axis(1, 1));
        EXPECT_DOUBLE_EQ(im.get_axis(1, 2), im_readed.get_axis(1, 2));

        EXPECT_DOUBLE_EQ(im.get_axis(2, 0), im_readed.get_axis(2, 0));
        EXPECT_DOUBLE_EQ(im.get_axis(2, 1), im_readed.get_axis(2, 1));
        EXPECT_DOUBLE_EQ(im.get_axis(2, 2), im_readed.get_axis(2, 2));

        EXPECT_FLOAT_EQ(im_readed(0, 0), 3.42);
        EXPECT_FLOAT_EQ(im_readed(0, 1), 3.42);
        EXPECT_FLOAT_EQ(im_readed(8, 2), 3.42);
        EXPECT_FLOAT_EQ(im_readed(11, 2), 3.42);

        EXPECT_EQ(im.header_.flags, im_readed.header_.flags);
        EXPECT_EQ(im.header_.patient_table_position[0], im_readed.header_.patient_table_position[0]);
        EXPECT_EQ(im.header_.patient_table_position[1], im_readed.header_.patient_table_position[1]);
        EXPECT_EQ(im.header_.patient_table_position[2], im_readed.header_.patient_table_position[2]);

        EXPECT_STREQ(attrib_content.c_str(), attrib_content_readed.c_str());
    }

    GDEBUG_STREAM("Test save and load item ... ");
    std::string workingdirectory = "/tmp/gadgetron";
    std::string session_id = "23456";

    Gadgetron::save_item(workingdirectory, session_id, array);

    std::vector<std::string> items_list(1, session_id);
    std::vector< Gadgetron::hoNDArray<ImageType> > arrays;
    Gadgetron::load_items(workingdirectory, items_list, arrays);
}
