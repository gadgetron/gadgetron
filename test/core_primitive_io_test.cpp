
#include <gtest/gtest.h>
#include "hoNDImage.h"
#include "hoMRImage.h"
#include "hoNDArray_elemwise.h"
#include "io/primitives.h"

TEST(TypeTests, readWritehoNDImageTest) {

    Gadgetron::hoNDImage<double, 2> im(12, 3);

    im.set_pixel_size(0, 1.2);
    im.set_pixel_size(1, 3.2);

    im.set_origin(0, 0.23);
    im.set_origin(1, 12.45);

    im.set_axis(0, 0, 0.8);
    im.set_axis(0, 1, 0.2);

    im.set_axis(1, 0, 0.1);
    im.set_axis(1, 1, 0.9);

    Gadgetron::fill(im, 3.42);

    std::ostringstream ostr;
    Gadgetron::Core::IO::write(ostr, im);

    std::string content = ostr.str();

    std::istringstream istr;
    istr.str(content);

    Gadgetron::hoNDImage<double, 2> im_read;
    Gadgetron::Core::IO::read(istr, im_read);

    EXPECT_EQ(im.get_size(0), im_read.get_size(0));
    EXPECT_EQ(im.get_size(1), im_read.get_size(1));

    EXPECT_DOUBLE_EQ(im.get_pixel_size(0), im_read.get_pixel_size(0));
    EXPECT_DOUBLE_EQ(im.get_pixel_size(1), im_read.get_pixel_size(1));

    EXPECT_DOUBLE_EQ(im.get_origin(0), im_read.get_origin(0));
    EXPECT_DOUBLE_EQ(im.get_origin(1), im_read.get_origin(1));

    EXPECT_DOUBLE_EQ(im.get_axis(0, 0), im_read.get_axis(0, 0));
    EXPECT_DOUBLE_EQ(im.get_axis(0, 1), im_read.get_axis(0, 1));

    EXPECT_DOUBLE_EQ(im.get_axis(1, 0), im_read.get_axis(1, 0));
    EXPECT_DOUBLE_EQ(im.get_axis(1, 1), im_read.get_axis(1, 1));

    EXPECT_DOUBLE_EQ(im_read(0, 0), 3.42);
    EXPECT_DOUBLE_EQ(im_read(0, 1), 3.42);
    EXPECT_DOUBLE_EQ(im_read(8, 2), 3.42);
    EXPECT_DOUBLE_EQ(im_read(11, 2), 3.42);
}
