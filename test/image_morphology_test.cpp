/** \file       image_morphology_test.cpp
    \brief      Test case for image processing morphological functions

    \author     Hui Xue
*/

#include "morphology.h"
#include <gtest/gtest.h>

using namespace Gadgetron;
using testing::Types;

template<typename T> class image_morphology_test : public ::testing::Test
{
protected:
    virtual void SetUp()
    {
    }
};

typedef Types<float> realImplementations;
TYPED_TEST_SUITE(image_morphology_test, realImplementations);

TYPED_TEST(image_morphology_test, bwlabel)
{
    // define an image with four regions
    size_t RO = 365;
    size_t E1 = 187;

    hoNDArray<TypeParam> aImage;
    aImage.create(RO, E1);
    Gadgetron::clear(aImage);

    // fill in four regions

    std::vector<size_t> start_ro(4), end_ro(4), start_e1(4), end_e1(4);

    start_ro[0] = 0;
    start_ro[1] = 25;
    start_ro[2] = 66;
    start_ro[3] = 287;

    end_ro[0] = 12;
    end_ro[1] = 56;
    end_ro[2] = 88;
    end_ro[3] = 364;

    start_e1[0] = 0;
    start_e1[1] = 38;
    start_e1[2] = 76;
    start_e1[3] = 121;

    end_e1[0] = 36;
    end_e1[1] = 72;
    end_e1[2] = 94;
    end_e1[3] = 186;

    size_t n, ro, e1;

    for (n = 0; n < start_ro.size(); n++)
    {
        for (e1 = start_e1[n]; e1 <= end_e1[n]; e1++)
        {
            for (ro = start_ro[n]; ro <= end_ro[n]; ro++)
            {
                aImage(ro, e1) = 1;
            }
        }
    }

    hoNDArray<unsigned int> label;
    bool is_8_connected = true;
    Gadgetron::bwlabel_2d(aImage, (TypeParam)1, label, is_8_connected);

    std::vector<unsigned int> labels;
    std::vector<unsigned int> areas;
    Gadgetron::bwlabel_area_2d(label, labels, areas);

    EXPECT_EQ(labels.size(), 4);
    EXPECT_EQ(areas.size(), 4);
    EXPECT_EQ(areas[0], (end_ro[0] - start_ro[0] + 1)*(end_e1[0] - start_e1[0] + 1));
    EXPECT_EQ(areas[1], (end_ro[1] - start_ro[1] + 1)*(end_e1[1] - start_e1[1] + 1));
    EXPECT_EQ(areas[2], (end_ro[2] - start_ro[2] + 1)*(end_e1[2] - start_e1[2] + 1));
    EXPECT_EQ(areas[3], (end_ro[3] - start_ro[3] + 1)*(end_e1[3] - start_e1[3] + 1));
}
