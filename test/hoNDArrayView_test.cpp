#include <gtest/gtest.h>
// Created by dchansen on 4/12/19.
//
#include "hoNDArray.h"
#include "hoNDArray_math.h"
using namespace Gadgetron;
TEST(hoNDArrayView,copy){

    hoNDArray<float> x(127,49,3);
    std::fill(x.begin(),x.end(),0.0f);

    hoNDArray<float> y(127,3);
    std::fill(y.begin(),y.end(),1.0f);

    x(slice,2,slice) = y;

    ASSERT_EQ(asum(x),asum(y));


}

TEST(hoNDArrayView,copy2){

    hoNDArray<float> x(127,49,3);
    std::fill(x.begin(),x.end(),0.0f);

    hoNDArray<float> y(127,3);
    std::fill(y.begin(),y.end(),1.0f);

    x(slice,1,slice) = y;

    hoNDArray<float> z(127,3);
    std::fill(z.begin(),z.end(),2.0f);

    x(slice,2,slice) = z;
    ASSERT_EQ(asum(x),asum(y)*3);


}
