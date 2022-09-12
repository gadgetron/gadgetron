#include <gtest/gtest.h>
// Created by dchansen on 4/12/19.
//
#include "hoNDArray.h"
#include "hoNDArray_math.h"
using namespace Gadgetron;
using namespace Gadgetron::Indexing;
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


TEST(hoNDArrayView,copy3){

    hoNDArray<float> x(2,4,3);
    std::fill(x.begin(),x.end(),0.0f);

    hoNDArray<float> y(2,4);
    std::fill(y.begin(),y.end(),1.0f);

    x(slice,slice,1) = y;

    hoNDArray<float> z(2,4);
    std::fill(z.begin(),z.end(),2.0f);

    x(slice,slice,2) = z;
    ASSERT_EQ(asum(x),asum(y)*3);
}

TEST(hoNDArrayView,copy4){

    hoNDArray<float> x(2,4,3);
    std::fill(x.begin(),x.end(),0.0f);

    hoNDArray<float> y(2);
    std::fill(y.begin(),y.end(),1.0f);

    x(slice,1,2) = y;

    hoNDArray<float> z(2);
    std::fill(z.begin(),z.end(),2.0f);

    x(slice,3,2) = z;
    ASSERT_EQ(asum(x),asum(y)*3);
}

TEST(hoNDArrayView,assignment3){
    hoNDArray<float> x(2,4,3);
    std::fill(x.begin(),x.end(),0.0f);


    hoNDArray<float> y(27,4,3);

    std::fill(y.begin(),y.end(),1.0f);
    x(1,slice,2) = y(15,slice,0);

    ASSERT_FLOAT_EQ(asum(x),4.0f);

}

TEST(hoNDArrayView,assignment){
    hoNDArray<float> x(2,4,3);
    std::fill(x.begin(),x.end(),0.0f);


    hoNDArray<float> y(2,4,5);

    std::fill(y.begin(),y.end(),2.0f);

    x(slice,slice,2) = y(slice,slice,0);

    ASSERT_FLOAT_EQ(asum(x),16.0f);

    x(slice,slice,1) = y(slice,slice,4);
    ASSERT_FLOAT_EQ(asum(x),32.0f);

}

TEST(hoNDArrayView,assignment2){
    hoNDArray<float> x(4,3);
    std::fill(x.begin(),x.end(),2.0f);


    hoNDArray<float> y(27,4,3);

    std::fill(y.begin(),y.end(),0.0f);

    y(3,slice,slice) = x;

    ASSERT_FLOAT_EQ(asum(y),24.0f);

    y(1,slice,slice) = x;

    ASSERT_FLOAT_EQ(asum(y),48.0f);

}

TEST(hoNDArrayView,assignment4){
    hoNDArray<float> x(4,27,3);
    std::fill(x.begin(),x.end(),2.0f);


    hoNDArray<float> y(27,4,3);

    std::fill(y.begin(),y.end(),0.0f);

    y(slice,1,slice) = x(3,slice,slice);

    ASSERT_FLOAT_EQ(asum(y),27*3*2.0f);


}
TEST(hoNDArrayView,failure){
    hoNDArray<float> x(4,27,3);
    std::fill(x.begin(),x.end(),2.0f);


    hoNDArray<float> y(27,4,3);

    std::fill(y.begin(),y.end(),0.0f);

    ASSERT_ANY_THROW(y(slice,1,slice) = x(slice,slice,3));
}