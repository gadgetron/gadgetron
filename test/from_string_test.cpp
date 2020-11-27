//
// Created by dchansen on 2/21/19.
//

#include <gtest/gtest.h>
#include "io/from_string.h"
using namespace Gadgetron::Core;
TEST(FromStringTest,FloatTest){


    ASSERT_EQ(42.0f, IO::from_string<float>("42"));
    ASSERT_EQ(42.0f, IO::from_string<float>("42.0"));
    ASSERT_EQ(42.0f, IO::from_string<float>("42.0f"));
    ASSERT_EQ(42.0f, IO::from_string<float>("42e0"));

}
TEST(FromStringTest,VectorTest){

    std::vector<int> vals = {42,43,44,45};

    auto result = IO::from_string<std::vector<int>>("42 43 44 45");

    ASSERT_EQ(vals,result);

    auto result2 = IO::from_string<std::vector<int>>("42, 43, 44, 45");

    ASSERT_EQ(vals,result2);
}

TEST(FromStringTest,BoolTest){
    ASSERT_EQ(true, IO::from_string<bool>("true"));
    ASSERT_EQ(false, IO::from_string<bool>("false"));
}
