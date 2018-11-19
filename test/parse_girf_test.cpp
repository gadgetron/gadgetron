//
// Created by dchansen on 9/26/18.
//

#include <gtest/gtest.h>
#include "mri_core_girf_correction.h"



TEST(parse,parse_girf_test){


    std::string test_data = "[ [2+5i, \n 3+1i],[3+8i , 2+5i], [1+1i,1+1i]]";


    auto data = Gadgetron::GIRF::load_girf_kernel(test_data);

    EXPECT_EQ(data.get_size(1),3);
    EXPECT_EQ(data.get_size(0),2);

}
