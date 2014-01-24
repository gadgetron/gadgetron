/*
 * tests.cpp
 *
 *  Created on: Feb 28, 2013
 *      Author: Dae
 */

#include <gtest/gtest.h>

int main(int argc, char **argv)
{
    //::testing::GTEST_FLAG(filter) = "*grappa*:*spirit*";

    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
