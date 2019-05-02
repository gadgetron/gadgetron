//
// Created by dchansen on 2/21/19.
//

#include <gtest/gtest.h>
#include "ThreadPool.h"

using namespace Gadgetron::Core;
TEST(ThreadPoolTest,VoidTest){
    ThreadPool pool{4};
    auto return_value = pool.async([](){});
    return_value.get();
    pool.join();

}


TEST(ThreadPoolTest,moveTest){
    ThreadPool pool{4};
    auto a = std::make_unique<int>(5);
    auto return_value = pool.async([](auto b){return *b;},std::move(a));
    auto val = return_value.get();
    EXPECT_EQ(val,5);
    pool.join();

}
