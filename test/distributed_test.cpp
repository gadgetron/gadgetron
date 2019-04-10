#include <gtest/gtest.h>
#include "Message.h"
#include "Channel.h"
#include "Types.h"



#include "../apps/gadgetron/connection/distributed/remote_workers.h"


#include <stdlib.h>

TEST(distributed,get_remote_workers){


    std::string test_string = R"(GADGETRON_REMOTE_WORKER_COMMAND= echo ["localhost:9002", "[::1]:9003"])";
    std::vector<char> string_copy(test_string.data(),test_string.data()+test_string.size());
    putenv(string_copy.data());

    auto workers = Gadgetron::Server::Distributed::get_remote_workers();

    EXPECT_EQ(workers.size(),2);
    EXPECT_EQ(workers[0].ip,"localhost");
    EXPECT_EQ(workers[0].port,"9002");
    EXPECT_EQ(workers[1].ip,"::1");
    EXPECT_EQ(workers[1].port,"9003");
}



