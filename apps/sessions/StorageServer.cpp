//
// Created by dch on 4/12/21.
//
#include "SessionServer.h"

int main(int argc, char** argv){

    auto server = Gadgetron::Sessions::start_session_server(9200);
    server.join();
}
