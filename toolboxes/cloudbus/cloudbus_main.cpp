#include <iostream>

#include "CloudBus.h"

int main(int argc, char** argv)
{
  std::cout << "CloudBus Main Program" << std::endl;

  int port = GADGETRON_DEFAULT_MULTICAST_PORT;
  const char* addr = GADGETRON_DEFAULT_MULTICAST_ADDR;
  bool query_only_mode = true;

  if (argc > 1) {
    addr = argv[1];
    std::cout << "Setting multicast address to: " << addr << std::endl;
  }

  if (argc > 2) {
    port = std::atoi(argv[2]);
    std::cout << "Setting multicast port to: " << port << std::endl;
  }

  if (argc > 3)
  {
    query_only_mode = false;
  }

  //Port and address must be set before grabbing the instance for the first time. 
  Gadgetron::CloudBus::set_mcast_address(addr);
  Gadgetron::CloudBus::set_mcast_port(port);
  Gadgetron::CloudBus::set_query_only(query_only_mode);
  Gadgetron::CloudBus* cb = Gadgetron::CloudBus::instance();
  cb->wait();
  return 0;
}
