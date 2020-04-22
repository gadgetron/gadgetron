#pragma once

#include <string>

#include "Context.h"
#include <boost/process/child.hpp>
#include "Storage.h"

namespace Gadgetron::Server {

struct StorageServer {
  using Address = std::string;
  boost::process::child process;
  Address address;

};

StorageServer
start_storage_server(const boost::filesystem::path &working_directory);



Gadgetron::Core::Storage setup_storage(const StorageServer::Address& address, const Core::Context::Header& header);

} // namespace Gadgetron::Server
