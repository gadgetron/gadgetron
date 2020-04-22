#pragma once

#include <iostream>
#include <memory>

#include "Context.h"
#include "storage_server.h"
#include "Server.h"

namespace Gadgetron::Server::Connection {
void handle(const Settings &settings, std::unique_ptr<std::iostream> stream);

struct StreamContext {
  Core::Context context;
  Settings settings;

  StreamContext(ISMRMRD::IsmrmrdHeader header, const Settings& settings_,
                Core::Storage storage)
      : context{std::move(header), settings_.paths, std::move(storage)}, settings(settings_) {}
};
} // namespace Gadgetron::Server::Connection

