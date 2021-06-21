#pragma once

#include <boost/program_options.hpp>
#include <optional>

#include "StorageServer.h"

namespace Gadgetron::Server {
    std::tuple<std::string, std::optional<Gadgetron::Storage::StorageServer>>
    ensure_storage_server(const boost::program_options::variables_map &args);
}
