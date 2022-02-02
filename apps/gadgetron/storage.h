#pragma once

#include <optional>

#include <ismrmrd/xml.h>

#include <boost/process.hpp>
#include <boost/program_options.hpp>

#include "Storage.h"

namespace Gadgetron::Server {

    std::tuple<std::string, std::optional<boost::process::child>>
    ensure_storage_server(const boost::program_options::variables_map &args);

    StorageSpaces
    setup_storage_spaces(const std::string& address, const ISMRMRD::IsmrmrdHeader& header);

    class StorageClient : public Gadgetron::Storage::StreamProvider {
      public:
        StorageClient(std::string address, std::string group);
        ~StorageClient() override = default;

        [[nodiscard]] std::vector<std::string>
        content(const std::string& identifier, const std::string &name) const override;

        [[nodiscard]] std::vector<char>
        fetch(const std::string &uri) const override;

        void store(const std::string& identifier,
                   const std::string &name,
                   const std::vector<char> &data,
                   boost::posix_time::time_duration duration) override;

      private:
        std::string address;
        std::string space;
    };
}
