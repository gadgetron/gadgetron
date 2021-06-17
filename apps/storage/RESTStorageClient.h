#pragma once

#include <ismrmrd/xml.h>
#include "Storage.h"
#include "Address.h"

namespace Gadgetron::Storage{


class RESTStorageClient : public StreamProvider {
    public:
        RESTStorageClient(std::string host, std::string service, std::string group);

        ~RESTStorageClient() override = default;

        std::vector<std::string> content(const std::string& subject, const std::string &key) const override;
        std::vector<char> fetch(const std::string &uuid) const override;

        void store(const std::string& subject, const std::string &key, const std::vector<char> &value,
                   boost::posix_time::time_duration duration) override;

    private:
        std::string host, service;
        std::string group;
    };

    StorageSpaces setup_storage(const std::string& address, const ISMRMRD::IsmrmrdHeader& header );
}

