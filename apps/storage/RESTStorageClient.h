#pragma once

#include <ismrmrd/xml.h>
#include "Storage.h"
#include "Address.h"

namespace Gadgetron::Storage{


class RESTStorageClient : public Core::StorageSpace::StreamProvider {
    public:
        RESTStorageClient(const Address& address, const std::string &group,
                          const std::string &subject);

        ~RESTStorageClient() override = default;

        std::vector<std::string> content(const std::string &key) const override;
        std::vector<char> fetch(const std::string &key) const override;

        void store(const std::string &key, const std::vector<char> &value,
                   boost::posix_time::time_duration duration) const override;

    private:
        std::string port;
        std::string server_address;
        std::string group;
        std::string subject;
    };


    Core::Storage setup_storage(const Address& address, const ISMRMRD::IsmrmrdHeader& header );

}

