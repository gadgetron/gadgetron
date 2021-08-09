#include "Storage.h"
#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/device/back_inserter.hpp>
#include <boost/iostreams/stream_buffer.hpp>
#include <boost/iostreams/stream.hpp>

namespace bio = boost::iostreams;
namespace Gadgetron::Storage {
    GenericStorageSpace::GenericStorageSpace(
        std::shared_ptr<StreamProvider> provider,
        const Core::optional<std::string>& subject,
        boost::posix_time::time_duration default_duration
    ) : provider(std::move(provider)), subject(subject), default_duration{default_duration} {}

    std::unique_ptr<std::istream> istream_from_data(const std::vector<char> &data) {

            return std::make_unique<bio::stream<bio::array_source>>(data.data(),data.size());

    }

    std::unique_ptr<std::ostream> ostream_view(std::vector<char> &data) {
        return std::make_unique<bio::stream<bio::back_insert_device<std::vector<char>>>>(data);
    }
}
