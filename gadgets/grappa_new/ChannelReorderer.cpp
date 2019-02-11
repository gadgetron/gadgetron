#include "ChannelReorderer.h"

#include <set>
#include <regex>
#include <vector>
#include <boost/algorithm/string.hpp>

#include "hoNDArray_utils.h"

namespace {
    using namespace Gadgetron::Core;
}

namespace Gadgetron::Grappa
{
    ChannelReorderer::ChannelReorderer(
            const Gadgetron::Core::Context &context,
            const std::unordered_map<std::string, std::string> &props
    ) : TypedPureGadget<AnnotatedAcquisition, Acquisition>(props),
            context(context),
            channel_labels(build_channel_label_map()),
            uncombined_indices(parse_uncombined_channels()) {}

    AnnotatedAcquisition ChannelReorderer::process_function(Core::Acquisition acquisition) const {

        auto header     = std::get<ISMRMRD::AcquisitionHeader>(acquisition);
        auto trajectory = std::get<optional<hoNDArray<float>>>(acquisition);
        auto data       = std::get<hoNDArray<std::complex<float>>>(acquisition);

        auto reordered_data = data;
        auto reordered_header = header;

        auto reordering = create_reordering(header.available_channels);


        // TODO: Actually reorder.

        GINFO_STREAM("Reordering: (" << header.available_channels << ")");
        for (auto r : reordering) GINFO_STREAM("\t" << r);



        return AnnotatedAcquisition{
            reordered_header,
            trajectory,
            reordered_data,
            ChannelAnnotation{
                0, 0,
                std::move(reordering)
            }
        };
    }

    std::map<std::string, size_t> ChannelReorderer::build_channel_label_map() {

        if (!context.header.userParameters) return std::map<std::string, size_t>{};

        std::map<std::string, size_t> labels{};
        std::regex pattern{"COIL_(.*)"}; std::smatch match;

        for (auto &pair : context.header.userParameters->userParameterString) {
            if (std::regex_search(pair.name, match, pattern)) {
                labels[match[0]] = size_t(std::stoi(pair.value));
            }
        }

        return std::move(labels);
    }

    std::vector<size_t> ChannelReorderer::parse_uncombined_channels() {

        auto raw = uncombined_channels;

        std::set<size_t> uncombined{};
        std::vector<std::string> tokens{};
        boost::split(tokens, raw, boost::is_any_of(","));

        for (auto token : tokens) {
            boost::trim(token);
            uncombined.insert(parse_uncombined_channel(token));
        }

        std::vector<size_t> sorted{uncombined.begin(), uncombined.end()};
        std::sort(sorted.begin(), sorted.end());

        return sorted;
    }

    size_t ChannelReorderer::parse_uncombined_channel(const std::string &token) {

        std::smatch result{};

        std::regex label_pattern("'(.*)'");
        std::regex index_pattern("(\\d+)");

        if (std::regex_match(token, result, label_pattern)) {
            return channel_labels.at(result[0]);
        }

        if (std::regex_match(token, result, index_pattern)) {
            return size_t(std::stoi(result[0]));
        }

        throw std::runtime_error("Unable to parse uncombined channel list token: " + token);
    }

    std::vector<size_t> ChannelReorderer::create_reordering(size_t number_of_channels) const {

        std::vector<size_t> reordering(number_of_channels);
        std::iota(reordering.begin(), reordering.end(), 0);

        reordering.erase(std::remove_if(reordering.begin(), reordering.end(), [&](auto index) {
            return std::find(uncombined_indices.begin(), uncombined_indices.end(), index) != uncombined_indices.end();
        }), reordering.end());

        reordering.insert(reordering.end(), uncombined_indices.begin(), uncombined_indices.end());

        return reordering;
    }


    GADGETRON_GADGET_EXPORT(ChannelReorderer)
}
