#include "ChannelReorderer.h"

#include <set>
#include <regex>
#include <vector>
#include <boost/algorithm/string.hpp>

#include "hoNDArray_iterators.h"
#include "hoNDArray_utils.h"

namespace {

    template<class T>
    std::vector<T> reorder(std::vector<T> v, const std::vector<uint64_t> &reordering) {
        std::vector<T> reordered; reordered.reserve(reordering.size());
        for (auto i : reordering) reordered.push_back(v[i]);
        return reordered;
    }

    mrd::AcquisitionHeader reorder(const mrd::AcquisitionHeader &header, const std::vector<uint64_t> &reordering) {
        auto reordered_header = header;

        for (uint16_t out_idx = 0; out_idx < reordering.size(); out_idx++) {
            auto in_idx = uint16_t(reordering[out_idx]);
        }

        return reordered_header;
    }
}

namespace Gadgetron::Grappa
{
    ChannelReorderer::ChannelReorderer(
            const Gadgetron::Core::Context &context,
            const std::unordered_map<std::string, std::string> &props
    ) : PureGadget<AnnotatedAcquisition, mrd::Acquisition>(context,props),
            context(context),
            channel_labels(build_channel_label_map()),
            uncombined_indices(parse_uncombined_channels()) {}

    AnnotatedAcquisition ChannelReorderer::process_function(mrd::Acquisition acquisition) const {

        auto header = acquisition.head;
        auto trajectory = acquisition.trajectory;
        auto data = acquisition.data;

        auto reordering = create_reordering(acquisition.Coils());

        auto channels = spans(data, 1);
        auto reordered_channels = reorder(
                std::vector<hoNDArray<std::complex<float>>>(channels.begin(), channels.end()),
                reordering
        );

        acquisition.head = reorder(header, reordering);
        acquisition.data = concat(reordered_channels);
        return AnnotatedAcquisition{acquisition, ChannelAnnotation{acquisition.Coils() - uncombined_indices.size(), uncombined_indices.size(), std::move(reordering)}};
    }

    std::map<std::string, size_t> ChannelReorderer::build_channel_label_map() {

        if (!context.header.user_parameters) return std::map<std::string, size_t>{};

        std::map<std::string, size_t> labels{};
        std::regex pattern{"COIL_(.*)"}; std::smatch match;

        for (auto &pair : context.header.user_parameters->user_parameter_string) {
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
            if (!token.empty()) uncombined.insert(parse_uncombined_channel(token));
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

    std::vector<uint64_t> ChannelReorderer::create_reordering(size_t number_of_channels) const {

        std::vector<uint64_t> reordering(number_of_channels);
        std::iota(reordering.begin(), reordering.end(), 0);

        reordering.erase(std::remove_if(reordering.begin(), reordering.end(), [&](auto index) {
            return std::find(uncombined_indices.begin(), uncombined_indices.end(), index) != uncombined_indices.end();
        }), reordering.end());

        reordering.insert(reordering.end(), uncombined_indices.begin(), uncombined_indices.end());

        return reordering;
    }

    GADGETRON_GADGET_EXPORT(ChannelReorderer)
}
