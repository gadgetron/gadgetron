#pragma once

#include "AnnotatedAcquisition.h"

namespace Gadgetron::Grappa {

    template<class T>
    bool is_last_in_slice(const T &acquisition) {
        return std::get<mrd::Acquisition>(acquisition).head.flags.HasFlags(mrd::AcquisitionFlags::kLastInSlice);
    }

    template<class T>
    uint16_t slice_of(const T &acquisition) {
        return std::get<mrd::Acquisition>(acquisition).head.idx.slice.value_or(0);
    }

    template<class T>
    uint16_t line_of(const T &acquisition) {
        return std::get<mrd::Acquisition>(acquisition).head.idx.kspace_encode_step_1.value_or(0);
    }

    template<class T>
    uint16_t samples_in(const T &acquisition) {
        return std::get<mrd::Acquisition>(acquisition).Samples();
    }

    template<class T>
    uint16_t channels_in(const T &acquisition) {
        return std::get<mrd::Acquisition>(acquisition).Coils();
    }

    template<class T>
    uint64_t combined_channels(const T &acquisition) {
        std::optional<ChannelAnnotation> optional_annotation = std::get<std::optional<ChannelAnnotation>>(acquisition);
        return optional_annotation ? optional_annotation->number_of_combined_channels : channels_in(acquisition);
    }

    template<class T>
    uint64_t uncombined_channels(const T &acquisition) {
        std::optional<ChannelAnnotation> optional_annotation = std::get<std::optional<ChannelAnnotation>>(acquisition);
        return optional_annotation ? optional_annotation->number_of_uncombined_channels : 0;
    }

    template<class Coll>
    std::string to_string(Coll collection) {

        if (collection.empty()) return "[]";

        std::stringstream stream;

        stream << "[";
        for (auto i : collection) stream << i << ", ";
        stream << "\b\b]";

        return stream.str();
    }
}