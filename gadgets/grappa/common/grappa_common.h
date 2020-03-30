#pragma once

#include "AnnotatedAcquisition.h"
#include "Types.h"

namespace Gadgetron::Grappa {

    template<class T>
    bool is_last_in_slice(const T &acquisition) {
        return std::get<ISMRMRD::AcquisitionHeader>(acquisition).isFlagSet(ISMRMRD::ISMRMRD_ACQ_LAST_IN_SLICE);
    }

    template<class T>
    uint16_t slice_of(const T &acquisition) {
        return std::get<ISMRMRD::AcquisitionHeader>(acquisition).idx.slice;
    }

    template<class T>
    uint16_t line_of(const T &acquisition) {
        return std::get<ISMRMRD::AcquisitionHeader>(acquisition).idx.kspace_encode_step_1;
    }

    template<class T>
    uint16_t samples_in(const T &acquisition) {
        return std::get<ISMRMRD::AcquisitionHeader>(acquisition).number_of_samples;
    }

    template<class T>
    uint16_t channels_in(const T &acquisition) {
        return std::get<ISMRMRD::AcquisitionHeader>(acquisition).active_channels;
    }

    template<class T>
    uint64_t combined_channels(const T &acquisition) {
        Core::optional<ChannelAnnotation> optional_annotation = std::get<Core::optional<ChannelAnnotation>>(acquisition);
        return optional_annotation ? optional_annotation->number_of_combined_channels : channels_in(acquisition);
    }

    template<class T>
    uint64_t uncombined_channels(const T &acquisition) {
        Core::optional<ChannelAnnotation> optional_annotation = std::get<Core::optional<ChannelAnnotation>>(acquisition);
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