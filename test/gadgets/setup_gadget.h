//
// Created by dchansen on 9/20/19.
//
#pragma once

#include "Channel.h"
#include "Context.h"
#include "Node.h"
#include "PropertyMixin.h"
#include <array>
#include <ismrmrd/ismrmrd.h>
#include <ismrmrd/xml.h>
#include <thread>

namespace Gadgetron { namespace Test {

    inline ISMRMRD::EncodingSpace generate_encodingspace(std::array<unsigned short, 3> matrix_size, std::array<float, 3> fov) {
        return ISMRMRD::EncodingSpace{ { matrix_size[0], matrix_size[1], matrix_size[2] }, { fov[0], fov[1], fov[2] } };
    }
    inline ISMRMRD::Encoding generate_encoding() {
        auto encoding         = ISMRMRD::Encoding{};
        encoding.encodedSpace = generate_encodingspace({ 192, 192, 1 }, { 256, 256, 10 });
        encoding.encodingLimits.kspace_encoding_step_0 = ISMRMRD::Limit();
        encoding.encodingLimits.kspace_encoding_step_1 = ISMRMRD::Limit();
        encoding.encodingLimits.kspace_encoding_step_0->minimum = 0;
        encoding.encodingLimits.kspace_encoding_step_0->center = 96;
        encoding.encodingLimits.kspace_encoding_step_0->maximum = 191;
        encoding.encodingLimits.kspace_encoding_step_1->minimum = 0;
        encoding.encodingLimits.kspace_encoding_step_1->center = 96;
        encoding.encodingLimits.kspace_encoding_step_1->maximum = 191;

        encoding.reconSpace   = generate_encodingspace({ 128, 128, 1 }, { 256, 256, 10 });
        return encoding;
    }
    inline ISMRMRD::IsmrmrdHeader generate_header() {
        auto header                                           = ISMRMRD::IsmrmrdHeader{};
        header.encoding                                       = { generate_encoding() };
        header.experimentalConditions.H1resonanceFrequency_Hz = 63.87 * 1e6;
        return header;
    }

    inline Core::Context generate_context() {
        Core::Context context;
        context.header = generate_header();
        return context;
    }

    template <class GADGET> struct GadgetChannels {
        Core::OutputChannel input;
        Core::GenericInputChannel output;
    };

    template <class GADGET>
    inline GadgetChannels<GADGET> setup_gadget(Core::GadgetProperties properties, Core::Context context = generate_context()) {

        auto channels  = Core::make_channel();
        auto channels2 = Core::make_channel();

        auto thread = std::thread(
            [](auto input, auto output, auto properties, auto context) {
                try {
                    GADGET gadget(context, properties);
                    Core::Node& node = gadget;
                    node.process(input, output);
                } catch (const Core::ChannelClosed&){}
            },
            std::move(channels.input), std::move(channels2.output), properties, context);

        thread.detach();
        return { std::move(channels.output), std::move(channels2.input) };
    }


    inline Core::Acquisition generate_acquisition(size_t number_of_samples, size_t channels, size_t measurement_uid = 42){
       auto header = ISMRMRD::AcquisitionHeader();
       header.number_of_samples = number_of_samples;
       header.active_channels = channels;
       header.available_channels = channels;
       header.center_sample = number_of_samples / 2;
       header.measurement_uid = measurement_uid;

       hoNDArray<std::complex<float>> data(number_of_samples,channels);
       std::fill(data.begin(),data.end(),std::complex<float>(1));

       return {header,data,Core::none};


    }

}}
