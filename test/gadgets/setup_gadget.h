//
// Created by dchansen on 9/20/19.
//
#pragma once

#include "Channel.h"
#include "Context.h"
#include "Node.h"
#include "PropertyMixin.h"
#include <mrd/types.h>
#include <array>
#include <thread>

namespace Gadgetron { namespace Test {

    inline mrd::EncodingSpaceType generate_encodingspace(std::array<unsigned short, 3> matrix_size, std::array<float, 3> fov) {
        return mrd::EncodingSpaceType{ { matrix_size[0], matrix_size[1], matrix_size[2] }, { fov[0], fov[1], fov[2] } };
    }
    inline mrd::EncodingType generate_encoding() {
        auto encoding         = mrd::EncodingType{};
        encoding.encoded_space = generate_encodingspace({ 192, 192, 1 }, { 256, 256, 10 });
        encoding.encoding_limits.kspace_encoding_step_0 = mrd::LimitType{};
        encoding.encoding_limits.kspace_encoding_step_1 = mrd::LimitType{};
        encoding.encoding_limits.kspace_encoding_step_0->minimum = 0;
        encoding.encoding_limits.kspace_encoding_step_0->center = 96;
        encoding.encoding_limits.kspace_encoding_step_0->maximum = 191;
        encoding.encoding_limits.kspace_encoding_step_1->minimum = 0;
        encoding.encoding_limits.kspace_encoding_step_1->center = 96;
        encoding.encoding_limits.kspace_encoding_step_1->maximum = 191;

        encoding.recon_space   = generate_encodingspace({ 128, 128, 1 }, { 256, 256, 10 });
        return encoding;
    }
    inline mrd::Header generate_header() {
        auto header                                           = mrd::Header{};
        header.encoding                                       = { generate_encoding() };
        header.experimental_conditions.h1resonance_frequency_hz = 63.87 * 1e6;
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


    inline mrd::Acquisition generate_acquisition(size_t number_of_samples, size_t channels, size_t measurement_uid = 42){
        mrd::Acquisition acq;
        acq.head.channel_order = std::vector<uint32_t>(channels, 1);
        acq.head.center_sample = number_of_samples / 2;
        acq.head.measurement_uid = measurement_uid;
        acq.data.create({number_of_samples, channels});
        std::fill(acq.data.begin(), acq.data.end(), std::complex<float>(1));
        return acq;
    }

}}
