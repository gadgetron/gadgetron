//
// Created by dchansen on 2/25/19.
//

#pragma once
#include <iostream>
#include <ismrmrd/meta.h>
#include <ismrmrd/waveform.h>
#include <ismrmrd/xml.h>

namespace Gadgetron::Core::IO {

    inline void read(std::istream& stream, ISMRMRD::IsmrmrdHeader& header);
    inline void write(std::ostream& stream, const ISMRMRD::IsmrmrdHeader& header);

    inline void write(std::ostream& stream, const ISMRMRD::MetaContainer& meta);
    inline void read(std::istream& stream, ISMRMRD::MetaContainer& meta);

    inline void write(std::ostream& stream, const ISMRMRD::Waveform& wave);
    inline void read(std::istream& stream, ISMRMRD::Waveform& wave);

}

#include "primitives.h"

void Gadgetron::Core::IO::write(std::ostream& stream, const ISMRMRD::MetaContainer& meta) {
    std::stringstream meta_stream;
    ISMRMRD::serialize(meta, meta_stream);
    write_string_to_stream(stream, meta_stream.str());
}
void Gadgetron::Core::IO::read(std::istream& stream, ISMRMRD::MetaContainer& meta) {
    auto meta_string = read_string_from_stream(stream);
    ISMRMRD::deserialize(meta_string.c_str(), meta);
}
void Gadgetron::Core::IO::write(std::ostream& stream, const ISMRMRD::Waveform& wave) {
        IO::write(stream,wave.head);
        IO::write(stream,wave.begin_data(),wave.size());
}

void Gadgetron::Core::IO::read(std::istream& stream, ISMRMRD::Waveform& wave) {

        auto header = IO::read<ISMRMRD::WaveformHeader>(stream);

        wave = ISMRMRD::Waveform(header.number_of_samples, header.channels);

        IO::read(stream, wave.begin_data(),wave.size());

        wave.head = header;
}

void Gadgetron::Core::IO::read(std::istream& stream, ISMRMRD::IsmrmrdHeader& header){
    auto xml = read_string_from_stream(stream);
    ISMRMRD::deserialize(xml.c_str(),header);
}

void Gadgetron::Core::IO::write(std::ostream& stream, const ISMRMRD::IsmrmrdHeader& header){
    std::stringstream output_stream;
    ISMRMRD::serialize(header,output_stream);
    write_string_to_stream(stream, output_stream.str());
}