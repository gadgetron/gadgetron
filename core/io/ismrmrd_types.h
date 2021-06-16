//
// Created by dchansen on 2/25/19.
//

#pragma once
#include <iostream>
#include <ismrmrd/meta.h>
#include <ismrmrd/waveform.h>
#include <ismrmrd/xml.h>
#include <ismrmrd/ismrmrd.h>
#include "complext.h"
#include "Types.h"
namespace Gadgetron::Core::IO {

    template<class T> inline constexpr uint16_t ismrmrd_data_type(){ return T::this_function_is_not_defined; }
    template<> inline constexpr uint16_t ismrmrd_data_type<unsigned short>(){return ISMRMRD::ISMRMRD_USHORT;}
    template<> inline constexpr uint16_t ismrmrd_data_type<short>(){return ISMRMRD::ISMRMRD_SHORT;}
    template<> inline constexpr uint16_t ismrmrd_data_type<unsigned int>(){return ISMRMRD::ISMRMRD_UINT;}
    template<> inline constexpr uint16_t ismrmrd_data_type<int>(){return ISMRMRD::ISMRMRD_INT;}
    template<> inline constexpr uint16_t ismrmrd_data_type<float>(){return ISMRMRD::ISMRMRD_FLOAT;}
    template<> inline constexpr uint16_t ismrmrd_data_type<double>(){return ISMRMRD::ISMRMRD_DOUBLE;}
    template<> inline constexpr uint16_t ismrmrd_data_type<std::complex<float>>(){return ISMRMRD::ISMRMRD_CXFLOAT;}
    template<> inline constexpr uint16_t ismrmrd_data_type<std::complex<double>>(){return ISMRMRD::ISMRMRD_CXDOUBLE;}
    template<> inline constexpr uint16_t ismrmrd_data_type<complext<float>>(){return ISMRMRD::ISMRMRD_CXFLOAT;}
    template<> inline constexpr uint16_t ismrmrd_data_type<complext<double>>(){return ISMRMRD::ISMRMRD_CXDOUBLE;}

    inline void read(std::istream& stream, ISMRMRD::IsmrmrdHeader& header);
    inline void write(std::ostream& stream, const ISMRMRD::IsmrmrdHeader& header);

    inline void write(std::ostream& stream, const ISMRMRD::MetaContainer& meta);
    inline void read(std::istream& stream, ISMRMRD::MetaContainer& meta);

    inline void write(std::ostream& stream, const ISMRMRD::Waveform& wave);
    inline void read(std::istream& stream, ISMRMRD::Waveform& wave);

    template<class T>
    inline void read(std::istream& stream, Image<T>& img);
    template<class T>
    inline void write(std::ostream& stream, const Image<T>& img);

}

#include "primitives.h"
template<class T>
void Gadgetron::Core::IO::read(std::istream &stream, Image <T> &img) {
    auto& [header, data, meta] = img;
    header = IO::read<ISMRMRD::ImageHeader>(stream);
    if (header.attribute_string_len > 0 ){
        meta = ISMRMRD::MetaContainer();
        read(stream,*meta);
    } else {
        read_string_from_stream(stream);
    }
    data = hoNDArray<T>(header.matrix_size[0],header.matrix_size[1],header.matrix_size[2],header.channels);
    Gadgetron::Core::IO::read(stream,data.data(),data.size());
}

template<class T>
void Gadgetron::Core::IO::write(std::ostream& stream, const Image<T>& img) {

    const auto& [header, data, meta] = img;
    std::string serialized_meta;
    uint64_t meta_size = 0;

    if (meta) {
        std::stringstream meta_stream;
        ISMRMRD::serialize(*meta, meta_stream);
        serialized_meta = meta_stream.str();
        meta_size = serialized_meta.size() + 1;
    }

    ISMRMRD::ImageHeader corrected_header = header;
    corrected_header.data_type = ismrmrd_data_type<T>();
    corrected_header.attribute_string_len = meta_size;
    auto header_dims = vector_td<uint64_t,3>(corrected_header.matrix_size);

    if ((prod(header_dims)*corrected_header.channels) != data.size()){
        if (data.dimensions().size() != 4) throw std::runtime_error("Trying to write ISMRMRD Image, but data and header do not match");
        for (int i =0; i < 3; i++) corrected_header.matrix_size[i] = data.get_size(i);
        corrected_header.channels = data.get_size(3);
    }

    IO::write(stream, corrected_header);
    IO::write(stream, meta_size);
    stream.write(serialized_meta.c_str(), meta_size);
    IO::write(stream, data.get_data_ptr(), data.get_number_of_elements());

}
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