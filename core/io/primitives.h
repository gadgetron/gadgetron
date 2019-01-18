#pragma once

#include <iostream>
#include <type_traits>

#include "hoNDArray.h"

namespace Gadgetron::Core::IO {

    template<class T, class V = std::enable_if<std::is_trivially_copyable_v<T>>>
    void read(std::istream& stream, T& value ){
        stream.read(reinterpret_cast<char*>(&value),sizeof(value));

        // TODO: Check error bits on stream and throw on failure.
    }

    template<class T, class V = std::enable_if<std::is_trivially_copyable_v<T>>>
    T read(std::istream& stream ){
        T value;
        read(stream, value);
        return value;
    }

    template<class T> void read(std::istream& stream, hoNDArray<T>& array){
        stream.read(reinterpret_cast<char*>(array.get_data_ptr()), array.get_number_of_bytes());
    }

    template<class T>
    std::string read_string_from_stream(std::istream &stream) {

        auto n = read<T>(stream);

        std::string str(n, 0);
        stream.read(str.data(), n);

        return str;
    }

    template<class T, class V = std::enable_if_t<std::is_trivially_copyable_v<T>>> void write(std::ostream& stream, const T& value){
        stream.write(reinterpret_cast<const char*>(&value), sizeof(value));
    }


    template<class T> void write(std::ostream& stream, const hoNDArray<T>& array ){
        stream.write(reinterpret_cast<const char*>(array.get_data_ptr()),array.get_number_of_bytes());
    }

    template<class SIZE_TYPE>
    void write_string_to_stream(std::ostream &stream, const std::string& str) {

        SIZE_TYPE string_length = static_cast<SIZE_TYPE>(str.size());
        write(stream,string_length);
        stream.write(str.data(),string_length);


    }
}