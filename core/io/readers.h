#pragma once

#include <istream>
#include <type_traits>

#include "hoNDArray.h"

namespace Gadgetron::Core::IO {

    template<class T, class V = std::enable_if<std::is_trivially_copyable_v<T>>>
    void read(std::istream& stream, T& value ){
        stream.read(reinterpret_cast<char*>(&value),sizeof(value));
    }

   template<class T, class V = std::enable_if<std::is_trivially_copyable_v<T>>>
    T read(std::istream& stream ){
        T value;
        stream.read(reinterpret_cast<char*>(&value),sizeof(value));
        return value;
    }

    template<class T> void read(std::istream& stream, hoNDArray<T>& array){
        stream.read(reinterpret_cast<char*>(array.get_data_ptr()), array.get_number_of_bytes());
    }
}