#pragma once

#include <ostream>
#include <type_traits>
#include "hoNDArray.h"

namespace Gadgetron::Core::IO {

    template<class T, class V = std::enable_if_t<std::is_trivially_copyable_v<T>>> void write(std::ostream& stream, const T& value){
        stream.write(reinterpret_cast<const char*>(&value), sizeof(value));
    }


    template<class T> void write(std::ostream& stream, const hoNDArray<T>& array ){
        stream.write(reinterpret_cast<const char*>(array.get_data_ptr()),array.get_number_of_bytes());
    }

}