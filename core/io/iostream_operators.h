#pragma once
#include "Types.h"
#include <ostream>


namespace std { 
    template<class T, class... TYPES>
    std::ostream& operator<<(std::ostream& stream, const Gadgetron::Core::variant<T,TYPES...>& value){
        visit([&stream](auto&& arg){
            stream << arg;
        }, value);
        return stream;
    }
}