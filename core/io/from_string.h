#pragma once

#include <sstream>
#include <string>
#include <vector>
namespace Gadgetron::Core::IO {


    template<class T> void from_string(const std::string& str, T& va);
    template<class T> void from_string(const std::string& str, std::vector<T>& va);


    template <class T> T from_string(const std::string& str) {
        T val;
        from_string(str, val);
        return val;
    }
}
