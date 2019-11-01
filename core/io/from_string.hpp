//
// Created by dchansen on 10/14/19.
//

#pragma once
#include <boost/filesystem/path.hpp>
#include <bitset>
namespace Gadgetron::Core::IO {


    template<size_t N>
    void from_string(const std::string& str, std::bitset<N>& bset){
         bset = std::bitset<N>(from_string<unsigned int>(str));
    }

}
