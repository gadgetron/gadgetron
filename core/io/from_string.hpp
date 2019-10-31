//
// Created by dchansen on 10/14/19.
//

#pragma once
#include <boost/filesystem/path.hpp>
#include <bitset>
namespace Gadgetron::Core::IO {

    void from_string(const std::string& str, float& val);
    void from_string(const std::string& str, double& val);
    void from_string(const std::string& str, double& val);
    void from_string(const std::string& str, int& val);
    void from_string(const std::string& str, unsigned int& val);
    void from_string(const std::string& str, char& val);
    void from_string(const std::string& str, unsigned short& val);
    void from_string(const std::string& str, short& val);
    void from_string(const std::string& str, bool& val);
    void from_string(const std::string& str, size_t& val);

    void from_string(const std::string& str, std::vector<float>& val);
    void from_string(const std::string& str, std::vector<double>& val);
    void from_string(const std::string& str, std::vector<double>& val);
    void from_string(const std::string& str, std::vector<int>& val);
    void from_string(const std::string& str, std::vector<unsigned int>& val);
    void from_string(const std::string& str, std::vector<char>& val);
    void from_string(const std::string& str, std::vector<unsigned short>& val);
    void from_string(const std::string& str, std::vector<short>& val);
    void from_string(const std::string& str, std::vector<bool>& val);
    void from_string(const std::string& str, std::vector<size_t>& val);

    void from_string(const std::string& str, boost::filesystem::path& path);

    template<size_t N>
    void from_string(const std::string& str, std::bitset<N>& bset){
         bset = std::bitset<N>(from_string<unsigned int>(str));
    }

}


template <class T> T Gadgetron::Core::IO::from_string(const std::string& str) {
    T val;
    from_string(str, val);
    return val;
}
