#pragma once

#include <sstream>
#include <string>
#include <vector>
#include <boost/filesystem/path.hpp>
namespace Gadgetron::Core::IO {



    void from_string(const std::string&, int&);
    void from_string(const std::string&, unsigned int&);
    void from_string(const std::string&, unsigned short&);
    void from_string(const std::string&, short&);
    void from_string(const std::string&, size_t&);
    void from_string(const std::string&, double&);
    void from_string(const std::string&, float&);
    void from_string(const std::string&, bool&);


    void from_string(const std::string&, std::vector<int>&);
    void from_string(const std::string&, std::vector<unsigned int>&);
    void from_string(const std::string&, std::vector<unsigned short>&);
    void from_string(const std::string&, std::vector<short>&);
    void from_string(const std::string&, std::vector<size_t>&);
    void from_string(const std::string&,std::vector< double>&);
    void from_string(const std::string&, std::vector<float>&);
    void from_string(const std::string&, std::vector<bool>&);

    void from_string(const std::string&, boost::filesystem::path&);

    template <class T> T from_string(const std::string& str) {
        T val;
        from_string(str, val);
        return val;
    }
}
