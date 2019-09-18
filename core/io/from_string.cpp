//
// Created by dchansen on 2/22/19.
//
#include "from_string.h"

#include <boost/spirit/include/qi.hpp>
#include <boost/filesystem/path.hpp>

template <class T> void Gadgetron::Core::IO::from_string(const std::string& str, T& val) {

    namespace qi    = boost::spirit::qi;
    namespace ascii = boost::spirit::ascii;
    using ascii::space;
    bool success = qi::phrase_parse(str.begin(), str.end(), qi::auto_, space, val);
    if (!success)
        throw std::runtime_error("Cannot parse value from string: " + str);
}

template <class T> void Gadgetron::Core::IO::from_string(const std::string& str, std::vector<T>& vals) {
    namespace qi    = boost::spirit::qi;
    namespace ascii = boost::spirit::ascii;
    using ascii::space;
    bool success = qi::phrase_parse(str.begin(), str.end(), (*(qi::auto_ >> (-qi::lit(',')))), space, vals);
    if (!success)
        throw std::runtime_error("Cannot parse vector from string: " + str);
}

template<>
void Gadgetron::Core::IO::from_string(const std::string& str, boost::filesystem::path& path){
    path = boost::filesystem::path(str);
}

template void Gadgetron::Core::IO::from_string<float>(const std::string& str, float& val);
template void Gadgetron::Core::IO::from_string<double>(const std::string& str, double& val);
template void Gadgetron::Core::IO::from_string<unsigned int>(const std::string& str, unsigned int& val);
template void Gadgetron::Core::IO::from_string<int>(const std::string& str, int& val);
template void Gadgetron::Core::IO::from_string<unsigned short>(const std::string& str, unsigned short& val);
template void Gadgetron::Core::IO::from_string<short>(const std::string& str, short& val);
template void Gadgetron::Core::IO::from_string<size_t>(const std::string& str, size_t& val);
template void Gadgetron::Core::IO::from_string<bool>(const std::string& str, bool& val);

template void Gadgetron::Core::IO::from_string<float>(const std::string& str, std::vector<float>& val);
template void Gadgetron::Core::IO::from_string<double>(const std::string& str, std::vector<double>& val);
template void Gadgetron::Core::IO::from_string<unsigned int>(const std::string& str, std::vector<unsigned int>& val);
template void Gadgetron::Core::IO::from_string<int>(const std::string& str, std::vector<int>& val);
template void Gadgetron::Core::IO::from_string<unsigned short>(const std::string& str, std::vector<unsigned short>& val);
template void Gadgetron::Core::IO::from_string<short>(const std::string& str, std::vector<short>& val);
template void Gadgetron::Core::IO::from_string<size_t>(const std::string& str, std::vector<size_t>& val);
template void Gadgetron::Core::IO::from_string<bool>(const std::string& str, std::vector<bool>& val);

template void Gadgetron::Core::IO::from_string<boost::filesystem::path>(const std::string& str, boost::filesystem::path& val);
