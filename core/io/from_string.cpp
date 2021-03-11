//
// Created by dchansen on 2/22/19.
//
#include "from_string.h"

#include <boost/spirit/include/qi.hpp>
#include <boost/filesystem/path.hpp>
namespace {
    template <class T> void from_string_impl(const std::string& str, T& val) {

        namespace qi    = boost::spirit::qi;
        namespace ascii = boost::spirit::ascii;
        using ascii::space;
        bool success = qi::phrase_parse(str.begin(), str.end(), qi::auto_, space, val);
        if (!success)
            throw std::runtime_error("Cannot parse value from string: " + str);
    }

    template <class T> void from_string_impl(const std::string& str, std::vector<T>& vals) {
        namespace qi    = boost::spirit::qi;
        namespace ascii = boost::spirit::ascii;
        using ascii::space;
        bool success = qi::phrase_parse(str.begin(), str.end(), (*(qi::auto_ >> (-qi::lit(',')))), space, vals);
        if (!success)
            throw std::runtime_error("Cannot parse vector from string: " + str);
    }
}

void Gadgetron::Core::IO::from_string(const std::string& str, boost::filesystem::path& path){
    path = boost::filesystem::path(str);
}

void Gadgetron::Core::IO::from_string(const std::string& str, long long& val) { from_string_impl(str, val); }

void Gadgetron::Core::IO::from_string(const std::string& str, double& val) { from_string_impl(str,val);}
void Gadgetron::Core::IO::from_string(const std::string& str, bool& val) { from_string_impl(str,val);}
void Gadgetron::Core::IO::from_string(const std::string& str, std::vector<long long>& val) { from_string_impl(str,val);}
void Gadgetron::Core::IO::from_string(const std::string& str, std::vector<double>& val) { from_string_impl(str,val);}
void Gadgetron::Core::IO::from_string(const std::string& str, std::vector<bool>& val) { from_string_impl(str,val);}
void Gadgetron::Core::IO::from_string(const std::string& str, std::vector<std::string>& val) { 
    std::stringstream stream(str);
    while(stream.good()){
        std::string phrase;
        stream >> phrase;
        val.emplace_back(std::move(phrase));
    }
}


