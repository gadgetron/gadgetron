#pragma once

#include <bitset>
#include <boost/filesystem/path.hpp>
#include <sstream>
#include <string>
#include <vector>
#include <set>

namespace Gadgetron::Core::IO {

    void from_string(const std::string&, long long& val);

    template <class T> std::enable_if_t<std::is_integral<T>::value> from_string(const std::string& str, T& val) {
        long long tmp;
        from_string(str, tmp);
        if (tmp < std::numeric_limits<T>::min() || tmp > std::numeric_limits<T>::max())
            throw std::runtime_error("Value does not fit in desired type");
        val = static_cast<T>(tmp);
    }

  void from_string(const std::string&, double&);

template <class T> auto from_string(const std::string& str, T& val) -> std::enable_if_t<std::is_floating_point<T>::value> {
        double tmp;
        from_string(str, tmp);
        if (tmp < std::numeric_limits<T>::lowest() || tmp > std::numeric_limits<T>::max())
            throw std::runtime_error("Value does not fit in desired type");
        val = static_cast<T>(tmp);
    }

    void from_string(const std::string& str, std::vector<double>& values);

    void from_string(const std::string& str, std::vector<long long>& values);

    template <class T>
    std::enable_if_t<std::is_floating_point<T>::value> from_string(const std::string& str, std::vector<T>& values) {
        std::vector<double> tmp;
        from_string(str, tmp);
        values.clear();
        auto check_conversion = [](double val) {
            if (val < std::numeric_limits<T>::min() || val > std::numeric_limits<T>::max())
                throw std::runtime_error("Value does not fit in desired type");
            return static_cast<T>(val);
        };
        for (auto& val : tmp)
            values.push_back(check_conversion(val));
    }
    template <class T>
    std::enable_if_t<std::is_integral<T>::value> from_string(const std::string& str, std::vector<T>& values) {
        std::vector<long long> tmp;
        from_string(str, tmp);
        values.clear();
        auto check_conversion = [](long long val) {
            if (val < std::numeric_limits<T>::min() || val > std::numeric_limits<T>::max())
                throw std::runtime_error("Value does not fit in desired type");
            return static_cast<T>(val);
        };
        for (auto& val : tmp)
            values.push_back(check_conversion(val));
    }

    void from_string(const std::string&, boost::filesystem::path&);
    void from_string(const std::string&, bool&);
    void from_string(const std::string&, std::vector<bool>&);

    void from_string(const std::string&, std::vector<std::string>&);
    void from_string(const std::string&, std::set<std::string>&);

    template<class T> T from_string(const std::string& str);

    template <size_t N> void from_string(const std::string& str, std::bitset<N>& bset) {
        bset = std::bitset<N>(from_string<unsigned int>(str));
    }

    template <class T> T from_string(const std::string& str) {
        T val;
        from_string(str, val);
        return val;
    }

}
