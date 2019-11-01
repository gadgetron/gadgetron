#pragma once

#include <sstream>
#include <string>
#include <vector>
#include <boost/filesystem/path.hpp>
namespace Gadgetron::Core::IO {

	template<class T>
	std::enable_if_t<std::is_integral<T>::value> from_string(const std::string&, T& val);

	template<>
	void from_string<long long>(const std::string&, long long& val);

	template<class T>
	std::enable_if_t<std::is_integral<T>::value> from_string<T>(const std::string& str, T& val) {
		long long tmp;
		from_string(str, tmp);
		if (tmp < std::numeric_limits<T>::min() || tmp > std::numeric_limits<T>::max()) throw std::runtime_error("Value does not fit in desired type");
		val = static_cast<T>(tmp);
	}


	template<class T>
	std::enable_if_t<std::is_floating_point<T>::value> from_string(const std::string& str, T& val);

	template<>
	void from_string<double>(const std::string&, double& val);

	template<class T>
	std::enable_if_t<std::is_floating_point<T>::value> from_string<T>(const std::string& str, T& val) {
		double tmp;
		from_string(str, tmp);
		if (tmp < std::numeric_limits<T>::min() || tmp > std::numeric_limits<T>::max()) throw std::runtime_error("Value does not fit in desired type");
		val = static_cast<T>(tmp);
	}


	template<class T>
	std::enable_if_t<std::is_floating_point<T>::value>
		from_string(const std::string& str, std::vector<T>& values);

	template<class T>
	std::enable_if_t<std::is_integral<T>::value>
		from_string(const std::string& str, std::vector<T>& values);


	template<>
	void from_string<double>(const std::string& str, std::vector<double>& values);

	template<>
	void from_string<long long>(const std::string& str, std::vector<long long>& values);


	template<class T>
	std::enable_if_t<std::is_floating_point<T>::value>
		from_string<T>(const std::string& str, std::vector<T>& values) {
		std::vector<double> tmp;
		from_string(str, tmp);
		values.clear();
		auto check_conversion = [](double val) {
			if (val < std::numeric_limits<T>::min() || val > std::numeric_limits<T>::max()) throw std::runtime_error("Value does not fit in desired type");
			return static_cast<T>(val);
		};
		for (auto& val : tmp) values.push_back(check_conversion(val));
		
	}
	template<class T>
	std::enable_if_t<std::is_integral<T>::value>
		from_string<T>(const std::string& str, std::vector<T>& values) {
		std::vector<long long> tmp;
		from_string(str, tmp);
		values.clear();
		auto check_conversion = [](long long val) {
			if (val < std::numeric_limits<T>::min() || val > std::numeric_limits<T>::max()) throw std::runtime_error("Value does not fit in desired type");
			return static_cast<T>(val);
		};
		for (auto& val : tmp) values.push_back(check_conversion(val));
		
	}

    void from_string(const std::string&, boost::filesystem::path&);
	void from_string(const std::string&, bool&);
	void from_string(const std::string&, std::vector<bool>&);

    template <class T> T from_string(const std::string& str) {
        T val;
        from_string(str, val);
        return val;
    }
}
