#pragma once

#include <boost/type_traits.hpp>
#include "complext.h"

namespace Gadgetron {
	template<class T> struct enable_operators : public boost::false_type{};
	template<> struct enable_operators<float> : public boost::true_type{};
	template<> struct enable_operators<Gadgetron::complext<float> > : public boost::true_type{};
	template<> struct enable_operators<double> : public boost::true_type{};
	template<> struct enable_operators<Gadgetron::complext<double> > : public boost::true_type{};
}
