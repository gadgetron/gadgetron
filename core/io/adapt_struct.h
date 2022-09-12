#pragma once
#include <boost/hana/adapt_struct.hpp>

#define GADGETRON_ADAPT_STRUCT(STRUCTNAME, ...) \
namespace boost {\
	namespace hana { \
		template <> \
		struct accessors_impl<STRUCTNAME> { \
			static BOOST_HANA_CONSTEXPR_LAMBDA auto apply() { \
				return make_tuple( __VA_ARGS__ ); \
			} \
		}; \
	}\
}

#define GADGETRON_ACCESS_ELEMENT(ELEMENT) (make_pair(BOOST_HANA_STRING(#ELEMENT), [](auto&& p) -> auto& { return p.ELEMENT;}))
