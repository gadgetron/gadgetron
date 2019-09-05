#pragma once

#include <unordered_map>
#include "io/from_string.h"

namespace Gadgetron::Core {

    using GadgetProperties = std::unordered_map<std::string,std::string>;

    class PropertyMixin {

    protected:

        template<class KEY_VALUE_STORAGE>
        explicit PropertyMixin(const KEY_VALUE_STORAGE &pairs) : properties(pairs.begin(), pairs.end()) {}

        template<class T>
        inline T get_property(const std::string &name, T default_value, const std::string &) {
            if (!properties.count(name)) return default_value;
            return IO::from_string<T>(properties.at(name));
        }

    private:
        const GadgetProperties properties;
    };

    template<>
    inline std::string PropertyMixin::get_property<std::string>(const std::string &name, std::string default_value, const std::string &) {
        if (!properties.count(name)) return default_value;
        return properties.at(name);
    }
}

#define NODE_PROPERTY(NAME,TYPE, DESCRIPTION, DEFAULT) const TYPE NAME = this->get_property<TYPE>(#NAME,DEFAULT,DESCRIPTION)
