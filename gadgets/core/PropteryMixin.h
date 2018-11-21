#pragma once

#include <unordered_map>
#include <sstream>

namespace Gadgetron::Core {

    class PropertyMixin {

    protected:

        template<class KEY_VALUE_STORAGE>
        PropertyMixin(const KEY_VALUE_STORAGE &pairs) : properties(pairs.begin(), pairs.end()) {

        }

        template<class T>
        T get_property(const std::string &name, T default_value, const std::string &description) {
            if (properties.count(name)) {
                T val;
                std::stringstream stream(properties.at(name));
                stream >> val;
                return val;
            } else {
                return default_value;
            }
        }

    private:

        const std::unordered_map<std::string, std::string> properties;

    };
}

#define GADGET_PROPERTY(TYPE, NAME, DEFAULT, DESCRIPTION) const TYPE NAME = this->get_property<TYPE>(NAME,DEFAULT,DESCRIPTION)
