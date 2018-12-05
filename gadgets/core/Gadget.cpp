#include "Gadget.h"
#include <boost/make_shared.hpp>

namespace Gadgetron {
    std::string Gadget::get_string_value(const char *name) {

        std::string str_val;
        GadgetPropertyBase *p = find_property(name);
        if (!p) {
            GERROR("Property %s\n", name);
            throw std::runtime_error("Attempting to access non existent property on Gadget");
        }
        return std::string(p->string_value());
    }

}
