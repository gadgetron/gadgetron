#include "Gadget.h"
#include <boost/make_shared.hpp>

namespace Gadgetron {
    std::string Gadget::get_string_value(const char *name) {
        std::string str_val;
        GadgetPropertyBase *p = find_property(name);
        if (!p) {
            GERROR("Property %s\n", name);
            throw std::runtime_error(
                    "Attempting to access non existent property on Gadget");
        }
        return std::string(p->string_value());
    }
    void Gadget::set_context(const Core::Context& input_context) { this->context = input_context;}

    LegacyGadgetNode::LegacyGadgetNode(std::unique_ptr<Gadget> gadget_ptr,
        const Core::Context& context,
        const std::unordered_map<std::string, std::string>& props)
        : gadget(std::move(gadget_ptr)) {

        gadget->set_context(context);
        for (auto &key_val : props) {
            gadget->set_parameter(key_val.first.c_str(), key_val.second.c_str());
        }

        gadget->process_config(context.header);
    }

    void LegacyGadgetNode::process(
            Core::GenericInputChannel& in,
            Core::OutputChannel& out) {

        gadget->next(std::make_shared<ChannelAdaptor>(out));

        for (auto message : in) {
            gadget->process(message.to_container_message());
        }
        gadget->close();
    }
}  // namespace Gadgetron
