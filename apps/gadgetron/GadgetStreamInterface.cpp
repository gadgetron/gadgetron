

#include "ace/Stream.h"
#include "ace/DLL.h"
#include "ace/DLL_Manager.h"

#include "gadgetron_home.h"
#include "gadgetron_xml.h"
#include "Gadget.h"
#include "EndGadget.h"

#include "GadgetStreamInterface.h"

namespace Gadgetron {

    GadgetStreamInterface::GadgetStreamInterface()
            : stream_configured_(false)
            , stream_(nullptr, nullptr, default_end_module())
    {
        gadgetron_home_ = get_gadgetron_home();
    }


    Gadget* GadgetStreamInterface::find_gadget(std::string gadget_name)
    {
        GadgetModule* gm = stream_.find(gadget_name.c_str());

        if (gm) {
            Gadget* g = dynamic_cast<Gadget*>(gm->writer());
            return g;
        } else {
            GDEBUG("Gadget with name %s not found! Returning null pointer\n", gadget_name.c_str());
        }
        return 0;
    }

    void GadgetStreamInterface::set_global_gadget_parameters(const std::map<std::string, std::string>& globalGadgetPara)
    {
        global_gadget_parameters_ = globalGadgetPara;
    }

    const GadgetronXML::GadgetStreamConfiguration& GadgetStreamInterface::get_stream_configuration()
    {
        return stream_configuration_;
    }


    GadgetModule *GadgetStreamInterface::create_gadget_module(const char* DLL, const char* gadget, const char* gadget_module_name)
    {

        Gadget* g = load_dll_component<Gadget>(DLL,gadget);

        if (!g) {
            GERROR("Failed to load gadget using factory\n");
            return 0;
        }

        g->set_controller(this);

        GadgetModule *module = 0;
        ACE_NEW_RETURN (module,
                        GadgetModule (gadget_module_name, g),
                        0);

        return module;
    }

    GadgetModule *GadgetStreamInterface::default_end_module(void)
    {
        Gadget *end_gadget = new EndGadget();
        end_gadget->set_controller(this);

        GadgetModule *end_module;
        ACE_NEW_RETURN(end_module, GadgetModule(ACE_TEXT("EndGadget"), end_gadget), nullptr);

        return end_module;
    }
}

