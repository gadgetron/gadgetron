#pragma once
#include "GadgetContainerMessage.h"
#include "GadgetronExport.h"
#include "Gadget.h"

#include <map>

namespace Gadgetron {

    enum GadgetronMessageID {
        GADGET_MESSAGE_INT_ID_MIN = 0,
        GADGET_MESSAGE_CONFIG_FILE = 1,
        GADGET_MESSAGE_CONFIG_SCRIPT = 2,
        GADGET_MESSAGE_PARAMETER_SCRIPT = 3,
        GADGET_MESSAGE_CLOSE = 4,
        GADGET_MESSAGE_TEXT = 5,
        GADGET_MESSAGE_INT_ID_MAX = 999
    };

    struct GadgetMessageIdentifier {
        uint16_t id;
    };

}
