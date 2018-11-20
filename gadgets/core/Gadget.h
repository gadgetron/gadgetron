#pragma once

#include <ismrmrd/xml.h>

#include "Channel.h"
namespace Gadgetron { namespace Core {

    template<class T> class Gadget {
    public:
        Gadget(const ISMRMRD::IsmrmrdHeader& header){};

        virtual void process(InputChannel<T>& in, OutputChannel& out) = 0;





    };
}}
