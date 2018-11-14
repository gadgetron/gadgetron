#pragma once

#include "ismrmrd.h"
namespace Gadgetron { namespace Core {

    class Gadget {
    public:
        Gadget(const ISMRMRD:IsmrmrdHeader& header){};

        virtual void process(Channel& in, Channel& out) = 0;





    };
}}