#ifndef GADGETRON_CONFIGPARSERS_H
#define GADGETRON_CONFIGPARSERS_H

#include <iostream>
#include "AST.h"

namespace Gadgetron {
    AST::Chain parse_stream_configuration(std::istream &stream);
}

#endif //GADGETRON_CONFIGPARSERS_H
