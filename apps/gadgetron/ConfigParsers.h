#ifndef GADGETRON_CONFIGPARSERS_H
#define GADGETRON_CONFIGPARSERS_H

#include "gadgetron_xml.h"
#include <iostream>


namespace Gadgetron {
    GadgetronXML::GadgetStreamConfiguration parse_stream_configuration(std::istream &stream);
}

#endif //GADGETRON_CONFIGPARSERS_H
