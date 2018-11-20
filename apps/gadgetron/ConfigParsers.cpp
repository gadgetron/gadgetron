#include "ConfigParsers.h"

#include "gadgetron_xml.h"

#include <iostream>
#include <memory>

namespace {

    ConfigParser EXPORTGADGETBASE select_config_parser(pugi::xml_document &raw_config) {

    }

    class ConfigParser {
        virtual bool accepts(pugi::xml_document &config);
        virtual GadgetronXML::GadgetStreamConfiguration parse(pugi::xml_document &config);
    };


    class Legacy : public ConfigParser {

    };

    class V2 : public ConfigParser {

    };
}

namespace Gadgetron {

    GadgetronXML::GadgetStreamConfiguration parse_stream_configuration(std::istream &stream) {

        pugi::xml_document doc;
        pugi::xml_parse_result result = doc.load(stream);

        if (result.status != pugi::status_ok) {
            GERROR("Loading config file failed with following error: %s (%d)\n", result.description(), result.status);
            throw std::runtime_error(result.description());
        }

        auto parser = select_config_parser(config)
        return parser.parse(config);
    }
}