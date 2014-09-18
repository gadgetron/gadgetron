#include "gadgetron_xml.h"
#include "pugixml.hpp"
#include <stdexcept>
#include <cstdlib>

namespace GadgetronXML
{

  void deserialize(const char* xml_config, GadgetronConfiguration& h)
  {
    pugi::xml_document doc;
    pugi::xml_parse_result result = doc.load(xml_config);
    pugi::xml_node root = doc.child("gadgetronConfiguration");

    if (!root) {
      throw std::runtime_error("gadgetronConfiguration element not found in configuration file");
    }
    
    pugi::xml_node port = root.child("port");
    if (!port) {
      throw std::runtime_error("Port not found in Gadgetron configuration");
    }

    h.port = port.child_value();

    pugi::xml_node p = root.child("globalGadgetParameter");
    while (p) {
      GadgetronParameter pp;
      pp.name = p.child_value("name");
      pp.value = p.child_value("value");
      h.globalGadgetParameter.push_back(pp);
      p = p.next_sibling("globalGadgetParameter");
    }
  }

  void deserialize(const char* xml_config, GadgetStreamConfiguration& cfg)
  {
    pugi::xml_document doc;
    pugi::xml_parse_result result = doc.load(xml_config);
    pugi::xml_node root = doc.child("gadgetronStreamConfiguration");

    if (!root) {
      throw std::runtime_error("gadgetronStreamConfiguration element not found in configuration file");
    }

    pugi::xml_node reader = root.child("reader");
    while (reader) {
      Reader r;
      r.slot = static_cast<unsigned short>(std::atoi(reader.child_value("slot")));
      r.dll = reader.child_value("dll");
      r.classname = reader.child_value("classname");
      cfg.reader.push_back(r);
      reader = reader.next_sibling("reader");
    }
    
    pugi::xml_node writer = root.child("writer");
    while (writer) {
      Writer w;
      w.slot = static_cast<unsigned short>(std::atoi(writer.child_value("slot")));
      w.dll = writer.child_value("dll");
      w.classname = writer.child_value("classname");
      cfg.writer.push_back(w);
      writer = writer.next_sibling("writer");
    }

    pugi::xml_node gadget = root.child("gadget");
    while (gadget) {
      Gadget g;
      g.name = gadget.child_value("name");
      g.dll = gadget.child_value("dll");
      g.classname = gadget.child_value("classname");
      
      pugi::xml_node property = gadget.child("property");
      while (property) {
	GadgetronParameter p;
	p.name = property.child_value("name");
	p.value = property.child_value("value");
	g.property.push_back(p);
	property = property.next_sibling("property");
      }

      cfg.gadget.push_back(g);
      gadget = gadget.next_sibling("gadget");
    }
  }
}
