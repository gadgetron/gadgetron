#include "gadgetron_xml.h"
#include "log.h"
#include "pugixml.hpp"
#include <stdexcept>
#include <cstdlib>
#include <iostream>

namespace GadgetronXML
{
  void deserialize(std::istream& stream, GadgetronConfiguration& h)
  {
    pugi::xml_document doc;
    pugi::xml_parse_result result = doc.load(stream);
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

    pugi::xml_node b = root.child("cloudBus");
    if (b) {
      CloudBus cb;
      cb.relayAddress = b.child_value("relayAddress");
      cb.port = static_cast<unsigned int>(std::atoi(b.child_value("port")));
      if ((cb.relayAddress.size() == 0) || (cb.port == 0))
      {
        throw std::runtime_error("Invalid CloudBus configuration.");
      }
      pugi::xml_node lbn = b.child("loadBalancedEndpoint");
      if (lbn) {
          cb.lbEndpoint = b.child_value("loadBalancedEndpoint");
      }
      h.cloudBus = cb;
    }

    pugi::xml_node r = root.child("rest");
    if (r) {
      ReST re;
      re.port = static_cast<unsigned int>(std::atoi(r.child_value("port")));
      if (re.port == 0) {
	throw std::runtime_error("Invalid ReST configuration.");
      }
      h.rest = re;
    }
  }

  void deserialize(const std::string& config_string, GadgetStreamConfiguration& cfg)
  {

    std::istringstream str(config_string);

    pugi::xml_document doc;
    pugi::xml_parse_result result = doc.load(str);

    if (result.status != pugi::status_ok) {
        GERROR("Loading config file failed with following error: %s (%d)\n", result.description(), result.status);
        throw std::invalid_argument(result.description());
    }

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


  void append_node(pugi::xml_node& n, const char* child, const std::string& v)
  {
    pugi::xml_node n2 = n.append_child(child);
    n2.append_child(pugi::node_pcdata).set_value(v.c_str());
  }

  std::string to_string_val(const unsigned short& v)
  {
    char buffer[256];
    sprintf(buffer,"%d",v);
    return std::string(buffer);
  }

  void serialize(const GadgetStreamConfiguration& cfg, std::ostream& o)
  {
    pugi::xml_document doc;
    pugi::xml_node root = doc.append_child();
    pugi::xml_node n1,n2,n3;
    pugi::xml_attribute a;

    root.set_name("gadgetronStreamConfiguration");

    a = root.append_attribute("xmlns");
    a.set_value("http://gadgetron.sf.net/gadgetron");

    a = root.append_attribute("xmlns:xsi");
    a.set_value("http://www.w3.org/2001/XMLSchema-instance");

    a = root.append_attribute("xmlns:xs");
    a.set_value("http://www.w3.org/2001/XMLSchema");

    a = root.append_attribute("xsi:schemaLocation");
    a.set_value("http://gadgetron.sf.net/gadgetron gadgetron.xsd");


    for (std::vector<Reader>::const_iterator it = cfg.reader.begin();
    it != cfg.reader.end(); it++)
    {
      n1 = root.append_child("reader");

      n2 = n1.append_child("slot");
      n2.append_child(pugi::node_pcdata).set_value(to_string_val(it->slot).c_str());

      n2 = n1.append_child("dll");
      n2.append_child(pugi::node_pcdata).set_value(it->dll.c_str());

      n2 = n1.append_child("classname");
      n2.append_child(pugi::node_pcdata).set_value(it->classname.c_str());
    }

    for (std::vector<Writer>::const_iterator it = cfg.writer.begin();
    it != cfg.writer.end(); it++)
    {
      n1 = root.append_child("writer");

      n2 = n1.append_child("slot");
      n2.append_child(pugi::node_pcdata).set_value(to_string_val(it->slot).c_str());

      n2 = n1.append_child("dll");
      n2.append_child(pugi::node_pcdata).set_value(it->dll.c_str());

      n2 = n1.append_child("classname");
      n2.append_child(pugi::node_pcdata).set_value(it->classname.c_str());
    }

    for (std::vector<Gadget>::const_iterator it = cfg.gadget.begin();
    it != cfg.gadget.end(); it++)
    {
      n1 = root.append_child("gadget");

      n2 = n1.append_child("name");
      n2.append_child(pugi::node_pcdata).set_value(it->name.c_str());

      n2 = n1.append_child("dll");
      n2.append_child(pugi::node_pcdata).set_value(it->dll.c_str());

      n2 = n1.append_child("classname");
      n2.append_child(pugi::node_pcdata).set_value(it->classname.c_str());

      for (std::vector<GadgetronParameter>::const_iterator it2 = it->property.begin();
      it2 != it->property.end(); it2++)
      {
        n2 = n1.append_child("property");
        n3 = n2.append_child("name");
        n3.append_child(pugi::node_pcdata).set_value(it2->name.c_str());
        n3 = n2.append_child("value");
        n3.append_child(pugi::node_pcdata).set_value(it2->value.c_str());
      }
    }

    doc.save(o);

  }

}
