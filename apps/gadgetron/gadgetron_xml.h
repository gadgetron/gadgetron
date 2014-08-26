#ifndef GADGETRON_XML_H
#define GADGETRON_XML_H

#include <string>
#include <vector>

namespace GadgetronXML
{
  struct GadgetronParameter
  {
    std::string name;
    std::string value;
  };
  
  struct GadgetronConfiguration
  {
    std::string port;
    std::vector<GadgetronParameter> globalGadgetParameter;
  };
  
  void deserialize(const char* xml_config, GadgetronConfiguration& h);
  
  struct Reader
  {
    unsigned short slot;
    std::string dll;
    std::string classname;
  };

  typedef Reader Writer;
  
  struct Gadget
  {
    std::string name;
    std::string dll;
    std::string classname;
    std::vector<GadgetronParameter> property;
  };

  struct GadgetStreamConfiguration
  {
    std::vector<Reader> reader;
    std::vector<Writer> writer;
    std::vector<Gadget> gadget;
  };

  void deserialize(const char* xml, GadgetStreamConfiguration& cfg);

};

#endif //GADGETRON_XML_H


