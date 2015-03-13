#include <iostream>
#include <string>
#include <fstream>
#include "gadgetron_xml.h"

int main(int argc, char** argv)
{
  if (argc != 2) {
    std::cout << "Two input arguments reguired" << std::endl;
  }

  GadgetronXML::GadgetStreamConfiguration cfg;
  std::ifstream t(argv[1]);
  std::string gcfg_text((std::istreambuf_iterator<char>(t)),
			std::istreambuf_iterator<char>());
  
  GadgetronXML::deserialize(gcfg_text.c_str(), cfg);


  std::cout << "Number of Gadgets in configuration: " << cfg.gadget.size() << std::endl;
  return 0;
}
