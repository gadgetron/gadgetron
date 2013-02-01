#include <octave/oct.h>
#include <iostream>

#include "pugixml.hpp"
     
DEFUN_DLD (XMLGetXPath, args, nargout,
	   "XMLGetXPath: Returns the text contents of the xml node with the given XPATH")
{
  int nargin = args.length ();

  octave_value retval;
     
  if (nargin != 2) {
    print_usage(); 
  } else {
    std::string xml(args(0).string_value());
    std::string xpath(args(1).string_value());

    pugi::xml_document doc;
    
    pugi::xml_parse_result result = doc.load_buffer_inplace(const_cast<char*>(xml.c_str()), xml.length());

    if (!result) {
      std::cout << "XML parsed with errors." << std::endl;
      std::cout << "Error description: " << result.description() << std::endl;
      return retval;
    }

    pugi::xpath_node target_node = doc.select_single_node(xpath.c_str());

    retval = octave_value(target_node.node().child_value());
  }

  return retval;
}
