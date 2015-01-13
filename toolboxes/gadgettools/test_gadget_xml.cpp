#include "GadgetXml.h"

#include <iostream>

int main(int argc, char** argv)
{
  GDEBUG_STREAM("GadgetXML Test Program" << std::endl);

  TiXmlDocument doc( "demo.xml" );
  doc.LoadFile();

  GadgetXMLNode n(&doc);

  std::vector<long> vals = n.get<long>(std::string("gadgetron.encoding.kspace.matrix_size.value"));

  GDEBUG_STREAM("Number of values: " << vals.size() << std::endl);
  for (unsigned int i = 0; i < vals.size(); i++) {
    GDEBUG_STREAM("   :" << vals[i] << std::endl);
  }


  //Let's add something to the document
  n.add(std::string("gadgetron.encoding.mysection.value"), 6.789);
  n.add(std::string("gadgetron.encoding.mysection.value"), 612);
  n.add(std::string("gadgetron.encoding.mysection.value"), 512);
  n.add(std::string("gadgetron.encoding.mysection.value"), vals);
  
  n.get_document()->Print();
  

  return 0;
}
