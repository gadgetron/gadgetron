#ifndef GADGETXML_H
#define GADGETXML_H

#include <string>

#include "tinyxml.h"

template <class T> inline TiXmlElement* AddParameterToXML(TiXmlNode* anchor, const char* section, 
							  const char* name, T value) 
{
  
  TiXmlElement* pSectionElem = anchor->FirstChildElement(section);

  if (!pSectionElem) {
    pSectionElem = new TiXmlElement(section);
    anchor->LinkEndChild(pSectionElem);
  }

  TiXmlElement* parameterElement = new TiXmlElement( "parameter" );
  parameterElement->SetAttribute("name", name);
  parameterElement->SetAttribute("value", value);

  pSectionElem->LinkEndChild(parameterElement);
  
  return parameterElement;
}

inline std::string GetStringParameterValueFromXML(TiXmlNode* anchor, const char* section, 
				     const char* name)
{
  TiXmlNode* child = 0;;
  TiXmlNode* parent = anchor->FirstChildElement(section);

  std::string ret("");
  while( (child = parent->IterateChildren( child )) ) {
    if ((child->Type() == TiXmlNode::ELEMENT) &&
	(ACE_OS::strncmp(child->ToElement()->Value(),"parameter",9) == 0)) 
      {
	if ((ACE_OS::strncmp(child->ToElement()->Attribute("name"),name,100) == 0)) {
	  ret = std::string(child->ToElement()->Attribute("value"));
	  break;
	}
      }
  }

  return ret;
}

inline int GetIntParameterValueFromXML(TiXmlNode* anchor, const char* section, 
				       const char* name)
{
  return ACE_OS::atoi(GetStringParameterValueFromXML(anchor, section, name).c_str());
}

inline double GetDoubleParameterValueFromXML(TiXmlNode* anchor, const char* section, 
					     const char* name)
{
  return ACE_OS::atof(GetStringParameterValueFromXML(anchor, section, name).c_str());
}

inline std::string XmlToString(TiXmlNode& base)
{
  TiXmlPrinter printer;
  base.Accept( &printer );
  std::string xmltext = printer.CStr();
  return xmltext;
}


#endif //GADGETXML_H
