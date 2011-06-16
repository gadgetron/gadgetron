#ifndef GADGETXML_H
#define GADGETXML_H

#include <string>

#include "tinyxml.h"

inline TiXmlElement* AddClassToXML(TiXmlNode* anchor, const char* section, const char* type, const char* class_name, const char* dll, int slot = 0, const char* name = 0)
{
  TiXmlElement* pSectionElem = anchor->FirstChildElement(section);
  if (!pSectionElem) {
    pSectionElem = new TiXmlElement(section);
    anchor->LinkEndChild(pSectionElem);
  }
  
  TiXmlElement* classElement = new TiXmlElement( type );
  classElement->SetAttribute("class", class_name);
  classElement->SetAttribute("dll", dll);
  if (slot) {
    classElement->SetAttribute("slot", slot);
  }
  if (name) {
    classElement->SetAttribute("name", name);
  }
  pSectionElem->LinkEndChild(classElement);
  return classElement;
}

inline TiXmlElement* AddWriterToXML(TiXmlNode* anchor, const char* class_name, const char* dll, int slot)
{
  return AddClassToXML(anchor, "writers", "writer", class_name, dll, slot);
}

inline TiXmlElement* AddReaderToXML(TiXmlNode* anchor, const char* class_name, const char* dll, int slot)
{
  return AddClassToXML(anchor, "readers", "reader", class_name, dll, slot);
}

inline TiXmlElement* AddGadgetToXML(TiXmlNode* anchor, const char* name, const char* class_name, const char* dll)
{
  return AddClassToXML(anchor, "stream", "gadget", class_name, dll,0,name);
}

template <class T> inline TiXmlElement* AddPropertyToXMLElement(TiXmlNode* anchor, const char* name, T value) 
{
  TiXmlElement* parameterElement = new TiXmlElement( "property" );
  parameterElement->SetAttribute("name", name);
  parameterElement->SetAttribute("value", value);
  anchor->LinkEndChild(parameterElement);
  return parameterElement;
}

 
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

inline TiXmlElement* AddDoubleParameterToXML(TiXmlNode* anchor, const char* section, 
					     const char* name, double value) 
{
  
  TiXmlElement* pSectionElem = anchor->FirstChildElement(section);

  if (!pSectionElem) {
    pSectionElem = new TiXmlElement(section);
    anchor->LinkEndChild(pSectionElem);
  }

  TiXmlElement* parameterElement = new TiXmlElement( "parameter" );
  parameterElement->SetAttribute("name", name);
  parameterElement->SetDoubleAttribute("value", value);

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
    if ((child->Type() == TiXmlNode::TINYXML_ELEMENT) &&
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
