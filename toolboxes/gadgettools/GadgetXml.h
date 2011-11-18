#ifndef GADGETXML_H
#define GADGETXML_H

#include <string>
#include <sstream>
#include <vector>

#include "tinyxml.h"

#include "ace/OS_NS_stdlib.h"
#include "ace/OS_NS_string.h"
#include "ace/OS_NS_stdio.h"
#include "ace/OS_NS_netdb.h"

//#include <boost/algorithm/string.hpp>

class GadgetXMLNode
{
 public:
  GadgetXMLNode(TiXmlNode* anchor) 
    : anchor_(anchor)
    { }

  template<typename T> std::vector<T> get(std::string name);

  GadgetXMLNode add(const std::string name) {
    //We have to make a copy of this string, the strtok function will eat it up
    char* tmp = new char[name.size()+1];
    memcpy(tmp,name.c_str(),name.size());
	tmp[name.size()] = 0;
    
    TiXmlHandle h( anchor_ );
    char* lname = strtok(tmp,".");
    TiXmlElement* child = 0;
    while (lname) {
      char* next = strtok(NULL,"."); //Get next token
      child = h.FirstChildElement(lname).ToElement();
      if (!child) {
	child = h.ToNode()->LinkEndChild(new TiXmlElement(lname))->ToElement();
      } else if (!next) {
	//We are on the last level, and we are finding the path, so there must tags with this name.
	//Let's add a sibling.
	child = child->Parent()->LinkEndChild(new TiXmlElement(lname))->ToElement();
      }
      lname = next;
      h = TiXmlHandle(child);
    }
    delete [] tmp; tmp = 0x0;
    return GadgetXMLNode(child);
  }


  template<typename T> GadgetXMLNode add(const std::string name, T contents) {
    GadgetXMLNode n = add(name);
    std::stringstream ss; ss << contents;
    n.get_node()->LinkEndChild(new TiXmlText(ss.str().c_str()));
    return n;
  }

  template<typename T> GadgetXMLNode add(std::string n, std::vector<T> c) {
    GadgetXMLNode ret(NULL);
    for (unsigned int i = 0; i < c.size(); i++) {
      ret = add<T>(n,c[i]);
    }
    return ret;
  }
  
  template<typename T> T contents();

  TiXmlNode* get_node() {
    return anchor_;
  }

  TiXmlDocument* get_document() {
    if (!anchor_) {
      return 0;
    }

    return anchor_->GetDocument();
  }

 protected:
  TiXmlNode* anchor_;

};


template<> inline std::string GadgetXMLNode::contents() {
  std::string ret_val("");
  TiXmlNode* e = anchor_->FirstChild();
  if (e && (e->Type() == TiXmlNode::TINYXML_TEXT)) {
    return std::string(e->Value());
  }
  return ret_val;
}

template<> inline long int GadgetXMLNode::contents() {
  std::string s = contents<std::string>();
  return atoi(s.c_str());
}

template<> inline double GadgetXMLNode::contents() {
  std::string s = contents<std::string>();
  return atof(s.c_str());
}

template<typename T> inline std::vector<T> contents(std::vector<GadgetXMLNode>& v)
{
  std::vector<T> ret(v.size());
  for (unsigned int i = 0; i < v.size(); i++) {
    ret[i] = v[i].contents<T>();
  }
  return ret;
}

template<> inline std::vector<GadgetXMLNode> GadgetXMLNode::get<GadgetXMLNode>(std::string name) {
  std::vector<GadgetXMLNode> ret;
  
  //We have to make a copy of this string, the strtok function will eat it up
  char* tmp = new char[name.size()+1];
  memcpy(tmp,name.c_str(),name.size());
  tmp[name.size()] = '\0';

  TiXmlHandle h( anchor_ );
  char* lname = strtok(tmp,".");
  TiXmlElement* child;
  while (lname && (child = h.FirstChildElement(lname).ToElement())) {
    h = TiXmlHandle(child);
    lname = strtok(NULL,"."); //Get next token
  }
  
  if (child) {
    ret.push_back(GadgetXMLNode(child));
    while ( (child = child->NextSiblingElement(child->Value())) ) {
      ret.push_back(GadgetXMLNode(child));
    }
  }
  
  delete [] tmp; tmp = 0x0;

  return ret;
}

template<> inline std::vector<long> GadgetXMLNode::get<long>(std::string name) {
  std::vector<GadgetXMLNode> vals = get<GadgetXMLNode>(name);
  return ::contents<long>(vals);
}

template<> inline std::vector<double> GadgetXMLNode::get<double>(std::string name) {
  std::vector<GadgetXMLNode> vals = get<GadgetXMLNode>(name);
  return ::contents<double>(vals);
}

template<> inline std::vector<std::string> GadgetXMLNode::get< std::string >(std::string name) {
  std::vector<GadgetXMLNode> vals = get<GadgetXMLNode>(name);
  return ::contents<std::string>(vals);
}


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
