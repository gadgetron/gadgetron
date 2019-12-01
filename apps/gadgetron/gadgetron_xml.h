#ifndef GADGETRON_XML_H
#define GADGETRON_XML_H

#include <string>
#include <vector>
#include <stdexcept>
#include "gadgetbase_export.h"

namespace GadgetronXML
{
  template <typename T> class Optional
  {
  public:
    Optional()
      : present_(false)
    {

    }

    Optional(const T&v) {
      present_ = true;
      value_ = v;      
    }

    const Optional& operator=(const T& v) {
      present_ = true;
      value_ = v;
      return *this;
    }

    const T* operator->() const {
      return &value_;
    }

    const T& operator*() const {
      return value_;
    }

    operator bool() const {
      return present_;
    }

    bool is_present() const {
      return present_;
    }

    T& get() {
      if (!present_) {
	throw std::runtime_error("Access optional value, which has not been set");
      }
      return value_;
    }
    
    T& operator()() {
      return get();
    }

    void set(const T& v) {
      present_ = true;
      value_ = v;
    }

  protected:
    bool present_;
    T value_;

  }; 


  struct GadgetronParameter
  {
    std::string name;
    std::string value;
  };

  
  struct CloudBus
  {
    std::string relayAddress;
    unsigned int port;
    Optional<std::string> lbEndpoint;
  };


  struct ReST
  {
    unsigned int port;
  };
  
  struct GadgetronConfiguration
  {
    std::string port;
    std::vector<GadgetronParameter> globalGadgetParameter;
    Optional<CloudBus> cloudBus;
    Optional<ReST> rest;
  };

  void EXPORTGADGETBASE deserialize(std::istream& stream, GadgetronConfiguration& h);

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

  void EXPORTGADGETBASE deserialize(const std::string& stream, GadgetStreamConfiguration& cfg);
  void EXPORTGADGETBASE serialize(const GadgetStreamConfiguration& cfg, std::ostream& o);

};

#endif //GADGETRON_XML_H


