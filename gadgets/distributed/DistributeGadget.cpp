#include "DistributeGadget.h"
#include "GadgetStreamInterface.h"
#include "gadgetron_xml.h"
#include "CloudBus.h"

namespace Gadgetron{

  const char* DistributeGadget::get_node_xml_config()
  {
    return node_xml_config_.c_str();
  }
  
  int DistributeGadget::process(ACE_Message_Block* m)
  {
    //It is enough to put the first one, since they are linked
    if (this->next()->putq(m) == -1) {
      m->release();
      GERROR("DistributeGadget::process, passing data on to next gadget");
      return -1;
    }
    
    return 0;
  }

  int DistributeGadget::process_config(ACE_Message_Block* m)
  {

    //Grab the original XML conifguration
    std::string xml = controller_->get_xml_configuration();

    GadgetronXML::GadgetStreamConfiguration cfg;
    GadgetronXML::deserialize(xml.c_str(),cfg);

    //Delete Gadgets up to this Gadget
    std::vector<GadgetronXML::Gadget>::iterator it = cfg.gadget.begin();
    while ((it->name != std::string(this->module()->name())) && (it != cfg.gadget.end())) it++; it++;
    cfg.gadget.erase(cfg.gadget.begin(),it);

    //Delete Gadgets after collector
    it = cfg.gadget.begin();
    while ((it->name != collector.value()) && (it != cfg.gadget.end())) it++;
    cfg.gadget.erase(it,cfg.gadget.end());
    
    std::stringstream o;
    GadgetronXML::serialize(cfg,o);

    node_xml_config_ = o.str();

    Gadget* tmp = this;
    while (tmp->next()) {
      if (std::string(tmp->module()->name()) == collector.value()) break;
      tmp = dynamic_cast<Gadget*>(tmp->next());
    }

    collect_gadget_ = tmp;

    if (!collect_gadget_) {
      GERROR("Failed to locate collector Gadget with name %s\n", collector.value().c_str());
      return GADGET_FAIL;
    }

    return GADGET_OK;
  }
  
  int DistributeGadget::node_index(ACE_Message_Block* m)
  {
    return 0;
  }
  
  GADGET_FACTORY_DECLARE(DistributeGadget)
}


