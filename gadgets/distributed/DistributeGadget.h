#ifndef DISTRIBUTEGADGET_H
#define DISTRIBUTEGADGET_H

#include "Gadget.h"
#include "gadgetron_distributed_gadgets_export.h"
#include "GadgetronConnector.h"

#include <complex>

namespace Gadgetron{

  class DistributeGadget;
  
  class DistributionConnector : public GadgetronConnector
  {
    
  public:
    DistributionConnector(DistributeGadget* g);
    virtual int process(size_t messageid, ACE_Message_Block* mb);

  protected:
    DistributeGadget* distribute_gadget_;
  };

  
  class EXPORTDISTRIBUTEDGADGETS DistributeGadget : public BasicPropertyGadget
    {
    public:
      GADGET_DECLARE(DistributeGadget);
      virtual int collector_putq(ACE_Message_Block* m);
      
    protected:
      GADGET_PROPERTY(collector, std::string, "Name of collection Gadget", "Collect");
      GADGET_PROPERTY(single_package_mode, bool, "Indicates that only one package is sent to each node", false);
      GADGET_PROPERTY(use_this_node_for_compute, bool, "This node can also be used for computation", true);
      
      virtual int process(ACE_Message_Block* m);
      virtual int process_config(ACE_Message_Block* m);
      virtual int close(unsigned long flags);
      /**
	 Returns the index of the node to process this package on. 

	 '0' means local
	 positive means networked node. 
	 negative is an error.
       */
      virtual int node_index(ACE_Message_Block* m);

      /**
	 Returns the message ID associated with this message
      */
      virtual int message_id(ACE_Message_Block* m)
      {
	return 0; //This is an invalid ID. 
      }

      const char* get_node_xml_config();
      
      Gadget* collect_gadget_;
      
    private:
      std::string node_xml_config_;
      std::string node_parameters_;
      std::map<int,GadgetronConnector*> node_map_;
      
    };
}
#endif //DISTRIBUTEGADGET_H
