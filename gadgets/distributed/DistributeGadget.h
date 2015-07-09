#ifndef DISTRIBUTEGADGET_H
#define DISTRIBUTEGADGET_H

#include "Gadget.h"
#include "gadgetron_distributed_gadgets_export.h"

#include <complex>

namespace Gadgetron{

  class EXPORTDISTRIBUTEDGADGETS DistributeGadget : public BasicPropertyGadget
    {
    public:
      GADGET_DECLARE(DistributeGadget);
      
    protected:
      GADGET_PROPERTY(collector, std::string, "Name of collection Gadget", "Collect");

      virtual int process(ACE_Message_Block* m);
      virtual int process_config(ACE_Message_Block* m);
      virtual int node_index(ACE_Message_Block* m);

      const char* get_node_xml_config();
      
      Gadget* collect_gadget_;
      
    private:
      std::string node_xml_config_;
      
    };
}
#endif //DISTRIBUTEGADGET_H
