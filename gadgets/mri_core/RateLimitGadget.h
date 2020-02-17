#ifndef RATELIMITGADGET_H
#define RATELIMITGADGET_H

#include "Gadget.h"
#include "hoNDArray.h"
#include "gadgetron_mricore_export.h"

#include <ismrmrd/ismrmrd.h>
#include <complex>
#include <chrono>

namespace Gadgetron{
  
  class EXPORTGADGETSMRICORE RateLimitGadget :
  public BasicPropertyGadget
    {
      
    public:
      GADGET_DECLARE(RateLimitGadget);
      
      RateLimitGadget();
      ~RateLimitGadget();
      
    protected:
      GADGET_PROPERTY(sleep_time_, int, "sleep_time", 0);

      virtual int process_config(ACE_Message_Block* mb);
        int process(ACE_Message_Block* mb);


        std::chrono::milliseconds sleep_time;

    };
}
#endif //ACCUMULATORGADGET_H
