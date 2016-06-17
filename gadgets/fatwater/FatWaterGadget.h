#ifndef FATWATERGADGET_H
#define FATWATERGADGET_H

#include "Gadget.h"
#include "hoNDArray.h"
#include "gadgetron_fatwater_export.h"

#include "mri_core_data.h"

namespace Gadgetron{

  class EXPORTGADGETFATWATER FatWaterGadget : 
  public Gadget1<IsmrmrdImageArray>
    {
    public:
      GADGET_DECLARE(FatWaterGadget)
      FatWaterGadget();
	
    protected:
      virtual int process(GadgetContainerMessage<IsmrmrdImageArray>* m1);
      virtual int process_config(ACE_Message_Block* mb);
      

    private:
      std::vector<float> echoTimes_;
      float fieldStrength_;
      
      
    };
}
#endif //FATWATERGADGET_H
