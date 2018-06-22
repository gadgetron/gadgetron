//
// Created by dchansen on 6/19/18.
//

#include <GadgetronTimer.h>
#include "DenoiseGadget.h"
#include "non_local_means.h"
#include "non_local_bayes.h"

namespace Gadgetron {
     int DenoiseGadget::process(GadgetContainerMessage<ISMRMRD::ImageHeader>* header_msg, GadgetContainerMessage<hoNDArray<float>>* image_msg) {


         auto& input = *image_msg->getObjectPtr();
         GadgetronTimer("Non local Bayes");
         input = Denoise::non_local_bayes(input,image_std,search_radius);

         this->next()->putq(header_msg);

         return GADGET_OK;

     }


    GADGET_FACTORY_DECLARE(DenoiseGadget)

}