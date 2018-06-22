//
// Created by dchansen on 6/20/18.
//
#include <mri_core_data.h>
#include <non_local_means.h>
#include <non_local_bayes.h>
#include "DenoiseImageArrayGadget.h"

namespace Gadgetron {
     int DenoiseImageArrayGadget::process( GadgetContainerMessage<IsmrmrdImageArray>* image_msg) {


         auto& input = image_msg->getObjectPtr()->data_;
         input = Denoise::non_local_bayes(input,image_std,search_radius);

         this->next()->putq(image_msg);

         return GADGET_OK;

     }


    GADGET_FACTORY_DECLARE(DenoiseImageArrayGadget)

}

