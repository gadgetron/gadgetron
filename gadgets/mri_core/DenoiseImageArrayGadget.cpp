//
// Created by dchansen on 6/20/18.
//
// Git blame <- More awesome than blame comments - CLion wasted 3 lines of everyone's life.

#include <mri_core_data.h>
#include <non_local_means.h>
#include <non_local_bayes.h>
#include "DenoiseImageArrayGadget.h"

namespace Gadgetron {
    // Namespace indented 5 spaces? You rebel! Again!
     int DenoiseImageArrayGadget::process( GadgetContainerMessage<IsmrmrdImageArray>* image_msg) {

        // This looks A LOT like code I've already seen. Header as well. What's with the duplication?
        // See the header for some serious comments on the duplication.
         auto& input = image_msg->getObjectPtr()->data_;

         if (denoiser == "non_local_bayes") {
             input = Denoise::non_local_bayes(input, image_std, search_radius);
         } else if (denoiser == "non_local_means") {
             input = Denoise::non_local_means(input, image_std, search_radius);
         } else {
             throw std::invalid_argument(std::string("DenoiseGadget: Unknown denoiser type: ") + std::string(denoiser));
         }


         this->next()->putq(image_msg);

         return GADGET_OK;

     }


    GADGET_FACTORY_DECLARE(DenoiseImageArrayGadget)

}

