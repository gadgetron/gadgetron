//
// Created by dchansen on 6/19/18.
//

#include <GadgetronTimer.h>
#include "DenoiseGadget.h"
#include "non_local_means.h"
#include "non_local_bayes.h"

namespace Gadgetron {
        int DenoiseGadget::process(ACE_Message_Block* mb)
        {

          auto m1 = AsContainerMessage<ISMRMRD::ImageHeader>(mb);

          if (m1){
              if(auto m2 = AsContainerMessage<hoNDArray<float>>(m1->cont())){
                  return this->process(m1,m2);
              } else if (auto m2 = AsContainerMessage<hoNDArray<std::complex<float>>>(m1->cont())){
                  return this->process(m1,m2);
              }
          }

          if (pass_on_undesired_data_){
              return this->next()->putq(mb);
          } else {
              GERROR("Denoise Gadget received wrong input data type");
              return -1;
          }


        }
     template<class T> int DenoiseGadget::process(GadgetContainerMessage<ISMRMRD::ImageHeader>* header_msg, GadgetContainerMessage<hoNDArray<T>>* image_msg) {


         auto& input = *image_msg->getObjectPtr();
         if (denoiser == "non_local_bayes") {
             input = Denoise::non_local_bayes(input, image_std, search_radius);
         } else if (denoiser == "non_local_means") {
             input = Denoise::non_local_means(input, image_std, search_radius);
         } else {
             throw std::invalid_argument(std::string("DenoiseGadget: Unknown denoiser type: ") + std::string(denoiser));
         }


         this->next()->putq(header_msg);

         return GADGET_OK;

     }


    GADGET_FACTORY_DECLARE(DenoiseGadget)

}