//
// Created by dchansen on 6/19/18.
//
// Above comment adds no value - it's just noise. Git remembers who did what far better than comments ever will.

#include <GadgetronTimer.h>
#include "DenoiseGadget.h"
#include "non_local_means.h"
#include "non_local_bayes.h"

namespace Gadgetron {
    // In general: Whitespace. Namespace is indented 8 spaces - why?
        int DenoiseGadget::process(ACE_Message_Block* mb)
        {
          auto m1 = AsContainerMessage<ISMRMRD::ImageHeader>(mb);

          // There's an awful lot of explicit null checking going on here. I'm aware it's idiomatic gadgetron, but
          // it feels needlessly complicated and messy. What you would like to do it have a good look at the
          // continuation, and call the appropriate process function, based on the types deduced.

          // auto message = AsContainerMessage<ISMRMRD::ImageHeader>(mb);
          // message.call_appropriate_process_function(this);

          if (m1){
              if(auto m2 = AsContainerMessage<hoNDArray<float>>(m1->cont())){
                  return this->process(m1,m2);
              } else if (auto m2 = AsContainerMessage<hoNDArray<std::complex<float>>>(m1->cont())){
                  return this->process(m1,m2);
              }
          }

          // This seems entirely generic - the flag lives on the Gadget superclass, why does the subclasses have to
          // worry about implementing this feature? Should they not inherit solid default behaviour from Gadget?
          if (pass_on_undesired_data_){
              return this->next()->putq(mb);
          } else {
              GERROR("Denoise Gadget received wrong input data type"); // Newline?
              return -1;
          }
        }
        // Newline here? I prefer two, one is okay, none is too few. I like the five-space indent though. Very rebellious ;)
     template<class T> int DenoiseGadget::process(GadgetContainerMessage<ISMRMRD::ImageHeader>* header_msg, GadgetContainerMessage<hoNDArray<T>>* image_msg) {

        // I don't like how 'determining what is to be done' is mixed in with 'doing what is to be done'. I'd prefer something like:
        // auto denoise_function = prepare_denoise_function(denoiser, image_std, search_radius);
        //
        // auto &input = *image_msg->getObjectPtr();
        // input = denoise_function(input);
        //
        // this->next()->putq(header_msg);
        //
        // return GADGET_OK;

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