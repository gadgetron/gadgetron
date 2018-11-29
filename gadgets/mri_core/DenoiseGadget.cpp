//
// Created by dchansen on 6/19/18.
//

#include "DenoiseGadget.h"
#include "GadgetronTimer.h"
#include "non_local_means.h"
#include "non_local_bayes.h"

namespace Gadgetron {
    int DenoiseGadget::process(ACE_Message_Block *mb) {


        if (auto m1 = AsContainerMessage<ISMRMRD::ImageHeader>(mb)) {
            if (auto m2 = AsContainerMessage<hoNDArray<float>>(m1->cont())) {
                return this->process(m1, m2);
            } else if (auto m2 = AsContainerMessage<hoNDArray<std::complex<float>>>(m1->cont())) {
                return this->process(m1, m2);
            }
        }


        if (auto m1 = AsContainerMessage<IsmrmrdImageArray>(mb)) {
            return this->process(m1);
        }

            return this->next()->putq(mb);
    }


    template<class T>
    int DenoiseGadget::process(GadgetContainerMessage<ISMRMRD::ImageHeader> *header_msg,
                               GadgetContainerMessage<hoNDArray<T>> *image_msg) {

        auto &input = *image_msg->getObjectPtr();
        input = denoise_function(input);
        this->next()->putq(header_msg);

        return GADGET_OK;

    }

    int DenoiseGadget::process(GadgetContainerMessage<IsmrmrdImageArray> *image_array) {

        auto &input = image_array->getObjectPtr()->data_;
        input = denoise_function(input);
        this->next()->putq(image_array);

        return GADGET_OK;

    }

    template<class T>
    hoNDArray<T> DenoiseGadget::denoise_function(const Gadgetron::hoNDArray<T> & input) {

        if (denoiser == "non_local_bayes") {
            return Denoise::non_local_bayes(input, image_std, search_radius);
        } else if (denoiser == "non_local_means") {
            return Denoise::non_local_means(input, image_std, search_radius);
        } else {
            throw std::invalid_argument(std::string("DenoiseGadget: Unknown denoiser type: ") + std::string(denoiser));
        }

    }


    GADGET_FACTORY_DECLARE(DenoiseGadget)

}