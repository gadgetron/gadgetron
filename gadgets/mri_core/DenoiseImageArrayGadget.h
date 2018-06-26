//
// Created by dchansen on 6/19/18.
//

#ifndef GADGETRON_DENOISEGADGET_H
#define GADGETRON_DENOISEGADGET_H


#include <ismrmrd/ismrmrd.h>
#include "Gadget.h"
#include "hoNDArray.h"
#include "gadgetron_mricore_export.h"

namespace Gadgetron {
    // What is up with the code duplication here. We need to have a fairly serious chat about why you need two different
    // classes to do the same thing. Why on earth doesn't the DenoiseGadget know how to denoise Image Arrays?

    class EXPORTGADGETSMRICORE DenoiseImageArrayGadget : public Gadget1<IsmrmrdImageArray> {

        using super = Gadget1<IsmrmrdImageArray>;

    public:
        GADGET_DECLARE(DenoiseImageArrayGadget)


        GADGET_PROPERTY(image_std,float,"Standard deviation of the noise in the produced image",1);
        GADGET_PROPERTY(search_radius,int,"Standard deviation of the noise in the produced image",25);
        GADGET_PROPERTY(denoiser,std::string,"Type of denoiser - non_local_means or non_local_bayes","non_local_bayes");


    protected:
        virtual int process(GadgetContainerMessage<IsmrmrdImageArray>*) override ;


    };

}


#endif //GADGETRON_DENOISEGADGET_H
