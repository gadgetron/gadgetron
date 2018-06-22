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
    class EXPORTGADGETSMRICORE DenoiseGadget : public Gadget2<ISMRMRD::ImageHeader, hoNDArray<float>> {

        using super = Gadget2<ISMRMRD::ImageHeader, hoNDArray<float>>;

    public:
        GADGET_DECLARE(DenoiseGadget)


        GADGET_PROPERTY(image_std,float,"Standard deviation of the noise in the produced image",1);
        GADGET_PROPERTY(search_radius,int,"Standard deviation of the noise in the produced image",25);


    protected:
        virtual int process(GadgetContainerMessage<ISMRMRD::ImageHeader>*, GadgetContainerMessage<hoNDArray<float>>*) override ;


    };

}


#endif //GADGETRON_DENOISEGADGET_H
