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
    class EXPORTGADGETSMRICORE DenoiseImageArrayGadget : public Gadget1<IsmrmrdImageArray> {

        using super = Gadget1<IsmrmrdImageArray>;

    public:
        GADGET_DECLARE(DenoiseImageArrayGadget)


        GADGET_PROPERTY(image_std,float,"Standard deviation of the noise in the produced image",1);
        GADGET_PROPERTY(search_radius,int,"Standard deviation of the noise in the produced image",5);


    protected:
        virtual int process(GadgetContainerMessage<IsmrmrdImageArray>*) override ;


    };

}


#endif //GADGETRON_DENOISEGADGET_H
