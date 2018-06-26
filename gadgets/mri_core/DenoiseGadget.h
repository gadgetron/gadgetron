//
// Created by dchansen on 6/19/18.
//

#ifndef GADGETRON_DENOISEGADGET_H
#define GADGETRON_DENOISEGADGET_H


#include <ismrmrd/ismrmrd.h>
#include <string>
#include "Gadget.h"
#include "hoNDArray.h"
#include "gadgetron_mricore_export.h"

namespace Gadgetron {
    class EXPORTGADGETSMRICORE DenoiseGadget : public BasicPropertyGadget{


    public:
        GADGET_DECLARE(DenoiseGadget)


        GADGET_PROPERTY(image_std,float,"Standard deviation of the noise in the produced image",1);
        GADGET_PROPERTY(search_radius,int,"Standard deviation of the noise in the produced image",25);
        GADGET_PROPERTY(denoiser,std::string,"Type of denoiser - non_local_means or non_local_bayes","non_local_bayes");



    protected:
       int process(ACE_Message_Block* mb);
       template<class T> int process(GadgetContainerMessage<ISMRMRD::ImageHeader>*, GadgetContainerMessage<hoNDArray<T>>*);


    };

}


#endif //GADGETRON_DENOISEGADGET_H
