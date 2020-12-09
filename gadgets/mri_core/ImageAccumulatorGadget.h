//
// Created by dchansen on 5/24/18.
//

#ifndef GADGETRON_IMAGEACCUMULATORGADGET_H
#define GADGETRON_IMAGEACCUMULATORGADGET_H

#include <mri_core_data.h>
#include "Gadget.h"
#include "hoNDArray.h"
#include "gadgetron_mricore_export.h"

namespace Gadgetron {
class EXPORTGADGETSMRICORE ImageAccumulatorGadget :  public Gadget1<IsmrmrdImageArray> {

    public:
        GADGET_DECLARE(ImageAccumulatorGadget)

        ImageAccumulatorGadget();

    protected:

        GADGET_PROPERTY(accumulate_dimension,std::string,"Dimension over which the images will be collected","contrast");
        GADGET_PROPERTY(combine_along,std::string,"Dimension used for stacking the images","N");
        virtual int process(GadgetContainerMessage<IsmrmrdImageArray>* m1);
        virtual int process_config(ACE_Message_Block* mb);

    private:

        template<class T> auto extract_value(T& val);

        IsmrmrdImageArray combine_images(std::vector<IsmrmrdImageArray>&);

        size_t encoding_spaces;
        std::vector<IsmrmrdImageArray> images;
        std::vector<uint16_t> required_values;
        std::set<uint16_t> seen_values;
};
}

#endif //GADGETRON_IMAGEACCUMULATORGADGET_H
