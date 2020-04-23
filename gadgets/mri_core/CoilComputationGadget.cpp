#include "CoilComputationGadget.h"
#include "gadgetron/mri_core_coil_map_estimation.h"

/*
*       CoilComputationGadget.cpp
*       Author: Johannes Mayer
*/



namespace Gadgetron
{
    CoilComputationGadget::CoilComputationGadget()
    {
    }

    CoilComputationGadget::~CoilComputationGadget()
    {
    }


    int CoilComputationGadget::process(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1, GadgetContainerMessage<ArrayType>* m2)
    {
        ArrayType* input_array = m2->getObjectPtr();

        std::vector<size_t> dims;
        input_array->get_dimensions(dims);

        if (input_array->get_number_of_dimensions() == 4)
        {
            ArrayType csm;

//            coil_map_Inati<ValueType>(*input_array, csm, ks_.value(), kz_.value(), power_.value());
            coil_map_Inati<ValueType>(*input_array, csm);
            *m2->getObjectPtr() = csm;
        }
        else
        {
            GERROR_STREAM("CoilComputationGadget, no 4D data was passed ... ");
            return GADGET_FAIL;
        }

        if (this->next()->putq(m1) < 0)
        {
            GERROR_STREAM("CoilComputationGadget, failed to pass images to next gadget ... ");
            return GADGET_FAIL;
        }

        return GADGET_OK;
    }

    GADGET_FACTORY_DECLARE(CoilComputationGadget)
}
