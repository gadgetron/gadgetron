#include "FatWaterGadget.h"

namespace Gadgetron{

FatWaterGadget::FatWaterGadget()
{

}

int FatWaterGadget::process( GadgetContainerMessage<IsmrmrdImageArray>* m1)
{

    GDEBUG("In FW process\n");
    
    //Grab a reference to the buffer containing the imaging data
    IsmrmrdImageArray & imagearr = *m1->getObjectPtr();

    //7D, fixed order [X, Y, Z, CHA, N, S, LOC]
    uint16_t X = imagearr.data_.get_size(0);
    uint16_t Y = imagearr.data_.get_size(1);
    uint16_t Z = imagearr.data_.get_size(2);
    uint16_t CHA = imagearr.data_.get_size(3);
    uint16_t N = imagearr.data_.get_size(4);
    uint16_t S = imagearr.data_.get_size(5);
    uint16_t LOC = imagearr.data_.get_size(6);

    //Pass the image down the chain
    if (this->next()->putq(m1) < 0) {
        return GADGET_FAIL;
    }

    return GADGET_OK;  

}

GADGET_FACTORY_DECLARE(FatWaterGadget)
}
