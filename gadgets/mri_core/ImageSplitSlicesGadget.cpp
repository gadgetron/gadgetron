#include "ImageSplitSlicesGadget.h"

namespace Gadgetron{

ImageSplitSlicesGadget::ImageSplitSlicesGadget()
{

}


int ImageSplitSlicesGadget::process(GadgetContainerMessage<ISMRMRD::ImageHeader> *m1, GadgetContainerMessage<hoNDArray<std::complex<float> > > *m2)
{

    //Grab a reference to the Image Header and some dimensions
    auto hdr = *m1->getObjectPtr();
    auto X = hdr.matrix_size[0];
    auto Y = hdr.matrix_size[1];
    auto Z = hdr.matrix_size[2];
    auto CHA = hdr.channels;

    //Grab a reference to the input img
    auto img_in  = *m2->getObjectPtr();

    //Get a handle to the meta attributes (if any)
    GadgetContainerMessage<ISMRMRD::MetaContainer>* m3 = AsContainerMessage<ISMRMRD::MetaContainer>(m2->cont());

    //Compute the slice thickness
    float thick = hdr.field_of_view[2] / hdr.matrix_size[2];

    //Size in bytes of a single slice single channel
    auto slice_size = X*Y*sizeof(std::complex<float>);

    //Loop over slices
    for (uint16_t sl=0; sl < Z; sl++) {

        //Copy the image header
        GadgetContainerMessage<ISMRMRD::ImageHeader>* cm1 = new GadgetContainerMessage<ISMRMRD::ImageHeader>();
        *cm1->getObjectPtr() = hdr;

        //Set the number of slices, slice index and position
        cm1->getObjectPtr()->matrix_size[2] = 1;
        cm1->getObjectPtr()->slice = sl;
        cm1->getObjectPtr()->position[0] = hdr.position[0] + sl*thick*hdr.slice_dir[0];
        cm1->getObjectPtr()->position[1] = hdr.position[1] + sl*thick*hdr.slice_dir[1];
        cm1->getObjectPtr()->position[2] = hdr.position[2] + sl*thick*hdr.slice_dir[2];

        //Create a new image for this slice - 4D data block [X,Y,1,CHA]
        GadgetContainerMessage< hoNDArray< std::complex<float> > >* cm2 = 
            new GadgetContainerMessage<hoNDArray< std::complex<float> > >();

        try{cm2->getObjectPtr()->create(X,Y,1,CHA);}
        catch (std::runtime_error &err){
            GEXCEPTION(err,"Unable to allocate new image\n");
            cm1->release();
            cm2->release();
            return GADGET_FAIL;
        }

        //Copy the data
        //TODO is there a way to do this without making copies?
        auto img_out = hoNDArray<std::complex<float> >();
        for (uint16_t c=0; c<CHA; c++) {
            memcpy(&(cm2->getObjectPtr()->at(c*slice_size)),&img_in(0,0,sl,c), slice_size);
        }

        //Chain the messages for the image header and image data
        cm1->cont(cm2);

        //Copy the meta container if needed
        if (m3) {
            GadgetContainerMessage< ISMRMRD::MetaContainer >* cm3 = 
                new GadgetContainerMessage< ISMRMRD::MetaContainer >();
            memcpy(cm3->getObjectPtr(), m3->getObjectPtr(), sizeof(ISMRMRD::MetaContainer));
            cm2->cont(cm3);
        }
 
        //Pass the image down the chain
        if (this->next()->putq(cm1) < 0) {
            return GADGET_FAIL;
        }
    }
    
    m1->release(); //We have copied all the data in this case
    return GADGET_OK;

}

GADGET_FACTORY_DECLARE(ImageSplitSlicesGadget)
}
