#include "ImageArraySplitGadget.h"

namespace Gadgetron{

ImageArraySplitGadget::ImageArraySplitGadget()
{

}

int ImageArraySplitGadget::process(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1)
{
    if (this->next()->putq(m1) < 0)
    {
        m1->release();
        return GADGET_FAIL;
    }

    return GADGET_OK;
}

int ImageArraySplitGadget::process( GadgetContainerMessage<IsmrmrdImageArray>* m1)
{
    
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

    //Each image will be [X,Y,Z,CHA] big
    std::vector<size_t> img_dims(4);
    img_dims[0] = X;
    img_dims[1] = Y;
    img_dims[2] = Z;
    img_dims[3] = CHA;

    //Loop over N, S and LOC
    for (uint16_t loc=0; loc < LOC; loc++) {
        for (uint16_t s=0; s < S; s++) {                
            for (uint16_t n=0; n < N; n++) {
        
                //Create a new image header and copy the header for this n, s and loc
                GadgetContainerMessage<ISMRMRD::ImageHeader>* cm1 = 
                        new GadgetContainerMessage<ISMRMRD::ImageHeader>();
                memcpy(cm1->getObjectPtr(), &imagearr.headers_(n,s,loc), sizeof(ISMRMRD::ImageHeader));

                //Create a new image image
                // and the 4D data block [X,Y,Z,CHA] for this n, s and loc
                GadgetContainerMessage< hoNDArray< std::complex<float> > >* cm2 = 
                        new GadgetContainerMessage<hoNDArray< std::complex<float> > >();

                try{cm2->getObjectPtr()->create(img_dims);}
                catch (std::runtime_error &err){
                    GEXCEPTION(err,"Unable to allocate new image\n");
                    cm1->release();
                    cm2->release();
                    return GADGET_FAIL;
                }
                memcpy(cm2->getObjectPtr()->get_data_ptr(), &imagearr.data_(0,0,0,0,n,s,loc), X*Y*Z*CHA*sizeof(std::complex<float>));
                //Chain them
                cm1->cont(cm2);

                //Create a new meta container if needed and copy
                if (imagearr.meta_.size()>0) {
                    GadgetContainerMessage< ISMRMRD::MetaContainer >* cm3 = 
                            new GadgetContainerMessage< ISMRMRD::MetaContainer >();
                    size_t mindex = loc*N*S + s*N + n;
                    *cm3->getObjectPtr() = imagearr.meta_[mindex];
                    cm2->cont(cm3);
                }

                //Pass the image down the chain
                if (this->next()->putq(cm1) < 0) {
                    return GADGET_FAIL;
                }
            }
        }
    }
    
    m1->release();
    return GADGET_OK;  

}

GADGET_FACTORY_DECLARE(ImageArraySplitGadget)
}
