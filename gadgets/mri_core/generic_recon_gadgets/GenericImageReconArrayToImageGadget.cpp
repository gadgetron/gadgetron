
#include "GenericImageReconArrayToImageGadget.h"
#include "hoNDArray_reductions.h"

namespace Gadgetron {

    GenericImageReconArrayToImageGadget::GenericImageReconArrayToImageGadget() : BaseClass()
    {
    }

    GenericImageReconArrayToImageGadget::~GenericImageReconArrayToImageGadget()
    {
    }

    int GenericImageReconArrayToImageGadget::process_config(ACE_Message_Block* mb)
    {
        GADGET_CHECK_RETURN(BaseClass::process_config(mb) == GADGET_OK, GADGET_FAIL);
        return GADGET_OK;
    }

    int GenericImageReconArrayToImageGadget::processImageBuffer(ImageBufferType& ori)
    {
        GDEBUG_CONDITION_STREAM(this->verbose.value(), "GenericImageReconArrayToImageGadget::processImageBuffer(...) starts ... ");

        std::vector<std::string> processStr;
        std::vector<std::string> dataRole;

        std::vector<size_t> dims = ori.dimensions();
        GDEBUG_CONDITION_STREAM(this->verbose.value(), "GenericImageReconArrayToImageGadget, [Cha Slice Con Phase Rep Set Ave] = [" << dims[0] << " " << dims[1] << " " << dims[2] << " "
            << dims[3] << " " << dims[4] << " " << dims[5] << " "
            << dims[6] << "]");

        // sending out image array
        size_t CHA = ori.get_size(0);
        size_t SLC = ori.get_size(1);
        size_t CON = ori.get_size(2);
        size_t PHS = ori.get_size(3);
        size_t REP = ori.get_size(4);
        size_t SET = ori.get_size(5);
        size_t AVE = ori.get_size(6);

        size_t ave(0), set(0), rep(0), phs(0), con(0), slc(0), cha(0);
        std::vector<size_t> dim3D(3);

        for (ave = 0; ave < AVE; ave++)
        {
            for (set = 0; set < SET; set++)
            {
                for (rep = 0; rep < REP; rep++)
                {
                    for (phs = 0; phs < PHS; phs++)
                    {
                        for (con = 0; con < CON; con++)
                        {
                            for (slc = 0; slc < SLC; slc++)
                            {
                                for (cha = 0; cha < CHA; cha++)
                                {
                                    Image2DType& img = ori(cha, slc, con, phs, rep, set, ave);
                                    if (img.get_number_of_elements()>0)
                                    {
                                        T v = Gadgetron::nrm2(img);
                                        if (v < FLT_EPSILON) continue; // do not send out empty image

                                        Gadgetron::GadgetContainerMessage<ISMRMRD::ImageHeader>* cm1 = new Gadgetron::GadgetContainerMessage<ISMRMRD::ImageHeader>();
                                        Gadgetron::GadgetContainerMessage<ImgArrayType>* cm2 = new Gadgetron::GadgetContainerMessage<ImgArrayType>();
                                        Gadgetron::GadgetContainerMessage<ISMRMRD::MetaContainer>* cm3 = new Gadgetron::GadgetContainerMessage<ISMRMRD::MetaContainer>();

                                        cm1->cont(cm2);
                                        cm2->cont(cm3);

                                        // set the ISMRMRD image header
                                        memcpy(cm1->getObjectPtr(), &img.header_, sizeof(ISMRMRD::ISMRMRD_ImageHeader));

                                        // set the image data
                                        size_t RO = img.get_size(0);
                                        size_t E1 = img.get_size(1);
                                        size_t E2 = img.get_size(2);

                                        dim3D[0] = RO;
                                        dim3D[1] = E1;
                                        dim3D[2] = E2;

                                        cm2->getObjectPtr()->create(dim3D);
                                        memcpy(cm2->getObjectPtr()->get_data_ptr(), img.get_data_ptr(), img.get_number_of_bytes());

                                        // set the attributes
                                        *cm3->getObjectPtr() = img.attrib_;

                                        if (this->next()->putq(cm1) < 0)
                                        {
                                            cm1->release();
                                            return false;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        GADGET_CHECK_RETURN(this->releaseImageBuffer(ori), GADGET_FAIL);

        GDEBUG_CONDITION_STREAM(this->verbose.value(), "GenericImageReconArrayToImageGadget::process(...) ends ... ");

        return GADGET_OK;
    }

    int GenericImageReconArrayToImageGadget::close(unsigned long flags)
    {
        GDEBUG_CONDITION_STREAM(true, "GenericImageReconArrayToImageGadget - close(flags) : " << flags);

        if (BaseClass::close(flags) != GADGET_OK) return GADGET_FAIL;

        return GADGET_OK;
    }

    GADGET_FACTORY_DECLARE(GenericImageReconArrayToImageGadget)
}
