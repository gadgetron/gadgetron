/*
*       ComplexToFloatGadget.cpp
*       Author: Hui Xue
*/

#include "GadgetIsmrmrdReadWrite.h"
#include "ComplexToFloatGadget.h"
#include "hoNDArray_elemwise.h"

namespace Gadgetron
{
    ComplexToFloatGadget::ComplexToFloatGadget()
    {
      this->msg_queue()->high_water_mark((size_t)(12.0 * 1024 * 1024 * 1024));
    }

    ComplexToFloatGadget::~ComplexToFloatGadget()
    {
    }

    int ComplexToFloatGadget::process(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1, GadgetContainerMessage< hoNDArray< ValueType > >* m2)
    {
        GadgetContainerMessage<hoNDArray< float > > *cm2 = new GadgetContainerMessage<hoNDArray< float > >();

        boost::shared_ptr< std::vector<size_t> > dims = m2->getObjectPtr()->get_dimensions();

        try
        {
            cm2->getObjectPtr()->create(dims);
        }
        catch (std::runtime_error &err)
        {
            GEXCEPTION(err,"Unable to create float storage in ComplexToFloatGadget");
            return GADGET_FAIL;
        }

        switch (m1->getObjectPtr()->image_type)
        {
            case ISMRMRD::ISMRMRD_IMTYPE_MAGNITUDE:
            {
                GADGET_CHECK_EXCEPTION_RETURN(Gadgetron::abs(*m2->getObjectPtr(), *cm2->getObjectPtr()), GADGET_FAIL);
            }
            break;

            case ISMRMRD::ISMRMRD_IMTYPE_REAL:
            {
                GADGET_CHECK_EXCEPTION_RETURN(Gadgetron::complex_to_real(*m2->getObjectPtr(), *cm2->getObjectPtr()), GADGET_FAIL);
            }
            break;

            case ISMRMRD::ISMRMRD_IMTYPE_IMAG:
            {
                GADGET_CHECK_EXCEPTION_RETURN(Gadgetron::complex_to_imag(*m2->getObjectPtr(), *cm2->getObjectPtr()), GADGET_FAIL);
            }
            break;

            case ISMRMRD::ISMRMRD_IMTYPE_PHASE:
            {
                GADGET_CHECK_EXCEPTION_RETURN(Gadgetron::argument(*m2->getObjectPtr(), *cm2->getObjectPtr()), GADGET_FAIL);
            }
            break;

            default:
                GDEBUG("Unknown image type %d, bailing out\n",m1->getObjectPtr()->image_type);
                m1->release();
                cm2->release();
                return GADGET_FAIL;
        }

        GadgetContainerMessage<ISMRMRD::MetaContainer>* m3 = AsContainerMessage<ISMRMRD::MetaContainer>(m2->cont());

        m1->cont(cm2);
        if(m3) cm2->cont(m3);

        m2->cont(NULL);
        m2->release();

        m1->getObjectPtr()->data_type = ISMRMRD::ISMRMRD_FLOAT;

        if (this->next()->putq(m1) == -1)
        {
            m1->release();
            GDEBUG("Unable to put unsigned short magnitude image on next gadgets queue");
            return GADGET_FAIL;
        }

        return GADGET_OK;
    }

    GADGET_FACTORY_DECLARE(ComplexToFloatGadget)
}
