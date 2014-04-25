/*
*       ComplexToFloatAttribGadget.cpp
*       Author: Hui Xue
*/

#include "GadgetIsmrmrdReadWrite.h"
#include "ComplexToFloatAttribGadget.h"
#include "hoNDArray_math_util.h"

namespace Gadgetron
{
    ComplexToFloatAttribGadget::ComplexToFloatAttribGadget()
    {
    }

    ComplexToFloatAttribGadget::~ComplexToFloatAttribGadget()
    {
    }

    int ComplexToFloatAttribGadget::process(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1, GadgetContainerMessage< hoNDArray< ValueType > >* m2, GadgetContainerMessage<GtImageAttribType>* m3)
    {
        GadgetContainerMessage<hoNDArray< float > > *cm2 = new GadgetContainerMessage<hoNDArray< float > >();

        boost::shared_ptr< std::vector<size_t> > dims = m2->getObjectPtr()->get_dimensions();

        try
        {
            cm2->getObjectPtr()->create(dims);
        }
        catch (std::runtime_error &err)
        {
            GADGET_DEBUG_EXCEPTION(err,"Unable to create float storage in ComplexToFloatAttribGadget");
            return GADGET_FAIL;
        }

        switch (m1->getObjectPtr()->image_type)
        {
            case ISMRMRD::TYPE_MAGNITUDE:
            {
                GADGET_CHECK_RETURN(Gadgetron::absolute(*m2->getObjectPtr(), *cm2->getObjectPtr()), GADGET_FAIL);
            }
            break;

            case ISMRMRD::TYPE_REAL:
            {
                GADGET_CHECK_RETURN(Gadgetron::complex_to_real(*m2->getObjectPtr(), *cm2->getObjectPtr()), GADGET_FAIL);
            }
            break;

            case ISMRMRD::TYPE_IMAG:
            {
                GADGET_CHECK_RETURN(Gadgetron::complex_to_imag(*m2->getObjectPtr(), *cm2->getObjectPtr()), GADGET_FAIL);
            }
            break;

            case ISMRMRD::TYPE_PHASE:
            {
                GADGET_CHECK_RETURN(Gadgetron::argument(*m2->getObjectPtr(), *cm2->getObjectPtr()), GADGET_FAIL);
            }
            break;

            default:
                GADGET_DEBUG2("Unknown image type %d, bailing out\n",m1->getObjectPtr()->image_type);
                m1->release();
                cm2->release();
                return GADGET_FAIL;
        }

        m1->cont(cm2);
        cm2->cont(m3);

        m2->cont(NULL);
        m2->release();

        m1->getObjectPtr()->image_data_type = ISMRMRD::DATA_FLOAT;

        if (this->next()->putq(m1) == -1)
        {
            m1->release();
            GADGET_DEBUG1("Unable to put unsigned short magnitude image on next gadgets queue");
            return GADGET_FAIL;
        }

        return GADGET_OK;
    }

    GADGET_FACTORY_DECLARE(ComplexToFloatAttribGadget)
}
