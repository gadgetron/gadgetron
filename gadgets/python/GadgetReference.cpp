#include "Gadget.h"
#include "GadgetReference.h"
#include "GadgetContainerMessage.h"
#include "hoNDArray.h"
#include <ismrmrd/ismrmrd.h>
#include <ismrmrd/meta.h>

/* #include <boost/preprocessor/stringize.hpp> */
#include <boost/python.hpp>

namespace Gadgetron {

    GadgetReference::GadgetReference()
        : gadget_(nullptr)
    {
    }

    GadgetReference::~GadgetReference()
    {
    }

    int GadgetReference::return_recondata(boost::python::object rec) {
        auto m1 = new GadgetContainerMessage<IsmrmrdReconData>(boost::python::extract<IsmrmrdReconData>(rec)());
        if (gadget_) {
            if (gadget_->next()->putq(m1 ) == -1) {
                m1->release();
                return GADGET_FAIL;
            }
            else
                return GADGET_OK;
        }
        else {
            GDEBUG("IsmrmrdReconData returned from python, but no next gadget in chain");
            m1->release();
            return GADGET_OK;
        }
    }

    int GadgetReference::return_ismrmrd_image_array(boost::python::object rec)
    {
        auto m1 = new GadgetContainerMessage<IsmrmrdImageArray>(boost::python::extract<IsmrmrdImageArray>(rec)());
        if (gadget_)
        {
            if (gadget_->next()->putq(m1 ) == -1)
            {
                m1->release();
                return GADGET_FAIL;
            }
            else
                return GADGET_OK;
        }
        else
        {
            GDEBUG("IsmrmrdImageArray returned from python, but no next gadget in chain");
            m1->release();
            return GADGET_OK;
        }
    }

    template<class TH, class TD>
    int GadgetReference::return_data(TH header, boost::python::object arr, const char* meta)
    {
        GadgetContainerMessage< TH >* m1 = new GadgetContainerMessage< TH >;
        memcpy(m1->getObjectPtr(), &header, sizeof(TH));

        // this works because the python converter for hoNDArray<std::complex<float>>
        // is registered in the python_toolbox
        GadgetContainerMessage< hoNDArray< TD > >* m2;
        m2 = new GadgetContainerMessage< hoNDArray< TD > >(
            boost::python::extract<hoNDArray < TD > >(arr)());
        m1->cont(m2);

        if (meta) {
            GadgetContainerMessage< ISMRMRD::MetaContainer >* m3 =
                new GadgetContainerMessage< ISMRMRD::MetaContainer >;

            ISMRMRD::deserialize(meta, *m3->getObjectPtr());
            m2->cont(m3);
        }

        if (gadget_) {
            //GDEBUG("Returning data (%s)\n", gadget_->module()->name());
            if (gadget_->next()->putq(m1) == -1) {
                m1->release();
                //if (gadget_->next()->putq(m1) == -1) {
                /*
                  GDEBUG("Putting message on Queue failed (%s)\n", gadget_->module()->name());
                  GDEBUG("Message Q: low mark %d, high mark %d, message bytes %d, message count %d\n",
                  gadget_->next()->msg_queue()->low_water_mark(), gadget_->next()->msg_queue()->high_water_mark(),
                  gadget_->next()->msg_queue()->message_bytes(),gadget_->next()->msg_queue()->message_count());
                */
                //GDEBUG("FAIL Returning data (%s)\n", gadget_->module()->name());
                return GADGET_FAIL;
            }
            else {
                //GDEBUG("SUCCESS Returning data (%s)\n", gadget_->module()->name());

                return GADGET_OK;
            }
            //return gadget_->next()->putq(m1);
        }
        else {
            GDEBUG("Data received from python, but no Gadget registered for output\n");
            m1->release();
            return GADGET_OK;
        }

        return GADGET_OK;
    }

    int GadgetReference::return_acquisition(ISMRMRD::AcquisitionHeader acq, boost::python::object arr)
    {
        return return_data<ISMRMRD::AcquisitionHeader, std::complex<float> >(acq, arr, 0);
    }

    int GadgetReference::return_image_cplx(ISMRMRD::ImageHeader img, boost::python::object arr)
    {
        return return_data<ISMRMRD::ImageHeader, std::complex<float> >(img, arr, 0);
    }

    int GadgetReference::return_image_cplx_attr(ISMRMRD::ImageHeader img, boost::python::object arr, const char* meta)
    {
        return return_data<ISMRMRD::ImageHeader, std::complex<float> >(img, arr, meta);
    }


    int GadgetReference::return_image_float(ISMRMRD::ImageHeader img, boost::python::object arr)
    {
        return return_data<ISMRMRD::ImageHeader, float>(img, arr, 0);
    }

    int GadgetReference::return_image_float_attr(ISMRMRD::ImageHeader img, boost::python::object arr, const char* meta)
    {
        return return_data<ISMRMRD::ImageHeader, float>(img, arr, meta);
    }

    int GadgetReference::return_image_ushort(ISMRMRD::ImageHeader img, boost::python::object arr)
    {
        return return_data<ISMRMRD::ImageHeader, unsigned short>(img, arr, 0);
    }

    int GadgetReference::return_image_ushort_attr(ISMRMRD::ImageHeader img, boost::python::object arr, const char* meta)
    {
        return return_data<ISMRMRD::ImageHeader, unsigned short>(img, arr, meta);
    }

}
