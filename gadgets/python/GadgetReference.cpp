#include "Gadget.h"
#include "GadgetReference.h"
#include "GadgetContainerMessage.h"
#include "hoNDArray.h"
#include <ismrmrd/ismrmrd.h>
#include <ismrmrd/meta.h>

/* #include <boost/preprocessor/stringize.hpp> */
#include <boost/python.hpp>

namespace Gadgetron{

  GadgetReference::GadgetReference()
    : gadget_(nullptr)
  {
  }

  GadgetReference::~GadgetReference()
  {
  }

  template<class T>
  int GadgetReference::return_data(T header, boost::python::object arr, const char* meta)
  {
    GadgetContainerMessage< T >* m1 = new GadgetContainerMessage< T >;
    memcpy(m1->getObjectPtr(), &header, sizeof(T));

    // this works because the python converter for hoNDArray<std::complex<float>>
    // is registered in the python_toolbox
    GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2;
    m2 = new GadgetContainerMessage< hoNDArray< std::complex<float> > >(
            boost::python::extract<hoNDArray <std::complex<float> > >(arr)());
    m1->cont(m2);

    if (meta) {
      GadgetContainerMessage< ISMRMRD::MetaContainer >* m3 = 
	new GadgetContainerMessage< ISMRMRD::MetaContainer >;
      
      ISMRMRD::deserialize(meta, *m3->getObjectPtr());
      m2->cont(m3);
    }

    if (gadget_) {
      //ACE_Time_Value wait = ACE_OS::gettimeofday() + ACE_Time_Value(0,1000); //1ms from now
      ACE_Time_Value nowait (ACE_OS::gettimeofday ());
      //GDEBUG("Returning data (%s)\n", gadget_->module()->name());
      if (gadget_->next()->putq(m1,&nowait) == -1) {
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
      } else {
	//GDEBUG("SUCCESS Returning data (%s)\n", gadget_->module()->name());

	return GADGET_OK;
      }
      //return gadget_->next()->putq(m1);
    } else {
      GDEBUG("Data received from python, but no Gadget registered for output\n");
      m1->release();
      return GADGET_OK;
    }

    return GADGET_OK;
  }

  int GadgetReference::return_acquisition(ISMRMRD::AcquisitionHeader acq, boost::python::object arr)
  {
    return return_data<ISMRMRD::AcquisitionHeader>(acq, arr, 0);
  }

  int GadgetReference::return_image(ISMRMRD::ImageHeader img, boost::python::object arr)
  {
    return return_data<ISMRMRD::ImageHeader>(img, arr, 0);
  }

  int GadgetReference::return_image_attr(ISMRMRD::ImageHeader img, boost::python::object arr, const char* meta)
  {
    return return_data<ISMRMRD::ImageHeader>(img, arr, meta);
  }

  template int GadgetReference::return_data<ISMRMRD::AcquisitionHeader>(ISMRMRD::AcquisitionHeader, boost::python::object, const char*);
  template int GadgetReference::return_data<ISMRMRD::ImageHeader>(ISMRMRD::ImageHeader, boost::python::object, const char*);
}
