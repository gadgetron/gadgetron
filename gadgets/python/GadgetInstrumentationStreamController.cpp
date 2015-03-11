#include "GadgetInstrumentationStreamController.h"
#include "EndGadget.h"
#include "AcquisitionFinishGadget.h"
#include "ImageFinishGadget.h"
#include "log.h"
#include "hoNDArray.h"
#include "GadgetMessageInterface.h"
#include "GadgetMRIHeaders.h"

namespace Gadgetron
{
  GadgetInstrumentationStreamController::GadgetInstrumentationStreamController()
  {
    if (this->open() != GADGET_OK) {
      throw std::runtime_error("Unable to initialize GadgetInstrumentationStreamController");
    }
  }

  GadgetInstrumentationStreamController::~GadgetInstrumentationStreamController()
  {
    if (this->close() != GADGET_OK) {
      throw std::runtime_error("Unable to shut down sream in  GadgetInstrumentationStreamController");
    }
  }
  
  int GadgetInstrumentationStreamController::open()
  {
    GadgetModule *head = 0;
    GadgetModule *tail = 0;
    
    if (tail == 0) {
      Gadget* eg = new EndGadget();
      if (eg) {
	eg->set_controller(this);
      }
      
      ACE_NEW_RETURN(tail,
		     ACE_Module<ACE_MT_SYNCH>( ACE_TEXT("EndGadget"),
					       eg ),
		     -1);
      
      stream_.open(0,head,tail);
    }

    //Adding some gadgets to "capture data and return to the stream"
    this->prepend_gadget("ImageFinishFloat","gadgetron_mricore","ImageFinishGadgetFLOAT");
    this->set_parameter("ImageFinishFloat","pass_on_undesired_data","true");
    this->prepend_gadget("AcquisitionFinish","gadgetron_mricore","AcquisitionFinishGadget");
    this->set_parameter("AcquisitionFinish","pass_on_undesired_data","true");

    return GADGET_OK;
  }

  int GadgetInstrumentationStreamController::close()
  {
    stream_.close(1); //Shutdown gadgets and wait for them
    return GADGET_OK;
  }

  int GadgetInstrumentationStreamController::prepend_gadget(const char* gadgetname,
							   const char* dllname, 
							   const char* classname)
  {
    GadgetModule* m = create_gadget_module(dllname,
					   classname,
					   gadgetname);
      
    if (!m) {
      GERROR("Failed to create GadgetModule from %s:%s\n",
	     classname,
	     dllname);
      return GADGET_FAIL;
    }
    
    if (stream_.push(m) < 0) {
      GERROR("Failed to push Gadget %s onto stream\n", gadgetname);
      delete m;
      return GADGET_FAIL;
    }
    return GADGET_OK;
  }

  template <class T1, class T2> int GadgetInstrumentationStreamController::return_data(ACE_Message_Block* mb)
  {
    GadgetContainerMessage<T1>* m1 = AsContainerMessage<T1>(mb);
    GadgetContainerMessage<T2>* m2 = AsContainerMessage<T2>(mb->cont());

    if (!m1 || !m2) {
      GERROR("Unable to convert input container messages");
      return GADGET_FAIL;
    }
    
    GILLock lock;
    try {
      python_gadget_.attr("put_next")(*m1->getObjectPtr(),m2->getObjectPtr());
    } catch(boost::python::error_already_set const &) {
      GERROR("Passing data on to python wrapper gadget failed\n");
      PyErr_Print();
      return GADGET_FAIL;
    }

    return GADGET_OK;
  }

  int GadgetInstrumentationStreamController::output_ready(ACE_Message_Block* mb)
  {
    GadgetContainerMessage<GadgetMessageIdentifier>* m0 = AsContainerMessage<GadgetMessageIdentifier>(mb);
    if (!m0) {
      GERROR("Unable to extract GadgetMessageIdentifier\n");
      mb->release();
      return GADGET_FAIL;
    }

    switch (m0->getObjectPtr()->id)
      {
      case (GADGET_MESSAGE_ACQUISITION):
	if (0 != this->return_data<ISMRMRD::AcquisitionHeader, hoNDArray< std::complex<float> > >(m0->cont()) )
	  {
	    GERROR("Unable to convert and return GADGET_MESSAGE_ACQUISITON\n");
	    m0->release();
	    return GADGET_FAIL;
	  }
	break;
      case (GADGET_MESSAGE_IMAGE_REAL_USHORT):
	if (0 != this->return_data<ISMRMRD::ImageHeader, hoNDArray< ACE_UINT16 > >(m0->cont()) )
	  {
	    GERROR("Unable to convert and return GADGET_MESSAGE_IMAGE_REAL_SHORT\n");
	    m0->release();
	    return GADGET_FAIL;
	  }
	break;
      case (GADGET_MESSAGE_IMAGE_REAL_FLOAT):
	if (0 != this->return_data<ISMRMRD::ImageHeader, hoNDArray< float > >(m0->cont()) )
	  {
	    GERROR("Unable to convert and return GADGET_MESSAGE_IMAGE_REAL_FLOAT");
	    m0->release();
	    return GADGET_FAIL;
	  }
	break;
      case (GADGET_MESSAGE_IMAGE_CPLX_FLOAT):
	if (0 != this->return_data<ISMRMRD::ImageHeader, hoNDArray< std::complex<float> > >(m0->cont()) )
	  {
	    GERROR("Unable to convert and return GADGET_MESSAGE_IMAGE_CPLX_FLOAT\n");
	    m0->release();
	    return GADGET_FAIL;
	  }
	break;
      default:
	GERROR("Unsupported message ID (%d) encountered\n", m0->getObjectPtr()->id);
	mb->release();
	return GADGET_FAIL;
      }

    mb->release();
    return GADGET_OK;
  }

  void GadgetInstrumentationStreamController::set_parameter(const char* gadgetname, const char* parameter, const char* value)
  {
    Gadget* g = this->find_gadget(gadgetname);
    if (!g) {
      throw std::runtime_error("Unable to find Gadget for setting parameter");
    }
    g->set_parameter(parameter,value,false);
  }


  template<class T>
  int GadgetInstrumentationStreamController::put_data(T header, boost::python::object arr)
  {
    GadgetContainerMessage< T >* m1 = new GadgetContainerMessage< T >;
    memcpy(m1->getObjectPtr(), &header, sizeof(T));

    // this works because the python converter for hoNDArray<std::complex<float>>
    // is registered in the python_toolbox
    GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2;
    m2 = new GadgetContainerMessage< hoNDArray< std::complex<float> > >(
            boost::python::extract<hoNDArray <std::complex<float> > >(arr)());
    m1->cont(m2);

    ACE_Time_Value wait = ACE_OS::gettimeofday() + ACE_Time_Value(0,10000); //10ms from now
    if (stream_.put(m1) == -1) {
      GERROR("Failed to put stuff on stream, too long wait, %d\n",  ACE_OS::last_error () ==  EWOULDBLOCK);
      m1->release();
      return GADGET_FAIL;
    }
    return GADGET_OK;
  }

  int GadgetInstrumentationStreamController::put_acquisition(ISMRMRD::AcquisitionHeader acq, 
								boost::python::object arr)
  {
    return put_data<ISMRMRD::AcquisitionHeader>(acq, arr);
  }

  int GadgetInstrumentationStreamController::put_image(ISMRMRD::ImageHeader img, 
							 boost::python::object arr)
  {
    return put_data<ISMRMRD::ImageHeader>(img, arr);
  }

  template int GadgetInstrumentationStreamController::put_data<ISMRMRD::AcquisitionHeader>(ISMRMRD::AcquisitionHeader, boost::python::object);

  template int GadgetInstrumentationStreamController::put_data<ISMRMRD::ImageHeader>(ISMRMRD::ImageHeader, boost::python::object);

}
