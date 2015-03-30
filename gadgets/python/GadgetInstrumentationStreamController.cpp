#include "GadgetInstrumentationStreamController.h"
#include "EndGadget.h"
#include "AcquisitionFinishGadget.h"
#include "ImageFinishGadget.h"
#include "log.h"
#include "hoNDArray.h"
#include "GadgetMessageInterface.h"
#include "GadgetMRIHeaders.h"
#include <string.h>
#include <ismrmrd/ismrmrd.h>
#include <ismrmrd/meta.h>

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
    if (this->prepend_gadget("ImageFinishFloat","gadgetron_mricore","ImageFinishGadget") != GADGET_OK) return GADGET_FAIL;
    this->find_gadget("ImageFinishFloat")->pass_on_undesired_data(true);

    /*
    if (this->prepend_gadget("ImageFinishCplx","gadgetron_mricore","ImageFinishGadgetCPLX") != GADGET_OK) return GADGET_FAIL;
    this->find_gadget("ImageFinishCplx")->pass_on_undesired_data(true);

    if (this->prepend_gadget("ImageFinishUShort","gadgetron_mricore","ImageFinishGadgetUSHORT") != GADGET_OK) return GADGET_FAIL;
    this->find_gadget("ImageFinishUShort")->pass_on_undesired_data(true);
    */

    if (this->prepend_gadget("AcquisitionFinish","gadgetron_mricore","AcquisitionFinishGadget") != GADGET_OK) return GADGET_FAIL;
    this->find_gadget("AcquisitionFinish")->pass_on_undesired_data(true);


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
    
    Gadget* g = dynamic_cast<Gadget*>(m->writer());//Get the gadget out of the module

    //We will set this very high to prevent race conditions in "mixed environments" such as when using Python or Matlab in Gadgets
    g->msg_queue()->high_water_mark(ACE_Message_Queue_Base::DEFAULT_HWM*100000);
    
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

  template <class T1, class T2, class T3> int GadgetInstrumentationStreamController::return_data(ACE_Message_Block* mb)
  {
    static int counter = 0;
    GadgetContainerMessage<T1>* m1 = AsContainerMessage<T1>(mb);
    GadgetContainerMessage<T2>* m2 = AsContainerMessage<T2>(mb->cont());
    GadgetContainerMessage<T3>* m3 = AsContainerMessage<T3>(m2->cont());

    if (!m1 || !m2) {
      GERROR("Unable to convert input container messages");
      return GADGET_FAIL;
    }
    
    {
      GILLock lock;
      try {
	if (m3) {
	  std::stringstream str;
	  ISMRMRD::serialize(*m3->getObjectPtr(), str);
	  python_gadget_.attr("put_next")(*m1->getObjectPtr(),m2->getObjectPtr(),str.str());
	} else {
	  python_gadget_.attr("put_next")(*m1->getObjectPtr(),m2->getObjectPtr());
	}
      } catch(boost::python::error_already_set const &) {
	GERROR("Passing data on to python wrapper gadget failed\n");
	PyErr_Print();
	return GADGET_FAIL;
      }
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

    GadgetContainerMessage<ISMRMRD::ImageHeader>* m_tmp = 0;

    switch (m0->getObjectPtr()->id) {
    case (GADGET_MESSAGE_ACQUISITION):
      if (0 != this->return_data<ISMRMRD::AcquisitionHeader, hoNDArray< std::complex<float> >, ISMRMRD::MetaContainer >(m0->cont()) )
	{
	  GERROR("Unable to convert and return GADGET_MESSAGE_ACQUISITON\n");
	  m0->release();
	  return GADGET_FAIL;
	}
      break;

    case (GADGET_MESSAGE_ISMRMRD_IMAGE):
      m_tmp = AsContainerMessage<ISMRMRD::ImageHeader>(m0->cont());
      if (!m_tmp) {
	GERROR("Error converting header of GADGET_MESSAGE_ISMRMRD_IMAG\n");
	mb->release();
	return GADGET_FAIL;
      }
      switch (m_tmp->getObjectPtr()->data_type) {
	
      case (ISMRMRD::ISMRMRD_USHORT):
	if (0 != this->return_data<ISMRMRD::ImageHeader, hoNDArray< uint16_t >, ISMRMRD::MetaContainer >(m0->cont()) )
	  {
	    GERROR("Unable to convert and return GADGET_MESSAGE_ISMRMRD_IMAGE\n");
	    m0->release();
	    return GADGET_FAIL;
	  }
	break;
      case (ISMRMRD::ISMRMRD_SHORT):
	if (0 != this->return_data<ISMRMRD::ImageHeader, hoNDArray< int16_t >, ISMRMRD::MetaContainer >(m0->cont()) )
	  {
	    GERROR("Unable to convert and return GADGET_MESSAGE_ISMRMRD_IMAGE\n");
	    m0->release();
	    return GADGET_FAIL;
	  }
	break;
      case (ISMRMRD::ISMRMRD_UINT):
	if (0 != this->return_data<ISMRMRD::ImageHeader, hoNDArray< uint32_t >, ISMRMRD::MetaContainer >(m0->cont()) )
	  {
	    GERROR("Unable to convert and return GADGET_MESSAGE_ISMRMRD_IMAGE\n");
	    m0->release();
	    return GADGET_FAIL;
	  }
	break;
      case (ISMRMRD::ISMRMRD_INT):
	if (0 != this->return_data<ISMRMRD::ImageHeader, hoNDArray< int32_t >, ISMRMRD::MetaContainer >(m0->cont()) )
	  {
	    GERROR("Unable to convert and return GADGET_MESSAGE_ISMRMRD_IMAGE\n");
	    m0->release();
	    return GADGET_FAIL;
	  }
	break;
      case (ISMRMRD::ISMRMRD_FLOAT):
	if (0 != this->return_data<ISMRMRD::ImageHeader, hoNDArray< float >, ISMRMRD::MetaContainer >(m0->cont()) )
	  {
	    GERROR("Unable to convert and return GADGET_MESSAGE_ISMRMRD_IMAGE\n");
	    m0->release();
	    return GADGET_FAIL;
	  }
	break;
      case (ISMRMRD::ISMRMRD_DOUBLE):
	if (0 != this->return_data<ISMRMRD::ImageHeader, hoNDArray< double >, ISMRMRD::MetaContainer >(m0->cont()) )
	  {
	    GERROR("Unable to convert and return GADGET_MESSAGE_ISMRMRD_IMAGE\n");
	    m0->release();
	    return GADGET_FAIL;
	  }
	break;
      case (ISMRMRD::ISMRMRD_CXFLOAT):
	if (0 != this->return_data<ISMRMRD::ImageHeader, hoNDArray< std::complex<float> >, ISMRMRD::MetaContainer >(m0->cont()) )
	  {
	    GERROR("Unable to convert and return GADGET_MESSAGE_ISMRMRD_IMAGE\n");
	    m0->release();
	    return GADGET_FAIL;
	  }
	break;
	
      case (ISMRMRD::ISMRMRD_CXDOUBLE):
	if (0 != this->return_data<ISMRMRD::ImageHeader, hoNDArray< std::complex<double> >, ISMRMRD::MetaContainer >(m0->cont()) )
	  {
	    GERROR("Unable to convert and return GADGET_MESSAGE_ISMRMRD_IMAGE\n");
	    m0->release();
	    return GADGET_FAIL;
	  }
	break;
      }
      break;
    case (GADGET_MESSAGE_CLOSE):
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


  int GadgetInstrumentationStreamController::put_config(const char* config)
  {
    size_t l = std::strlen(config);
    ACE_Message_Block* mb = new ACE_Message_Block(l+1);
    memcpy(mb->wr_ptr(),config,l+1);
    mb->wr_ptr(l+1);
    mb->set_flags(Gadget::GADGET_MESSAGE_CONFIG);
    if (stream_.put(mb) == -1) {
      GERROR("Failed to put configuration on stream, too long wait, %d\n",  ACE_OS::last_error () ==  EWOULDBLOCK);
      mb->release();
      return GADGET_FAIL;
    }
    return GADGET_OK;
  }

  template<class TH, class TD>
  int GadgetInstrumentationStreamController::put_data(TH header, boost::python::object arr, const char* meta)
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


    ACE_Time_Value wait = ACE_OS::gettimeofday() + ACE_Time_Value(0,10000); //10ms from now
    if (stream_.put(m1) == -1) {
      GERROR("Failed to put stuff on stream, too long wait, %d\n",  ACE_OS::last_error () ==  EWOULDBLOCK);
      m1->release();
      return GADGET_FAIL;
    }
    return GADGET_OK;
  }

  int GadgetInstrumentationStreamController::put_acquisition(ISMRMRD::AcquisitionHeader acq, 
							     boost::python::object arr, const char* meta)
  {
    return put_data<ISMRMRD::AcquisitionHeader, std::complex<float> >(acq, arr);
  }

  int GadgetInstrumentationStreamController::put_image_cplx(ISMRMRD::ImageHeader img, 
							    boost::python::object arr, const char* meta)
  {
    return put_data<ISMRMRD::ImageHeader, std::complex<float> >(img, arr, meta);
  }

  int GadgetInstrumentationStreamController::put_image_float(ISMRMRD::ImageHeader img, 
							    boost::python::object arr, const char* meta)
  {
    return put_data<ISMRMRD::ImageHeader, float >(img, arr, meta);
  }

  int GadgetInstrumentationStreamController::put_image_ushort(ISMRMRD::ImageHeader img, 
							    boost::python::object arr, const char* meta)
  {
    return put_data<ISMRMRD::ImageHeader, unsigned short >(img, arr, meta);
  }

}
