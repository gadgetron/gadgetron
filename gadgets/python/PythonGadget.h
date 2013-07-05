#pragma once 

#include "Gadget.h"
#include "Gadgetron.h"
#include "hoNDArray.h"
#include "GadgetMRIHeaders.h"
#include "PythonCommunicator.h"
#include "gadgetronpython_export.h"

#include <ismrmrd.h>
#include <boost/python.hpp>
#include <boost/algorithm/string.hpp>
#include <stdio.h>
#include <stdlib.h>
#include <complex>

namespace Gadgetron{

  template <class T> class PythonGadget : 
  public Gadget2<T, hoNDArray< std::complex<float> > >
    {
    protected:

      int process_config(ACE_Message_Block* mb)
      {
	communicator_ = PythonCommunicatorSingleton::instance();

	boost::shared_ptr<std::string> pypath        = this->get_string_value("python_path");
	boost::shared_ptr<std::string> pymod         = this->get_string_value("python_module");
	boost::shared_ptr<std::string> pyreffunc     = this->get_string_value("gadget_reference_function");
	boost::shared_ptr<std::string> pydatafunc    = this->get_string_value("input_function");
	boost::shared_ptr<std::string> pyconfigfunc  = this->get_string_value("config_function");

	GADGET_DEBUG2("Python Module          : %s\n", pymod.get()->c_str());
	GADGET_DEBUG2("Python Ref Function    : %s\n", pyreffunc.get()->c_str());
	GADGET_DEBUG2("Python Data Function   : %s\n", pydatafunc.get()->c_str());
	GADGET_DEBUG2("Python Config Function : %s\n", pyconfigfunc.get()->c_str());

	if (communicator_->addPath(*pypath.get()) != GADGET_OK) {
	  GADGET_DEBUG2("Failed to add paths in Gadget %s\n", this->module()->name());
	  return GADGET_FAIL;
	}

	if (communicator_->registerGadget(this, *pymod.get(),
					  *pyreffunc.get(), *pyconfigfunc.get(),
					  *pydatafunc.get()) != GADGET_OK) {
	  GADGET_DEBUG2("Failed to register Gadget (%s) with PythonCommunicator\n", this->module()->name());
	  return GADGET_FAIL;
	}

	if (communicator_->processConfig(this, mb) != GADGET_OK) {
	  GADGET_DEBUG2("Failed to process config in Python module of Gadget (%s)\n", this->module()->name());
	  return GADGET_FAIL;
	}

	return GADGET_OK;
      }

      int process(GadgetContainerMessage<T>* m1,
		  GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2)
      {
    
	//We want to avoid a deadlock for the Python GIL if this python call results in an output that the GadgetReference will not be able to get rid of.
	//This is kind of a nasty busy wait, maybe we should add an event handler to the NotificationStrategy of the Q or something, but for now, this will do it.
	while (this->next()->msg_queue()->is_full()) {
	  //GADGET_DEBUG2("Gadget (%s) sleeping while downstream Gadget (%s) does some work\n", this->module()->name(), this->next()->module()->name());
	  ACE_Time_Value tv(0,10000); //Sleep for 10ms while the downstream Gadget does some work
	  ACE_OS::sleep(tv);
	}

	//GADGET_DEBUG2("Process called in Gadget (%s)\n", this->module()->name());
	if (communicator_->process(this,m1,m2) != GADGET_OK) {
	  GADGET_DEBUG2("Failed to process data for Gadget (%s)\n", this->module()->name());
	  return GADGET_FAIL;
	}

	//GADGET_DEBUG2("Process done in Gadget (%s)\n", this->module()->name());
	return GADGET_OK;
      }
  
    private:
      PythonCommunicator* communicator_;
    };
  
  class EXPORTGADGETSPYTHON AcquisitionPythonGadget :
  public PythonGadget<ISMRMRD::AcquisitionHeader>
  {
  public:
    GADGET_DECLARE(AcquisitionPythonGadget);
  
  };

  class EXPORTGADGETSPYTHON ImagePythonGadget :
  public PythonGadget<ISMRMRD::ImageHeader>
  {
  public:
    GADGET_DECLARE(ImagePythonGadget);    
  };
}
