#pragma once 

#include "Gadget.h"
#include "hoNDArray.h"
#include "GadgetMRIHeaders.h"
#include "PythonCommunicator.h"
#include "gadgetronpython_export.h"

#include <ismrmrd/ismrmrd.h>
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
        boost::shared_ptr<std::string> pyclass       = this->get_string_value("python_class");

	GDEBUG("Python Module          : %s\n", pymod.get()->c_str());
	GDEBUG("Python Class           : %s\n", pyclass.get()->c_str());

	if (communicator_->addPath(*pypath.get()) != GADGET_OK) {
	  GDEBUG("Failed to add paths in Gadget %s\n", this->module()->name());
	  return GADGET_FAIL;
	}

	if (communicator_->registerGadget(this, *pymod.get(), *pyclass.get()) != GADGET_OK) {
	  GDEBUG("Failed to register Gadget (%s) with PythonCommunicator\n", this->module()->name());
	  return GADGET_FAIL;
	}

	if (communicator_->processConfig(this, mb) != GADGET_OK) {
	  GDEBUG("Failed to process config in Python module of Gadget (%s)\n", this->module()->name());
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
	  //GDEBUG("Gadget (%s) sleeping while downstream Gadget (%s) does some work\n", this->module()->name(), this->next()->module()->name());
	  ACE_Time_Value tv(0,10000); //Sleep for 10ms while the downstream Gadget does some work
	  ACE_OS::sleep(tv);
	}

	//GDEBUG("Process called in Gadget (%s)\n", this->module()->name());
	if (communicator_->process(this,m1,m2) != GADGET_OK) {
	  GDEBUG("Failed to process data for Gadget (%s)\n", this->module()->name());
	  return GADGET_FAIL;
	}

	//GDEBUG("Process done in Gadget (%s)\n", this->module()->name());
	return GADGET_OK;
      }
  
    private:
      PythonCommunicator* communicator_;
    };
  
  class EXPORTGADGETSPYTHON AcquisitionPythonGadget :
  public PythonGadget<ISMRMRD::AcquisitionHeader> {};

  class EXPORTGADGETSPYTHON ImagePythonGadget :
  public PythonGadget<ISMRMRD::ImageHeader> {};
}
