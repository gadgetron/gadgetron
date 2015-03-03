#pragma once

#include "Gadget.h"
#include "hoNDArray.h"
#include "GadgetReference.h"
#include "gadgetronpython_export.h"
#include "python_toolbox.h"

#include <ismrmrd/ismrmrd.h>
#include <boost/python.hpp>

namespace Gadgetron{

  template <class T> class PythonGadget : 
  public Gadget2<T, hoNDArray< std::complex<float> > >
    {
    protected:

      int process_config(ACE_Message_Block* mb)
      {
          if (initialize_python() != GADGET_OK) {
            GDEBUG("Failed to initialize Python in Gadget %s\n", this->module()->name());
            return GADGET_FAIL;
          }

          // ensure boost can convert between hoNDArrays and NumPy arrays automatically
          register_converter<hoNDArray<std::complex<float> > >();
          // ensure boost can convert ISMRMRD headers automatically
          register_converter<ISMRMRD::ImageHeader>();
          register_converter<ISMRMRD::AcquisitionHeader>();

          boost::shared_ptr<std::string> pypath        = this->get_string_value("python_path");
          boost::shared_ptr<std::string> pymod         = this->get_string_value("python_module");
          boost::shared_ptr<std::string> pyclass       = this->get_string_value("python_class");

          GDEBUG("Python Path            : %s\n", pypath.get()->c_str());
          GDEBUG("Python Module          : %s\n", pymod.get()->c_str());
          GDEBUG("Python Class           : %s\n", pyclass.get()->c_str());

        if (add_python_path(*pypath.get()) != GADGET_OK) {
            GDEBUG("Failed to add paths in Gadget %s\n", this->module()->name());
            return GADGET_FAIL;
        }

        std::string module_name = *pymod.get();
        std::string class_name = *pyclass.get();

        if (module_name.size() == 0) {
            GDEBUG("Null module name received in Gadget %s\n", this->module()->name());
            return GADGET_FAIL;
        }
        if (class_name.size() == 0) {
            GDEBUG("Null class name received in Gadget %s\n", this->module()->name());
            return GADGET_FAIL;
        }

        GILLock lock;
        try {
            module_ = boost::python::import(module_name.c_str());

            // Reload the module so changes take place at Gadgetron runtime
            boost::python::import("__main__").attr("__dict__")[module_name.c_str()] = module_;
            std::string tmp = std::string("reload(") + std::string(module_name.c_str()) + std::string(")\n");

            //GDEBUG("Reloading with command: %s\n", tmp.c_str());
            boost::python::exec(tmp.c_str(), boost::python::import("__main__").attr("__dict__"));

            gadget_ref_ = boost::shared_ptr<GadgetReference>(new GadgetReference());
            gadget_ref_->set_gadget(this);

            // Create instance of class (passing gadget reference)
            class_ = module_.attr(class_name.c_str())(gadget_ref_.get());
            // Increment reference count of Python class so that both the C++
            // destructor and the interpreter can decrement its reference count
            boost::python::incref(class_.ptr());

        } catch (boost::python::error_already_set const &) {
            GDEBUG("Error loading python modules in Gadget %s\n", this->module()->name());
            PyErr_Print();
            return GADGET_FAIL;
        }

        try {
            // retrieve and call python gadget's process_config method
            boost::python::object process_config_fn = class_.attr("process_config");
            boost::python::object ignored = process_config_fn(
                    boost::python::object(std::string(mb->rd_ptr())));
        } catch (boost::python::error_already_set const &) {
            GDEBUG("Error calling process_config in Gadget %s\n", this->module()->name());
            PyErr_Print();
            return GADGET_FAIL;
        }

        return GADGET_OK;
      }

      int process(GadgetContainerMessage<T>* m1,
            GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2)
      {
          // We want to avoid a deadlock for the Python GIL if this python call
          // results in an output that the GadgetReference will not be able to
          // get rid of.
          // This is kind of a nasty busy wait, maybe we should add an event
          // handler to the NotificationStrategy of the Q or something, but
          // for now, this will do it.
          while (this->next()->msg_queue()->is_full()) {
              // GDEBUG("Gadget (%s) sleeping while downstream Gadget (%s) does some work\n",
              //        this->module()->name(), this->next()->module()->name());
              // Sleep for 10ms while the downstream Gadget does some work
              ACE_Time_Value tv(0,10000);
              ACE_OS::sleep(tv);
          }

          T head = *m1->getObjectPtr();
          hoNDArray<std::complex<float> > *data = m2->getObjectPtr();

          GILLock lock;
          try {
              boost::python::object process_fn = class_.attr("process");
              if (boost::python::extract<int>(process_fn(head, data)) != GADGET_OK) {
                  GDEBUG("Gadget (%s) Returned from python call with error\n",
                          this->module()->name());
                  return GADGET_FAIL;
              }
              //Else we are done with this now.
              m1->release();
          } catch(boost::python::error_already_set const &) {
              GDEBUG("Passing data on to python module failed\n");
              PyErr_Print();
              return GADGET_FAIL;
          }

          //GDEBUG("Process done in Gadget (%s)\n", this->module()->name());
          return GADGET_OK;
      }

    private:
      boost::python::object module_;
      boost::python::object class_;
      boost::shared_ptr<GadgetReference> gadget_ref_;
    };

  class EXPORTGADGETSPYTHON AcquisitionPythonGadget :
  public PythonGadget<ISMRMRD::AcquisitionHeader> {};

  class EXPORTGADGETSPYTHON ImagePythonGadget :
  public PythonGadget<ISMRMRD::ImageHeader> {};
}
