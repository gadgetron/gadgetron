#pragma once

#include <ace/OS_NS_unistd.h>

#include "Gadget.h"
#include "hoNDArray.h"
#include "GadgetReference.h"
#include "gadgetronpython_export.h"
#include "python_toolbox.h"

#include <ismrmrd/ismrmrd.h>
#include <ismrmrd/meta.h>
#include <boost/python.hpp>

namespace Gadgetron {

    /// This PythonGadget is the gateway for c++ to call python
    class EXPORTGADGETSPYTHON PythonGadget : public BasicPropertyGadget
    {
    public:
        GADGET_DECLARE(PythonGadget);

        /*
          We are overloading this function from the base class to be able to capture a copy
          of the properties that should be passed on to the Python class itself.
         */
        virtual int set_parameter(const char* name, const char* val, bool trigger = true) {
            GadgetPropertyBase* p = this->find_property(name);
            if (p) {
                //This is a property, pass it on to the Gadget base class
                return Gadget::set_parameter(name, val, trigger);
            }
            else {
                //This is probably information for the Python class itself
                this->parameters_python_[std::string(name)] = std::string(val);
            }
            return GADGET_OK;
        }

    protected:
        int process_config(ACE_Message_Block* mb)
        {
            // start python interpreter
            if (initialize_python() != GADGET_OK) {
                GDEBUG("Failed to initialize Python in Gadget %s\n", this->module()->name());
                return GADGET_FAIL;
            }

            // ensure boost can convert between hoNDArrays and NumPy arrays automatically
            register_converter<hoNDArray< std::complex<float> > >();
            register_converter<hoNDArray< float > >();
            register_converter<hoNDArray< unsigned short > >();
            register_converter<hoNDArray< uint32_t > >();
            register_converter<hoNDArray< ISMRMRD::AcquisitionHeader > >();
            register_converter<hoNDArray< ISMRMRD::ImageHeader > >();

            // ensure boost can convert ISMRMRD headers automatically
            register_converter<ISMRMRD::ImageHeader>();
            register_converter<ISMRMRD::AcquisitionHeader>();
            register_converter<ISMRMRD::ISMRMRD_WaveformHeader>();
            register_converter<ISMRMRD::Waveform>();

            register_converter<IsmrmrdReconData>();
            register_converter<IsmrmrdImageArray>();

            register_converter< std::vector< std::complex<float> > >();
            register_converter< std::vector< float > >();
            register_converter< std::vector< unsigned short > >();
            register_converter< std::vector<ISMRMRD::MetaContainer> >();
            register_converter< std::vector<ISMRMRD::Waveform> >();

            std::string pypath = python_path.value();
            std::string pymod = python_module.value();
            std::string pyclass = python_class.value();

            GDEBUG("Python Path            : %s\n", pypath.c_str());
            GDEBUG("Python Module          : %s\n", pymod.c_str());
            GDEBUG("Python Class           : %s\n", pyclass.c_str());

            if (add_python_path(pypath) != GADGET_OK) {
                GDEBUG("Failed to add paths in Gadget %s\n", this->module()->name());
                return GADGET_FAIL;
            }

            std::string module_name = pymod;
            class_name = pyclass;

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
#if defined PYVER && PYVER == 3
                // prefix reload call for Python 3
                tmp = std::string("from importlib import reload;") + tmp;
#endif
                //GDEBUG("Reloading with command: %s\n", tmp.c_str());
                boost::python::exec(tmp.c_str(), boost::python::import("__main__").attr("__dict__"));

                gadget_ref_ = boost::shared_ptr<GadgetReference>(new GadgetReference());
                gadget_ref_->set_gadget(this);

                // Create instance of class (passing gadget reference)
                class_ = module_.attr(class_name.c_str())(gadget_ref_.get());
                // Increment reference count of Python class so that both the C++
                // destructor and the interpreter can decrement its reference count
                boost::python::incref(class_.ptr());

            }
            catch (boost::python::error_already_set const &) {
                GDEBUG("Error loading python modules in Gadget %s\n", this->module()->name());
                PyErr_Print();
                return GADGET_FAIL;
            }

            //Transfer all properties/parameters to Python gadget
            std::map<std::string, std::string>::iterator it;
            it = parameters_python_.begin();
            while (it != parameters_python_.end()) {
                std::string var_name = it->first;
                std::string var_val = it->second;
                try {
                    boost::python::object set_parameter_fn = class_.attr("set_parameter");
                    boost::python::object ignored = set_parameter_fn(var_name, var_val);
                }
                catch (boost::python::error_already_set const &) {
                    GERROR("Error setting PythonGadget parameters in Gadget %s\n", this->module()->name());
                    PyErr_Print();
                    return GADGET_FAIL;
                }
                it++;
            }

            try {
                // retrieve and call python gadget's process_config method
                boost::python::object process_config_fn = class_.attr("process_config");
                boost::python::object ignored = process_config_fn(
                    boost::python::object(std::string(mb->rd_ptr())));
            }
            catch (boost::python::error_already_set const &) {
                GERROR("Error calling process_config in Gadget %s\n", this->module()->name());
                PyErr_Print();
                return GADGET_FAIL;
            }

            return GADGET_OK;
        }

        int process(GadgetContainerMessage<IsmrmrdReconData>* recon_data)
        {
            if (!recon_data) {
                GERROR("Received null pointer to data block");
                return GADGET_FAIL;
            }

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
                ACE_Time_Value tv(0, 10000);
                ACE_OS::sleep(tv);
            }

            GILLock lock;
            try {
                boost::python::object process_fn = class_.attr("process");
                auto pyrecon_data = boost::python::object(*recon_data->getObjectPtr());
                int res = boost::python::extract<int>(process_fn(pyrecon_data));
                if (res != GADGET_OK) {
                    GDEBUG("Gadget (%s) Returned from python call with error\n",
                        this->module()->name());
                    return GADGET_FAIL;
                }
                //Else we are done with this now.
                recon_data->release();
            }
            catch (boost::python::error_already_set const &) {
                GDEBUG("Passing data on to python module failed\n");
                PyErr_Print();
                return GADGET_FAIL;
            }
            return GADGET_OK;
        }

        int process(GadgetContainerMessage<IsmrmrdImageArray>* recon_data)
        {
            if (!recon_data)
            {
                GERROR("Received null pointer to data block");
                return GADGET_FAIL;
            }

            while (this->next()->msg_queue()->is_full())
            {
                ACE_Time_Value tv(0, 10000);
                ACE_OS::sleep(tv);
            }

            GILLock lock;
            try
            {
                boost::python::object process_fn = class_.attr("process");
                auto pyrecon_data = boost::python::object(*recon_data->getObjectPtr());
                int res = boost::python::extract<int>(process_fn(pyrecon_data));
                if (res != GADGET_OK)
                {
                    GDEBUG("Gadget (%s) Returned from python call with error\n",
                        this->module()->name());
                    return GADGET_FAIL;
                }
                recon_data->release();
            }
            catch (boost::python::error_already_set const &)
            {
                GDEBUG("Passing IsmrmrdImageArray on to python module failed\n");
                PyErr_Print();
                return GADGET_FAIL;
            }

            return GADGET_OK;
        }

        template <typename H, typename D> int process(GadgetContainerMessage<H>* hmb,
            GadgetContainerMessage< hoNDArray< D > >* dmb,
            GadgetContainerMessage< ISMRMRD::MetaContainer>* mmb = nullptr)
        {
            if (!dmb) {
                GERROR("Received null pointer to data block");
                return GADGET_FAIL;
            }

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
                ACE_Time_Value tv(0, 10000);
                ACE_OS::sleep(tv);
            }

            H head = *hmb->getObjectPtr();
            hoNDArray< D > *data = dmb->getObjectPtr();
            ISMRMRD::MetaContainer* meta = 0;
            if (mmb) {
                meta = mmb->getObjectPtr();
            }

            GILLock lock;
            try {
                boost::python::object process_fn = class_.attr("process");
                int res;
                if (meta) {
                    std::stringstream str;
                    ISMRMRD::serialize(*meta, str);
                    res = boost::python::extract<int>(process_fn(head, data, str.str()));
                }
                else {
                    res = boost::python::extract<int>(process_fn(head, data));
                }
                if (res != GADGET_OK) {
                    GDEBUG("Gadget (%s) Returned from python call with error\n",
                        this->module()->name());
                    return GADGET_FAIL;
                }
                //Else we are done with this now.
                hmb->release();
            }
            catch (boost::python::error_already_set const &) {
                GDEBUG("Passing data on to python module failed\n");
                PyErr_Print();
                return GADGET_FAIL;
            }
            return GADGET_OK;
        }



        virtual int process(ACE_Message_Block* mb);

    protected:
        GADGET_PROPERTY(python_module, std::string, "Python module containing the Python Gadget class to be loaded", "");
        GADGET_PROPERTY(python_class, std::string, "Python class to load from python module", "");
        GADGET_PROPERTY(python_path, std::string, "Path(s) to add to the to the Python search path", "");

    private:
        boost::python::object module_;
        boost::python::object class_;
        boost::shared_ptr<GadgetReference> gadget_ref_;
        std::string class_name;
        int process_image(GadgetContainerMessage<ISMRMRD::ImageHeader>* hmi);
        /*
          We are going to keep a copy of the parameters in this gadget that are not properties.
          They should be passed on to the Python class.
         */
        std::map< std::string, std::string> parameters_python_;
    };
}
