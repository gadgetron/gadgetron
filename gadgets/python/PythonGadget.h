#pragma once 

#include <Python.h>
#include <boost/python.hpp>
#include <numpy/arrayobject.h>
#include <boost/algorithm/string.hpp>

#include "gadgetron_export.h"
#include "Gadget.h"
#include "Gadgetron.h"
#include "hoNDArray.h"
#include "../core/GadgetMRIHeaders.h"
#include "GadgetReference.h"

#include <stdio.h>
#include <stdlib.h>

#include <complex>

template <class T> class PythonGadget : 
public Gadget2<T, hoNDArray< std::complex<float> > >
{
 public:
  //GADGET_DECLARE(PythonGadget);
  //virtual ~PythonGadget();

 protected:
  bool last_python_call_success_;
  PyThreadState* py_interpreter_state_;


  int process_config(ACE_Message_Block* mb)
  {
    GADGET_DEBUG2("Initializing Threads (%s)\n", this->module()->name());

    if (PyEval_ThreadsInitialized()) {
      GADGET_DEBUG2("Acquire Lock (%s)\n", this->module()->name());
      PyEval_AcquireLock();
    } else {
      PyEval_InitThreads();
      GADGET_DEBUG2("Initialize (%s)\n", this->module()->name());
      Py_Initialize();
    }

    _import_array();

    GADGET_DEBUG2("New Interpreter (%s)\n", this->module()->name());

    
    boost::shared_ptr<std::string> pypath        = this->get_string_value("python_path");
    boost::shared_ptr<std::string> pymod         = this->get_string_value("python_module");
    boost::shared_ptr<std::string> pyreffunc     = this->get_string_value("gadget_reference_function");
    boost::shared_ptr<std::string> pydatafunc    = this->get_string_value("input_function");
    boost::shared_ptr<std::string> pyconfigfunc  = this->get_string_value("config_function");
    
    GADGET_DEBUG2("Python Module          : %s\n", pymod.get()->c_str());
    GADGET_DEBUG2("Python Ref Function    : %s\n", pyreffunc.get()->c_str());
    GADGET_DEBUG2("Python Data Function   : %s\n", pydatafunc.get()->c_str());
    GADGET_DEBUG2("Python Config Function : %s\n", pyconfigfunc.get()->c_str());
    
    try {
      //Let's first get the path set
      const char* gadgetron_home = ACE_OS::getenv("GADGETRON_HOME");
      
      std::string path_name;
      std::string path_cmd;
      if (gadgetron_home != 0) {
	path_name = std::string(gadgetron_home) + std::string("/lib"); 
	path_cmd = std::string("import sys;\nif (sys.path.count(\"") + path_name +
	  std::string("\") == 0):\n\tsys.path.append(\"") + path_name + std::string("\")\n"); 
	//GADGET_DEBUG2("Executing path command:\n%s\n", path_cmd.c_str());
	boost::python::exec(path_cmd.c_str(),boost::python::import("__main__").attr("__dict__"));
      }
      
      if (pypath.get()->size() > 0) { 
	std::vector<std::string> paths;
	boost::split(paths, *(pypath.get()), boost::is_any_of(";"));
	for (unsigned int i = 0; i < paths.size(); i++) {
	  path_cmd = std::string("import sys;\nif (sys.path.count(\"") + paths[i] +
	    std::string("\") == 0):\n\tsys.path.append(\"") + paths[i] + std::string("\")\n"); 
	  //GADGET_DEBUG2("Executing path command:\n%s\n", path_cmd.c_str());
	  boost::python::exec(path_cmd.c_str(),boost::python::import("__main__").attr("__dict__"));
	}
      }
      
      //Now Paths should be set. Let's load the module. 
      python_module_ = boost::python::import(pymod.get()->c_str());
      
      /* We will try t force a reload of the module */
      boost::python::import("__main__").attr("__dict__")[pymod.get()->c_str()] = python_module_;
      std::string tmp = std::string("reload(") + std::string(pymod.get()->c_str()) + std::string(")\n");
      //GADGET_DEBUG2("Reloading with command: %s\n", tmp.c_str());
      boost::python::exec(tmp.c_str(),boost::python::import("__main__").attr("__dict__"));
      
      python_set_gadget_reference_function_ = python_module_.attr(pyreffunc.get()->c_str());
      gadget_reference_.set_gadget(this);
      python_set_gadget_reference_function_(gadget_reference_);
      
      python_input_function_ =  python_module_.attr(pydatafunc.get()->c_str());
      
      if (pyconfigfunc.get()->size() > 0) {
	python_config_function_ =  python_module_.attr(pyconfigfunc.get()->c_str());
	
	python_config_function_(boost::python::handle<>(PyString_FromString(mb->rd_ptr())));
      }

      py_interpreter_state_ = PyEval_SaveThread();//
      //PyEval_ReleaseLock(); //Release Python Lock

    } catch(boost::python::error_already_set const &) {
      py_interpreter_state_ = PyEval_SaveThread();//PyEval_ReleaseLock();
      GADGET_DEBUG1("Error loading python modules\n");
      PyErr_Print();
      return GADGET_FAIL;
    }
    
    last_python_call_success_ = true;
    return GADGET_OK;
  }


  int process(GadgetContainerMessage<T>* m1,
	      GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2)
  {
   
    if (!last_python_call_success_) {
      GADGET_DEBUG2("Gadget (%s) last python call failed. Bailing out!\n", this->module()->name());
      //m1->release();
      return GADGET_FAIL;
    }
 
    try {
      PyEval_RestoreThread(py_interpreter_state_);

      std::vector<unsigned int> dims = (*(m2->getObjectPtr()->get_dimensions().get()));
      std::vector<int> dims2(dims.size());
      for (unsigned int i = 0; i < dims.size(); i++) dims2[dims.size()-i-1] = static_cast<int>(dims[i]);
      
      boost::python::object obj(boost::python::handle<>(PyArray_FromDims(dims2.size(), &dims2[0], PyArray_CFLOAT)));
      boost::python::object data = boost::python::extract<boost::python::numeric::array>(obj);
      
      //Copy data
      memcpy(PyArray_DATA(data.ptr()), m2->getObjectPtr()->get_data_ptr(), m2->getObjectPtr()->get_number_of_elements()*sizeof(float)*2);
      
      //Get Header
      T acq = *m1->getObjectPtr();      
      m1->release(); 
     
      python_input_function_(acq, data);

      py_interpreter_state_ = PyEval_SaveThread();//PyEval_ReleaseLock();

    } catch(boost::python::error_already_set const &) { 
      py_interpreter_state_ = PyEval_SaveThread();//PyEval_ReleaseLock();
      GADGET_DEBUG2("Passing data on to python module failed  (Gadget: %s)\n", this->module()->name());
      PyErr_Print();
      last_python_call_success_ = false;
      return GADGET_OK;//We will bail out next time we enter
    }
    
    //It is enough to put the first one, since they are linked
    /*
    if (this->next()->putq(m1) == -1) {
      m1->release();
      ACE_ERROR_RETURN( (LM_ERROR,
			 ACE_TEXT("%p\n"),
			 ACE_TEXT("PythonGadget::process, passing data on to next gadget")),
			-1);
    }
    */

    
    return GADGET_OK;
  }
  
 private:
  boost::python::object python_module_;
  boost::python::object python_set_gadget_reference_function_;
  boost::python::object python_input_function_;
  boost::python::object python_config_function_;
  
  GadgetReference gadget_reference_;

};



class EXPORTGADGETSCORE AcquisitionPythonGadget : 
public PythonGadget<GadgetMessageAcquisition>
{
 public:
  GADGET_DECLARE(AcquisitionPythonGadget);
  
};

class EXPORTGADGETSCORE ImagePythonGadget : 
public PythonGadget<GadgetMessageImage>
{
 public:
  GADGET_DECLARE(ImagePythonGadget);

};
