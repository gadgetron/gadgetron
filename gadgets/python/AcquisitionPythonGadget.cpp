#include <Python.h>
#include <numpy/arrayobject.h>

#include "AcquisitionPythonGadget.h"
#include "Gadgetron.h"

AcquisitionPythonGadget::~AcquisitionPythonGadget()
{
  Py_Finalize();
}

int AcquisitionPythonGadget::process_config(ACE_Message_Block* mb)
{
  Py_Initialize();
  _import_array();

  boost::shared_ptr<std::string> pymod       = get_string_value("python_module");
  boost::shared_ptr<std::string> pyreffunc   = get_string_value("gadget_reference_function");
  boost::shared_ptr<std::string> pydatafunc  = get_string_value("input_function");

  GADGET_DEBUG2("Python Module        : %s\n", pymod.get()->c_str());
  GADGET_DEBUG2("Python Ref Function  : %s\n", pyreffunc.get()->c_str());
  GADGET_DEBUG2("Python Data Function : %s\n", pydatafunc.get()->c_str());

  try {
    
    python_module_ = boost::python::import(pymod.get()->c_str());
    
    python_set_gadget_reference_function_ = python_module_.attr(pyreffunc.get()->c_str());
    gadget_reference_.set_gadget(this);
    python_set_gadget_reference_function_(gadget_reference_);
    
    python_input_function_ =  python_module_.attr(pydatafunc.get()->c_str());

  } catch(boost::python::error_already_set const &) { 
    GADGET_DEBUG1("Error loading python modules\n");
    PyErr_Print();
    return GADGET_FAIL;
  }
  
  return GADGET_OK;
}


int AcquisitionPythonGadget
::process(GadgetContainerMessage<GadgetMessageAcquisition>* m1,
	  GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2)
{
  
  try {
    std::vector<unsigned int> dims = (*(m2->getObjectPtr()->get_dimensions().get()));
    std::vector<int> dims2(dims.size());
    for (unsigned int i = 0; i < dims.size(); i++) dims2[i] = static_cast<int>(dims[i]);
  
    boost::python::object obj(boost::python::handle<>(PyArray_FromDims(dims2.size(), &dims2[0], PyArray_CFLOAT)));
    boost::python::object data = boost::python::extract<boost::python::numeric::array>(obj);

    //Copy data
    memcpy(PyArray_DATA(data.ptr()), m2->getObjectPtr()->get_data_ptr(), m2->getObjectPtr()->get_number_of_elements()*sizeof(float)*2);
  
    //Get Acquisition Header
    GadgetMessageAcquisition acq = *m1->getObjectPtr();

    python_input_function_(acq, data);

  } catch(boost::python::error_already_set const &) { 
    m1->release();
    GADGET_DEBUG1("Passing data on to python module failed\n");
    PyErr_Print();
    return GADGET_FAIL;
  }

  //It is enough to put the first one, since they are linked
  if (this->next()->putq(m1) == -1) {
    m1->release();
    ACE_ERROR_RETURN( (LM_ERROR,
		       ACE_TEXT("%p\n"),
		       ACE_TEXT("AcquisitionPythonGadget::process, passing data on to next gadget")),
		      -1);
  }

  return 0;
}


GADGET_FACTORY_DECLARE(AcquisitionPythonGadget)
