#include "PythonMath.h"
#include "log.h"

#include <boost/python.hpp>

#include <numpy/numpyconfig.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

namespace Gadgetron
{
  PythonMath* PythonMath::instance_ = 0;

  PythonMath* PythonMath::instance()
  {
    if (!instance_)
      {
	instance_ = new PythonMath();
      }
    return instance_;
  }

  PythonMath::PythonMath()
  {
    Py_Initialize();
    _import_array();
    
    PyEval_InitThreads();
    
    //Swap out and return current thread state and release the GIL
    //Must be done, otherwise subsequent calls to PyGILState_Ensure() will not be guaranteed to acuire lock
    PyThreadState* tstate = PyEval_SaveThread();
    if (!tstate) {
      GDEBUG("Error occurred returning lock to Python\n");
    }

  }

  void PythonMath::calculate_grappa_unmixing(hoNDArray< std::complex<float> >* source_data, unsigned int acc_factor, hoNDArray< std::complex<float> >* unmix_out,
					    unsigned int* kernel_size, hoNDArray< unsigned int >* data_mask, hoNDArray< std::complex<float> >* csm,
					    float regularization_factor, hoNDArray< std::complex<float> >* target_data, 
					    hoNDArray< float >* gmap_out)
  {
    GDEBUG("Called calculate grappa unmixing\n");
    PyGILState_STATE gstate;
    gstate = PyGILState_Ensure();
    
    boost::python::object mod = boost::python::import("ismrmrdtools.grappa"); 
    boost::python::object fcn = mod.attr("calculate_grappa_unmixing");
    
    size_t ndim = source_data->get_number_of_dimensions();
    std::vector<int> dims2(ndim);
    for (unsigned int i = 0; i < ndim; i++) dims2[ndim-i-1] = static_cast<int>(source_data->get_size(i));
    
    boost::python::object obj(boost::python::handle<>(PyArray_FromDims(dims2.size(), &dims2[0], NPY_COMPLEX64)));
    //boost::python::object data = boost::python::extract<boost::python::numeric::array>(obj);
    
    //Copy data
    memcpy(PyArray_DATA((PyArrayObject*)obj.ptr()), source_data->get_data_ptr(), source_data->get_number_of_elements()*sizeof(std::complex<float>));
  
    boost::python::object ret = fcn(obj,acc_factor);
    
    PyGILState_Release(gstate);

  }
  
}
