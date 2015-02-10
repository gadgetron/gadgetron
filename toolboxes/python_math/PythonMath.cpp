#include "PythonMath.h"
#include "log.h"

#include <stdexcept>

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <boost/python.hpp>
//#include <boost/numpy.hpp>
#include <numpy/numpyconfig.h>
#include <numpy/arrayobject.h>

namespace bp = boost::python;

namespace Gadgetron
{

  template <typename T> void bp_object_to_hondarray(bp::object& obj, hoNDArray< T >* arr)
  {
    if (!arr) {
      throw std::runtime_error("bp_object_to_hondarray received null pointer array");
    }

    if (sizeof(T) != PyArray_ITEMSIZE((PyArrayObject*)obj.ptr())) {
      GERROR("sizeof(T): %d, ITEMSIZE: %d\n", sizeof(T), PyArray_ITEMSIZE((PyArrayObject*)obj.ptr()));
      throw std::runtime_error("bp_object_to_hondarray: python object and array data type sizes do not match");
    }

    size_t ndim = PyArray_NDIM((PyArrayObject*)obj.ptr());
    std::vector<size_t> dims(ndim);
    for (size_t i = 0; i < ndim; i++) {
      dims[ndim-i-1] = PyArray_DIM((PyArrayObject*)obj.ptr(),i);
    }

    if (!arr->dimensions_equal(&dims)) {
      arr->create(dims);
    }

    memcpy(arr->get_data_ptr(), PyArray_DATA((PyArrayObject*)obj.ptr()), sizeof(T)*arr->get_number_of_elements());
  }

  template <typename T> int get_numpy_type();
  template <> int get_numpy_type< std::complex<float> >() {return NPY_COMPLEX64;}
  template <> int get_numpy_type< std::complex<double> >() {return NPY_COMPLEX128;}
  template <> int get_numpy_type< float >() {return NPY_FLOAT32;}
  template <> int get_numpy_type< unsigned int >() {return NPY_UINT32;}

  template <typename T> bp::object hondarray_to_bp_object(hoNDArray< T >* arr)
  {
    if (!arr) {
      throw std::runtime_error("hondarray_to_bp_object: null pointer passed for arr");
    }

    size_t ndim = arr->get_number_of_dimensions();
    std::vector<int> dims2(ndim);
    for (unsigned int i = 0; i < ndim; i++) dims2[ndim-i-1] = static_cast<int>(arr->get_size(i));
    bp::object obj(bp::handle<>(PyArray_FromDims(dims2.size(), &dims2[0], get_numpy_type<T>())));
    
    if (sizeof(T) != PyArray_ITEMSIZE((PyArrayObject*)obj.ptr())) {
      GERROR("sizeof(T): %d, ITEMSIZE: %d\n", sizeof(T), PyArray_ITEMSIZE((PyArrayObject*)obj.ptr()));
      throw std::runtime_error("hondarray_to_bp_object: python object and array data type sizes do not match");
    }

    //Copy data
    memcpy(PyArray_DATA((PyArrayObject*)obj.ptr()), arr->get_data_ptr(), arr->get_number_of_elements()*sizeof(std::complex<float>));

    return obj;
  }
  
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
    //boost::numpy::initialize();
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
    GDEBUG("Called calculate grappa unmixing.\n");
    
    if (!source_data) {
      throw std::runtime_error("Source data is NULL");
    }

    if (source_data->get_number_of_dimensions() != 3) {
      throw std::runtime_error("source_data input array has wrong number of dimensions");
    }

    PyGILState_STATE gstate;
    gstate = PyGILState_Ensure();

    try {
      bp::object mod = bp::import("ismrmrdtools.grappa"); 
      bp::object fcn = mod.attr("calculate_grappa_unmixing");
      bp::object source_data_obj = hondarray_to_bp_object(source_data);
      bp::tuple kernel_size_obj;
      bp::object data_mask_obj; //None
      bp::object csm_obj; //None
      bp::object target_data_obj; //None
      
      if (kernel_size) {
	kernel_size_obj = bp::make_tuple(kernel_size[0],kernel_size[1]);
      } else {
	kernel_size_obj = bp::make_tuple(4,5);
      }

      if (data_mask) data_mask_obj = hondarray_to_bp_object(data_mask);
      if (csm) csm_obj = hondarray_to_bp_object(csm);
      if (target_data) target_data_obj = hondarray_to_bp_object(target_data);

      bp::object result = fcn(source_data_obj,acc_factor,kernel_size_obj, data_mask_obj, csm_obj, regularization_factor, target_data_obj);

      if (bp::len(result) != 2) {
	throw std::runtime_error("Wrong number of results from python function");
      }

      bp::tuple result_ex = bp::extract<bp::tuple>(result);      
      bp::object umx = bp::extract<bp::object>(result_ex[0]);
      bp::object gmap = bp::extract<bp::object>(result_ex[1]);

      if (PyArray_NDIM((PyArrayObject*)umx.ptr()) != 3) {
	throw std::runtime_error("Returned unmixing coefficients have wrong number of dimensions");
      } 

      if (PyArray_NDIM((PyArrayObject*)gmap.ptr()) != 2) {
	throw std::runtime_error("Returned gmap coefficients have wrong number of dimensions");
      } 

      bp_object_to_hondarray(umx, unmix_out);

      if (gmap_out) {
	bp_object_to_hondarray(gmap, gmap_out);
      }

      GDEBUG("Grappa unmixing calculation returned.\n");
    } catch(bp::error_already_set const &) {
      GERROR("Error calling python code in PythonMath\n");
      PyErr_Print();
      PyGILState_Release(gstate);
    } catch (std::exception const & e) {
      GERROR("Run time error in PythonMath\n");
      PyGILState_Release(gstate);
      throw;
    }
    PyGILState_Release(gstate);

  }
  
}
