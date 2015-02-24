#include "python_toolbox.h"

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/numpyconfig.h>
#include <numpy/arrayobject.h>


namespace Gadgetron
{

static bool initialized = false;

void initialize_python(void)
{
    if (!initialized) {
        Py_Initialize();
        _import_array();    // import NumPy

        PyEval_InitThreads();

        //Swap out and return current thread state and release the GIL
        //Must be done, otherwise subsequent calls to PyGILState_Ensure()
        //will not be guaranteed to acquire lock
        PyThreadState* tstate = PyEval_SaveThread();
        if (!tstate) {
            GDEBUG("Error occurred returning lock to Python\n");
        }
        initialized = true;
    }
}

/// Adapted from http://stackoverflow.com/a/6576177/1689220
std::string pyerr_to_string(void)
{
    PyObject *exc, *val, *tb;
    bp::object formatted_list, formatted;
    PyErr_Fetch(&exc, &val, &tb);
    // wrap exception, value, traceback with bp::handle for auto memory management
    bp::handle<> hexc(exc), hval(bp::allow_null(val)), htb(bp::allow_null(tb));
    // import "traceback" module
    bp::object traceback(bp::import("traceback"));
    if (!tb) {
        bp::object format_exception_only(traceback.attr("format_exception_only"));
        formatted_list = format_exception_only(hexc, hval);
    } else {
        bp::object format_exception(traceback.attr("format_exception"));
        formatted_list = format_exception(hexc, hval, htb);
    }
    formatted = bp::str("\n").join(formatted_list);
    return bp::extract<std::string>(formatted);
}

/// Wraps PyArray_NDIM
int NumPyArray_NDIM(PyObject* obj)
{
    return PyArray_NDIM((PyArrayObject*)obj);
}

/// Wraps PyArray_DIM
npy_intp NumPyArray_DIM(PyObject* obj, int i)
{
    return PyArray_DIM((PyArrayObject*)obj, i);
}

/// Wraps PyArray_DATA
void* NumPyArray_DATA(PyObject* obj)
{
    return PyArray_DATA((PyArrayObject*)obj);
}

/// Wraps PyArray_ITEMSIZE
int NumPyArray_ITEMSIZE(PyObject* obj)
{
    return PyArray_ITEMSIZE((PyArrayObject*)obj);
}

/// Wraps PyArray_SimpleNew
PyObject* NumPyArray_SimpleNew(int nd, npy_intp* dims, int typenum)
{
    return PyArray_SimpleNew(nd, dims, typenum);
}

}
