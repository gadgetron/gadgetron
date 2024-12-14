#include "python_toolbox.h"

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/numpyconfig.h>
#include <numpy/arrayobject.h>

#include <boost/thread/mutex.hpp>
#include <boost/algorithm/string.hpp>

// #include "Gadget.h"             // for GADGET_OK/FAIL


namespace Gadgetron
{

static bool python_initialized = false;
static bool numpy_initialized = false;
static boost::mutex python_initialize_mtx;
static boost::mutex numpy_initialize_mtx;

void initialize_python(void)
{
    // lock here so only one thread can initialize/finalize Python
    boost::mutex::scoped_lock lock(python_initialize_mtx);

    if (!python_initialized) {
        Py_Initialize();
        initialize_numpy();

        //Swap out and return current thread state and release the GIL
        //Must be done, otherwise subsequent calls to PyGILState_Ensure()
        //will not be guaranteed to acquire lock
//        PyThreadState* tstate = PyEval_SaveThread();
        PyThreadState* tstate = PyThreadState_Get();
        if (!tstate) {
            throw std::runtime_error("Error occurred returning lock to Python\n");
        }

        PyEval_ReleaseThread(tstate);
        python_initialized = true; // interpreter successfully initialized
    }
}

void initialize_numpy(void)
{
    // lock here so only one thread can initialize NumPy
    boost::mutex::scoped_lock lock(numpy_initialize_mtx);

    if (!numpy_initialized) {
        _import_array();    // import NumPy
        numpy_initialized = true; // numpy successfully initialized
    }
}

void finalize_python(void)
{
    // lock here so only one thread can initialize/finalize Python
    boost::mutex::scoped_lock lock(python_initialize_mtx);

    if (python_initialized) {
        Py_Finalize();
        python_initialized = false;
    }
}

void add_python_path(const std::string& path)
{
    GILLock lock;   // Lock the GIL

    std::string add_path_cmd;
    if (path.size() > 0) {
        std::vector<std::string> paths;
        boost::split(paths, path, boost::is_any_of(";"));
        for (unsigned int i = 0; i < paths.size(); i++) {
            add_path_cmd = std::string("import sys;\nif sys.path.count(\"") +
                    paths[i] + std::string("\") == 0:\n\tsys.path.insert(0, \"") +
                    paths[i] + std::string("\")\n");
            boost::python::exec(add_path_cmd.c_str(),
                    boost::python::import("__main__").attr("__dict__"));
        }
    }

}


/// Adapted from http://stackoverflow.com/a/6576177/1689220
std::string pyerr_to_string(void)
{
    PyObject *exc, *val, *tb;
    bp::object formatted_list, formatted;
    PyErr_Fetch(&exc, &val, &tb);
    PyErr_NormalizeException(&exc, &val, &tb);

    // wrap exception, value, traceback with bp::handle for auto memory management
    bp::handle<> hexc(exc), hval(bp::allow_null(val)), htb(bp::allow_null(tb));

    bp::object traceback(bp::import("traceback"));
    bp::object format_exception(traceback.attr("format_exception"));
    formatted_list = format_exception(hexc, hval, htb);
    formatted = bp::str("").join(formatted_list);
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

npy_intp* NumPyArray_STRIDES(PyObject* obj)
{
    return PyArray_STRIDES((PyArrayObject*)obj);
}

PyObject* NumPyArray_FromAny(PyObject* op, PyArray_Descr* dtype, int min_depth, int max_depth, int requirements, PyObject* context){
  return PyArray_FromAny(op, dtype, min_depth, max_depth, requirements, context);
}

/// Wraps PyArray_ITEMSIZE
int NumPyArray_ITEMSIZE(PyObject* obj)
{
    return PyArray_ITEMSIZE((PyArrayObject*)obj);
}

npy_intp NumPyArray_SIZE(PyObject* obj)
{
    return PyArray_SIZE((PyArrayObject*)obj);
}

/// Wraps PyArray_SimpleNew
PyObject* NumPyArray_SimpleNew(int nd, npy_intp* dims, int typenum)
{
    return PyArray_SimpleNew(nd, dims, typenum);
}

/// Wraps PyArray_SimpleNew
PyObject* NumPyArray_EMPTY(int nd, npy_intp* dims, int typenum, int fortran)
{
    return PyArray_EMPTY(nd, dims, typenum,fortran);
}

} // namespace Gadgetron
