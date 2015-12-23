#include "python_toolbox.h"

#include "Gadget.h"             // for GADGET_OK/FAIL
#include "gadgetron_paths.h"    // for get_gadgetron_home()
#include "gadgetron_config.h"   // for GADGETRON_PYTHON_PATH

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/numpyconfig.h>
#include <numpy/arrayobject.h>

#include <boost/thread/mutex.hpp>
#include <boost/algorithm/string.hpp>


namespace Gadgetron
{

static bool python_initialized = false;
static bool numpy_initialized = false;
static boost::mutex python_initialize_mtx;
static boost::mutex numpy_initialize_mtx;

int initialize_python(void)
{
    // lock here so only one thread can initialize Python
    boost::mutex::scoped_lock lock(python_initialize_mtx);

    if (!python_initialized) {
        Py_Initialize();
        initialize_numpy();

        PyEval_InitThreads();

        //Swap out and return current thread state and release the GIL
        //Must be done, otherwise subsequent calls to PyGILState_Ensure()
        //will not be guaranteed to acquire lock
        PyThreadState* tstate = PyEval_SaveThread();
        if (!tstate) {
            GDEBUG("Error occurred returning lock to Python\n");
            return GADGET_FAIL;
        }

        python_initialized = true; // interpreter successfully initialized

        //Let's first get the path set for the library folder
        std::string gadgetron_home = get_gadgetron_home();
        std::string path_name = gadgetron_home + std::string("/") + std::string(GADGETRON_PYTHON_PATH);

        if (gadgetron_home.size() != 0) {
            if (add_python_path(path_name) == GADGET_FAIL) {
                GDEBUG("python_toolbox failed to add path %s\n", path_name.c_str());
                return GADGET_FAIL;
            }
        }
    }
    return GADGET_OK;
}

int initialize_numpy(void)
{
    // lock here so only one thread can initialize NumPy
    boost::mutex::scoped_lock lock(numpy_initialize_mtx);

    if (!numpy_initialized) {
        _import_array();    // import NumPy
        numpy_initialized = true; // numpy successfully initialized
    }
    return GADGET_OK;
}

int add_python_path(const std::string& path)
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
            //GDEBUG("Executing path command:\n%s\n", path_cmd.c_str());
            boost::python::exec(add_path_cmd.c_str(),
                    boost::python::import("__main__").attr("__dict__"));
        }
    }

    return GADGET_OK;
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

}
