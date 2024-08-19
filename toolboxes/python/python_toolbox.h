#pragma once

#include "python_export.h"
#include "log.h"

#include <boost/python.hpp>
namespace bp = boost::python;

namespace Gadgetron
{

/// Initialize Python and NumPy. Called by each PythonFunction constructor
EXPORTPYTHON void initialize_python(void);
/// Initialize NumPy
EXPORTPYTHON void initialize_numpy(void);
/// Finalize Python, Called by user expclictly
EXPORTPYTHON void finalize_python(void);
/// Add a path to the PYTHONPATH
EXPORTPYTHON void add_python_path(const std::string& path);

/// Extracts the exception/traceback to build and return a std::string
EXPORTPYTHON std::string pyerr_to_string(void);

}

// Include converters after declaring above functions
namespace Gadgetron {
/// Utility class for RAII handling of the Python GIL. Usage:
///
///    GILLock lg;  // at the top of a block
///
    class GILLock {
    public:
        GILLock() { gstate_ = PyGILState_Ensure(); }

        ~GILLock() { PyGILState_Release(gstate_); }

    private:
        // noncopyable
        GILLock(const GILLock &);

        GILLock &operator=(const GILLock &);

        PyGILState_STATE gstate_;

    };

}

#include "python_converters.h"


namespace Gadgetron {
/// Base class for templated PythonFunction class. Do not use directly.
class PythonFunctionBase
{
protected:
    PythonFunctionBase(const std::string& module, const std::string& funcname)
    {
        initialize_python(); // ensure Python and NumPy are initialized
        GILLock lg; // Lock the GIL, releasing at the end of constructor
        try {
            // import the module and load the function
            bp::object mod(bp::import(module.c_str()));
            fn_ = mod.attr(funcname.c_str());
        } catch (const bp::error_already_set&) {
            std::string err = pyerr_to_string();
            GERROR(err.c_str());
            throw std::runtime_error(err);
        }
    }

    bp::object fn_;
};

/// PythonFunction for multiple return types (std::tuple)
template <typename... ReturnTypes>
class PythonFunction : public PythonFunctionBase
{
public:
    typedef std::tuple<ReturnTypes...> TupleType;

    PythonFunction(const std::string& module, const std::string& funcname)
      : PythonFunctionBase(module, funcname)
    {
        // register the tuple return type converter
        register_converter<TupleType>();
    }

    template <typename... TS>
    TupleType operator()(const TS&... args)
    {
        // register type converter for each parameter type
        register_converter<TS...>();
        GILLock lg; // lock GIL and release at function exit
        try {
            bp::object res = fn_(args...);
            return bp::extract<TupleType>(res);
        } catch (bp::error_already_set const &) {
            std::string err = pyerr_to_string();
            GERROR(err.c_str());
            throw std::runtime_error(err);
        }
    }
};

/// PythonFunction for a single return type
template <typename RetType>
class PythonFunction<RetType> : public PythonFunctionBase
{
public:
    PythonFunction(const std::string& module, const std::string& funcname)
      : PythonFunctionBase(module, funcname)
    {
        // register the return type converter
        register_converter<RetType>();
    }

    template <typename... TS>
    RetType operator()(const TS&... args)
    {
        // register type converter for each parameter type
        register_converter<TS...>();
        GILLock lg; // lock GIL and release at function exit
        try {
            bp::object res = fn_(args...);
            return bp::extract<RetType>(res);
        } catch (bp::error_already_set const &) {
            std::string err = pyerr_to_string();
            GERROR(err.c_str());
            throw std::runtime_error(err);
        }
    }
};

/// PythonFunction for a single return type, special for bp::object type
template <>
class PythonFunction<bp::object> : public PythonFunctionBase
{
public:
    PythonFunction(const std::string& module, const std::string& funcname)
        : PythonFunctionBase(module, funcname)
    {
    }

    template <typename... TS>
    bp::object operator()(const TS&... args)
    {
        // register type converter for each parameter type
        register_converter<TS...>();
        GILLock lg; // lock GIL and release at function exit
        try {
            bp::object res = fn_(args...);
            return res;
        }
        catch (bp::error_already_set const &) {
            std::string err = pyerr_to_string();
            GERROR(err.c_str());
            throw std::runtime_error(err);
        }
    }
};


/// PythonFunction returning nothing
template <>
class PythonFunction<>  : public PythonFunctionBase
{
public:
    PythonFunction(const std::string& module, const std::string& funcname)
      : PythonFunctionBase(module, funcname) {}

    template <typename... TS>
    void operator()(const TS&... args)
    {
        // register type converter for each parameter type
        register_converter<TS...>();
        GILLock lg; // lock GIL and release at function exit
        try {
            bp::object res = fn_(args...);
        } catch (bp::error_already_set const &) {
            std::string err = pyerr_to_string();
            GERROR(err.c_str());
            throw std::runtime_error(err);
        }
    }
};

}

namespace boost { namespace python {
    EXPORTPYTHON bool hasattr(object o, const char* name);
} }

