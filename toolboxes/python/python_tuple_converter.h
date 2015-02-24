// Source: https://gist.github.com/niwibe/3729459
#ifndef GADGETRON_PYTHON_TUPLE_CONVERTER_H
#define GADGETRON_PYTHON_TUPLE_CONVERTER_H

#include <boost/python.hpp>
namespace bp = boost::python;

namespace Gadgetron {

/// indices trick
template<int ...> struct seq{};
template<int N, int ...S> struct gens : gens<N-1, N-1, S...>{};
template<int ...S> struct gens<0, S...> {typedef seq<S...> type;};

/// Used for expanding a C++ std::tuple into a boost::python::tuple
template <typename ...Args>
struct cpptuple2pytuple_wrapper {
    std::tuple<Args...> params;
    cpptuple2pytuple_wrapper(const std::tuple<Args...>& _params):params(_params){}

    bp::tuple delayed_dispatch() {
        return callFunc(typename gens<sizeof...(Args)>::type());
    }

    template<int ...S>
    bp::tuple callFunc(seq<S...>) {
        return bp::make_tuple(std::get<S>(params) ...);
    }
};

/// Used for expanding a boost::python::tuple into a C++ std::tuple
template <typename ...Args>
struct pytuple2cpptuple_wrapper {
    bp::tuple params;
    pytuple2cpptuple_wrapper(const bp::tuple& _params):params(_params){}

    std::tuple<Args...> delayed_dispatch() {
        return callFunc(typename gens<sizeof...(Args)>::type());
    }

    template<int ...S>
    std::tuple<Args...> callFunc(seq<S...>) {
        return std::make_tuple((static_cast<Args>(bp::extract<Args>(params[S])))...);
    }
};

/// Convert C++ std::tuple to boost::python::tuple as PyObject*.
template<typename ... Args> PyObject* cpptuple2pytuple(const std::tuple<Args...>& t) {
    cpptuple2pytuple_wrapper<Args...> wrapper(t);
    bp::tuple bpt = wrapper.delayed_dispatch();
    return bp::incref(bp::object(bpt).ptr());
}

/// Convert boost::python::tuple to C++ std::tuple.
template<typename ... Args> std::tuple<Args...> pytuple2cpptuple(PyObject* obj) {
    bp::tuple tup(bp::borrowed(obj));
    pytuple2cpptuple_wrapper<Args...> wrapper(tup);
    std::tuple<Args...> bpt = wrapper.delayed_dispatch();
    return bpt;
}

/// To-Python converter used by Boost
template<typename ... Args>
struct cpptuple_to_python_tuple {
    static PyObject* convert(const std::tuple<Args...>& t) {
        return cpptuple2pytuple<Args...>(t);
    }
};

/// From-Python converter used by Boost
template<typename ... Args>
struct cpptuple_from_python_tuple {
    cpptuple_from_python_tuple() {
        // actually register this converter
        bp::converter::registry::push_back(&convertible, &construct, bp::type_id<std::tuple<Args...> >());
    }

    /// Returns NULL if the bp::tuple is not convertible
    static void* convertible(PyObject* obj_ptr) {
        if (!PyTuple_CheckExact(obj_ptr)) {
            return NULL;
        }
        return obj_ptr;
    }

    /// Construct the std::tuple in place
    static void construct(PyObject* obj_ptr, bp::converter::rvalue_from_python_stage1_data* data) {
        void* storage = ((bp::converter::rvalue_from_python_storage<std::tuple<Args...> >*)data)->storage.bytes;
        // Use placement-new to make std::tuple in memory provided by Boost
        new (storage) std::tuple<Args...>(pytuple2cpptuple<Args...>(obj_ptr));
        data->convertible = storage;
    }
};

/// Create and register tuple converter as necessary
template <typename ...TS> void create_tuple_converter() {
    bp::type_info info = bp::type_id<std::tuple<TS...> >();
    const bp::converter::registration* reg = bp::converter::registry::query(info);
    // only register if not already registered!
    if (nullptr == reg || nullptr == (*reg).m_to_python) {
        bp::to_python_converter<std::tuple<TS...>, cpptuple_to_python_tuple<TS...> >();
        cpptuple_from_python_tuple<TS...>();
    }
}

/// Partial specialization of `python_converter` for std::tuple
template <typename ...TS>
struct python_converter<std::tuple<TS...> > {
    static void create()
    {
        // register tuple converter
        create_tuple_converter<TS...>();
        // register converter for each type in the tuple
        register_converter<TS...>();
    }
};

}

#endif // GADGETRON_PYTHON_TUPLE_CONVERTER_H
