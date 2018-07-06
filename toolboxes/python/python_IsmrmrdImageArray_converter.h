/**
\file   python_IsmrmrdImageArray_converter.h
\brief  Implementation converter from/to python for IsmrmrdImageArray data structure
\author Hui Xue
*/

#pragma once
#include "python_toolbox.h"
#include "python_numpy_wrappers.h"

#include "hoNDArray.h"
#include "mri_core_data.h"
#include "log.h"

#include <boost/python.hpp>
namespace bp = boost::python;

namespace Gadgetron {

    class IsmrmrdImageArray_to_python_object
    {
    public:
        static PyObject* convert(const IsmrmrdImageArray & arrayData)
        {
            GILLock lock;
            bp::object pygadgetron = bp::import("gadgetron");

            auto data = bp::object(arrayData.data_);
            /*auto pyHeaders = bp::list();

            size_t n;
            for (n=0; n<arrayData.headers_.get_number_of_elements(); n++)
            {
                auto headers = boost::python::object(arrayData.headers_(n));
                pyHeaders.append(headers);
            }

            auto pyMeta = bp::list();
            for (n = 0; n<arrayData.meta_.size(); n++)
            {
                auto meta = boost::python::object(arrayData.meta_[n]);
                pyMeta.append(meta);
            }*/

            auto pyHeaders = boost::python::object(arrayData.headers_);
            auto pyMeta = boost::python::object(arrayData.meta_);
            auto pyWav = arrayData.waveform_ ? boost::python::object(*arrayData.waveform_) : boost::python::object();
            auto pyAcqHeaders = arrayData.acq_headers_ ? boost::python::object(*arrayData.acq_headers_) : boost::python::object();

            bp::incref(data.ptr());
            bp::incref(pyHeaders.ptr());
            bp::incref(pyMeta.ptr());
            bp::incref(pyWav.ptr());
            bp::incref(pyAcqHeaders.ptr());
            auto buffer = pygadgetron.attr("IsmrmrdImageArray")(data, pyHeaders, pyMeta, pyWav, pyAcqHeaders);

            // increment the reference count so it exists after `return`
            return bp::incref(buffer.ptr());
        }
    };

    // ------------------------------------------------------------------------
    struct IsmrmrdImageArray_from_python_object
    {
        IsmrmrdImageArray_from_python_object()
        {
            // actually register this converter with Boost
            bp::converter::registry::push_back(
                &convertible,
                &construct,
                bp::type_id<IsmrmrdImageArray>());
        }

        /// Returns NULL if the object is not convertible. Or well.... it should
        static void* convertible(PyObject* obj)
        {
            return obj;
        }

        /// Construct an hoNDArray in-place
        static void construct(PyObject* obj, bp::converter::rvalue_from_python_stage1_data* data)
        {

            void* storage = ((bp::converter::rvalue_from_python_storage<IsmrmrdImageArray >*)data)->storage.bytes;
            IsmrmrdImageArray* reconData = new (storage) IsmrmrdImageArray;
            data->convertible = storage;

            try
            {
                bp::object pyImageArray((bp::handle<>(bp::borrowed(obj))));

                reconData->data_ = bp::extract<hoNDArray<std::complex<float>>>(pyImageArray.attr("data"));
                reconData->headers_ = bp::extract<hoNDArray<ISMRMRD::ImageHeader>>(pyImageArray.attr("headers"));
                reconData->meta_ = bp::extract<std::vector<ISMRMRD::MetaContainer>>(pyImageArray.attr("meta"));

                if (PyObject_HasAttrString(pyImageArray.ptr(), "waveform"))
                    reconData->waveform_ = bp::extract<std::vector<ISMRMRD::Waveform>>(pyImageArray.attr("waveform"));

                if (PyObject_HasAttrString(pyImageArray.ptr(), "acq_headers"))
                    reconData->acq_headers_ = bp::extract<hoNDArray<ISMRMRD::AcquisitionHeader>>(pyImageArray.attr("acq_headers"));
            }
            catch (const bp::error_already_set&)
            {
                std::string err = pyerr_to_string();
                GERROR(err.c_str());
                throw std::runtime_error(err);
            }
        }
    };

    // ------------------------------------------------------------------------
    template<> struct python_converter<IsmrmrdImageArray>
    {
        static void create()
        {
            bp::type_info info = bp::type_id<IsmrmrdImageArray >();
            const bp::converter::registration* reg = bp::converter::registry::query(info);
            // only register if not already registered!
            if (nullptr == reg || nullptr == (*reg).m_to_python)
            {
                bp::to_python_converter<IsmrmrdImageArray, IsmrmrdImageArray_to_python_object >();
                IsmrmrdImageArray_from_python_object();
            }

        }
    };
}
