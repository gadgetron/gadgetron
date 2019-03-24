#ifndef GADGETINSTRUMENTATIONSTREAMCONTROLLER_H
#define GADGETINSTRUMENTATIONSTREAMCONTROLLER_H

#include "GadgetStreamInterface.h"
#include "Gadget.h"
#include <ismrmrd/ismrmrd.h>
#include <ismrmrd/waveform.h>
#include "python_toolbox.h"
#include <boost/python.hpp>

namespace Gadgetron {

    class GadgetInstrumentationStreamController
        : public GadgetStreamInterface
    {
    public:
        GadgetInstrumentationStreamController();
        int open();
        int close();
        int prepend_gadget(const char* gadgetname,
            const char* dllname,
            const char* classname);

        virtual ~GadgetInstrumentationStreamController();

        template<class TH, class TD> int put_data(TH header, boost::python::object arr, const char* meta = 0);
        int put_config(const char* config);
        int put_acquisition(ISMRMRD::AcquisitionHeader acq, boost::python::object arr, const char* meta = 0);
        int put_image_cplx(ISMRMRD::ImageHeader img, boost::python::object arr, const char* meta = 0);
        int put_image_float(ISMRMRD::ImageHeader img, boost::python::object arr, const char* meta = 0);
        int put_image_ushort(ISMRMRD::ImageHeader img, boost::python::object arr, const char* meta = 0);
        int put_recondata(boost::python::object rec);
        int put_ismrmrd_image_array(boost::python::object rec);
        int set_python_gadget(boost::python::object g)
        {
            python_gadget_ = g;
            boost::python::incref(python_gadget_.ptr());
            return GADGET_OK;
        }

        virtual int output_ready(ACE_Message_Block* mb);
        void set_parameter(const char* gadgetname, const char* parameter, const char* value);

    protected:
        boost::python::object python_gadget_;
        template <class T1, class T2, class T3> int return_data(ACE_Message_Block* mb);
        int return_recondata(ACE_Message_Block* mb);
        int return_ismrmrd_image_array(ACE_Message_Block* mb);
    };

    class GadgetInstrumentationStreamControllerWrapper
    {
    public:
        GadgetInstrumentationStreamControllerWrapper()
        {
            // ensure boost can convert between hoNDArrays and NumPy arrays automatically
            register_converter<hoNDArray<std::complex<float> > >();
            register_converter<hoNDArray< float > >();
            register_converter<hoNDArray< uint32_t > >();
            register_converter<hoNDArray< unsigned short > >();
            // ensure boost can convert ISMRMRD headers automatically
            register_converter<ISMRMRD::ImageHeader>();
            register_converter<ISMRMRD::AcquisitionHeader>();
            register_converter<ISMRMRD::ISMRMRD_WaveformHeader>();
            register_converter<ISMRMRD::Waveform>();
            register_converter<ISMRMRD::MetaContainer>();
            // ensure other types are converted
            register_converter<hoNDArray<ISMRMRD::AcquisitionHeader>>();
            register_converter<hoNDArray<ISMRMRD::ImageHeader> >();
            register_converter<std::vector<ISMRMRD::MetaContainer> >();
            register_converter<std::vector<ISMRMRD::Waveform> >();

            register_converter<IsmrmrdReconData>();
            register_converter<IsmrmrdImageArray>();

            cntrl_ = new GadgetInstrumentationStreamController;
        }

        ~GadgetInstrumentationStreamControllerWrapper()
        {
            delete cntrl_;
        }

        int prepend_gadget(const char* gadgetname,
            const char* dllname,
            const char* classname)
        {
            return cntrl_->prepend_gadget(gadgetname, dllname, classname);
        }

        int put_config(const char* config)
        {
            return cntrl_->put_config(config);
        }

        int put_acquisition(ISMRMRD::AcquisitionHeader acq, boost::python::object arr)
        {
            return cntrl_->put_acquisition(acq, arr);
        }


        int put_image_cplx(ISMRMRD::ImageHeader img, boost::python::object arr)
        {
            return cntrl_->put_image_cplx(img, arr);
        }

        int put_image_cplx_attr(ISMRMRD::ImageHeader img, boost::python::object arr, const char* meta = 0)
        {
            return cntrl_->put_image_cplx(img, arr, meta);
        }

        int put_image_float(ISMRMRD::ImageHeader img, boost::python::object arr)
        {
            return cntrl_->put_image_float(img, arr);
        }

        int put_image_float_attr(ISMRMRD::ImageHeader img, boost::python::object arr, const char* meta = 0)
        {
            return cntrl_->put_image_float(img, arr, meta);
        }

        int put_image_ushort(ISMRMRD::ImageHeader img, boost::python::object arr)
        {
            return cntrl_->put_image_ushort(img, arr);
        }

        int put_image_ushort_attr(ISMRMRD::ImageHeader img, boost::python::object arr, const char* meta = 0)
        {
            return cntrl_->put_image_ushort(img, arr, meta);
        }

        int put_recondata(boost::python::object rec) {
            return cntrl_->put_recondata(rec);
        }

        int put_ismrmrd_image_array(boost::python::object rec) {
            return cntrl_->put_ismrmrd_image_array(rec);
        }

        int close()
        {
            // allow other threads to finish returning data to Python
            Py_BEGIN_ALLOW_THREADS;
            cntrl_->close();
            Py_END_ALLOW_THREADS;
            return 0;
        }

        int set_python_gadget(boost::python::object g)
        {
            return cntrl_->set_python_gadget(g);
        }

        void set_parameter(const char* gadgetname, const char* parameter, const char* value)
        {
            cntrl_->set_parameter(gadgetname, parameter, value);
        }

    protected:
        GadgetInstrumentationStreamController* cntrl_;
    };

}
#endif //GADGETINSTRUMENTATIONSTREAMCONTROLLER_H
