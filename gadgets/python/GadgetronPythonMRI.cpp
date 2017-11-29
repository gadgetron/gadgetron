#include "GadgetReference.h"
#include "GadgetInstrumentationStreamController.h"
#include <boost/python.hpp>
#include "../mri_core/GadgetMRIHeaders.h"
#include <ismrmrd/ismrmrd.h>

using namespace boost::python;

BOOST_PYTHON_MODULE(GadgetronPythonMRI)
{
    // boost::python::numeric::array::set_module_and_type("numpy", "ndarray");

    class_<Gadgetron::GadgetReference>("GadgetReference")
      .def("return_acquisition", &Gadgetron::GadgetReference::return_acquisition)
      .def("return_recondata",&Gadgetron::GadgetReference::return_recondata)
      .def("return_image_cplx", &Gadgetron::GadgetReference::return_image_cplx)
      .def("return_image_cplx_attr", &Gadgetron::GadgetReference::return_image_cplx_attr)
      .def("return_image_float", &Gadgetron::GadgetReference::return_image_float)
      .def("return_image_float_attr", &Gadgetron::GadgetReference::return_image_float_attr)
      .def("return_image_ushort", &Gadgetron::GadgetReference::return_image_ushort)
      .def("return_image_ushort_attr", &Gadgetron::GadgetReference::return_image_ushort_attr)
      ;

    class_<Gadgetron::GadgetInstrumentationStreamControllerWrapper>("GadgetInstrumentationStreamController")
      .def("put_config", &Gadgetron::GadgetInstrumentationStreamControllerWrapper::put_config)
      .def("put_acquisition", &Gadgetron::GadgetInstrumentationStreamControllerWrapper::put_acquisition)
      .def("put_recondata", &Gadgetron::GadgetInstrumentationStreamControllerWrapper::put_recondata)
      .def("put_image_cplx", &Gadgetron::GadgetInstrumentationStreamControllerWrapper::put_image_cplx)
      .def("put_image_cplx_attr", &Gadgetron::GadgetInstrumentationStreamControllerWrapper::put_image_cplx_attr)
      .def("put_image_float", &Gadgetron::GadgetInstrumentationStreamControllerWrapper::put_image_float)
      .def("put_image_float_attr", &Gadgetron::GadgetInstrumentationStreamControllerWrapper::put_image_float_attr)
      .def("put_image_ushort", &Gadgetron::GadgetInstrumentationStreamControllerWrapper::put_image_ushort)
      .def("put_image_ushort_attr", &Gadgetron::GadgetInstrumentationStreamControllerWrapper::put_image_ushort_attr)
      .def("prepend_gadget", &Gadgetron::GadgetInstrumentationStreamControllerWrapper::prepend_gadget)
      .def("close", &Gadgetron::GadgetInstrumentationStreamControllerWrapper::close)
      .def("set_python_gadget", &Gadgetron::GadgetInstrumentationStreamControllerWrapper::set_python_gadget)
      .def("set_parameter", &Gadgetron::GadgetInstrumentationStreamControllerWrapper::set_parameter)
      ;    
}
