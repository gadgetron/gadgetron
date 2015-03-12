#include "GadgetReference.h"
#include "GadgetInstrumentationStreamController.h"
#include <boost/python.hpp>
#include "../mri_core/GadgetMRIHeaders.h"
#include <ismrmrd/ismrmrd.h>

using namespace boost::python;

BOOST_PYTHON_MODULE(GadgetronPythonMRI)
{
    boost::python::numeric::array::set_module_and_type("numpy", "ndarray");

    class_<Gadgetron::GadgetReference>("GadgetReference")
      .def("return_acquisition", &Gadgetron::GadgetReference::return_data<ISMRMRD::AcquisitionHeader>)
      .def("return_image", &Gadgetron::GadgetReference::return_data<ISMRMRD::ImageHeader>)
      ;

    class_<Gadgetron::GadgetInstrumentationStreamControllerWrapper>("GadgetInstrumentationStreamController")
      .def("put_config", &Gadgetron::GadgetInstrumentationStreamControllerWrapper::put_config)
      .def("put_acquisition", &Gadgetron::GadgetInstrumentationStreamControllerWrapper::put_acquisition)
      .def("put_image", &Gadgetron::GadgetInstrumentationStreamControllerWrapper::put_image)
      .def("prepend_gadget", &Gadgetron::GadgetInstrumentationStreamControllerWrapper::prepend_gadget)
      .def("close", &Gadgetron::GadgetInstrumentationStreamControllerWrapper::close)
      .def("is_closed", &Gadgetron::GadgetInstrumentationStreamControllerWrapper::is_closed)
      .def("set_python_gadget", &Gadgetron::GadgetInstrumentationStreamControllerWrapper::set_python_gadget)
      ;
    
    enum_<Gadgetron::GadgetMessageID>("GadgetMessageID")
      .value("GADGET_MESSAGE_EXT_ID_MIN",Gadgetron::GADGET_MESSAGE_EXT_ID_MIN)
      .value("GADGET_MESSAGE_ACQUISITION",Gadgetron::GADGET_MESSAGE_ACQUISITION)
      .value("GADGET_MESSAGE_NEW_MEASUREMENT",Gadgetron::GADGET_MESSAGE_NEW_MEASUREMENT)
      .value("GADGET_MESSAGE_END_OF_SCAN",Gadgetron::GADGET_MESSAGE_END_OF_SCAN)
      .value("GADGET_MESSAGE_IMAGE_CPLX_FLOAT",Gadgetron::GADGET_MESSAGE_IMAGE_CPLX_FLOAT)
      .value("GADGET_MESSAGE_IMAGE_REAL_FLOAT",Gadgetron::GADGET_MESSAGE_IMAGE_REAL_FLOAT)
      .value("GADGET_MESSAGE_IMAGE_REAL_USHORT",Gadgetron::GADGET_MESSAGE_IMAGE_REAL_USHORT)
      .value("GADGET_MESSAGE_ISMRMRD_ACQUISITION", Gadgetron::GADGET_MESSAGE_ISMRMRD_ACQUISITION)
      .value("GADGET_MESSAGE_ISMRMRD_IMAGE_CPLX_FLOAT", Gadgetron::GADGET_MESSAGE_ISMRMRD_IMAGE_CPLX_FLOAT)
      .value("GADGET_MESSAGE_ISMRMRD_IMAGE_REAL_FLOAT", Gadgetron::GADGET_MESSAGE_ISMRMRD_IMAGE_REAL_FLOAT)
      .value("GADGET_MESSAGE_ISMRMRD_IMAGE_REAL_USHORT", Gadgetron::GADGET_MESSAGE_ISMRMRD_IMAGE_REAL_USHORT)
      .value("GADGET_MESSAGE_EMPTY",Gadgetron::GADGET_MESSAGE_EMPTY)
      .value("GADGET_MESSAGE_EXT_ID_MAX",Gadgetron::GADGET_MESSAGE_EXT_ID_MAX)
      ;
}
