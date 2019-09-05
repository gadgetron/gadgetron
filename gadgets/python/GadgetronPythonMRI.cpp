#include "GadgetReference.h"
#include <boost/python.hpp>
#include "../mri_core/GadgetMRIHeaders.h"

#include <ismrmrd/ismrmrd.h>

using namespace boost::python;

BOOST_PYTHON_MODULE(GadgetronPythonMRI)
{
    class_<Gadgetron::GadgetReference>("GadgetReference", no_init)
      .def("return_acquisition", &Gadgetron::GadgetReference::return_acquisition)
      .def("return_recondata",&Gadgetron::GadgetReference::return_recondata)
      .def("return_ismrmrd_image_array", &Gadgetron::GadgetReference::return_ismrmrd_image_array)
      .def("return_image_cplx", &Gadgetron::GadgetReference::return_image_cplx)
      .def("return_image_cplx_attr", &Gadgetron::GadgetReference::return_image_cplx_attr)
      .def("return_image_float", &Gadgetron::GadgetReference::return_image_float)
      .def("return_image_float_attr", &Gadgetron::GadgetReference::return_image_float_attr)
      .def("return_image_ushort", &Gadgetron::GadgetReference::return_image_ushort)
      .def("return_image_ushort_attr", &Gadgetron::GadgetReference::return_image_ushort_attr)
      ;


}
