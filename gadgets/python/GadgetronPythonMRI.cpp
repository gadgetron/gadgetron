#include <boost/python.hpp>
#include <numpy/arrayobject.h>

#include "../core/GadgetMRIHeaders.h"
#include "GadgetReference.h"

using namespace boost::python;

BOOST_PYTHON_MODULE(GadgetronPythonMRI)
{

  //import_array();
  boost::python::numeric::array::set_module_and_type("numpy", "ndarray");

  class_<LoopCounters>("LoopCounters")
   .def_readwrite("line", &LoopCounters::line)
   .def_readwrite("acquisition",&LoopCounters::acquisition)
   .def_readwrite("slice",&LoopCounters::slice)
   .def_readwrite("partition",&LoopCounters::partition)
   .def_readwrite("echo",&LoopCounters::echo)
   .def_readwrite("phase",&LoopCounters::phase)
   .def_readwrite("repetition",&LoopCounters::repetition)
   .def_readwrite("set",&LoopCounters::set)
   .def_readwrite("segment",&LoopCounters::segment)
   .def_readwrite("channel",&LoopCounters::channel)
    ;

  class_<GadgetMessageAcquisition>("GadgetMessageAcquisition")
    .def_readwrite("flags", &GadgetMessageAcquisition::flags)
    .def_readwrite("meas_uid", &GadgetMessageAcquisition::meas_uid)
    .def_readwrite("scan_counter", &GadgetMessageAcquisition::scan_counter)
    .def_readwrite("time_stamp", &GadgetMessageAcquisition::time_stamp)
    .def_readwrite("samples", &GadgetMessageAcquisition::samples)
    .def_readwrite("channels", &GadgetMessageAcquisition::channels)
    .def("get_position", &GadgetMessageAcquisition::get_position)
    .def("set_position", &GadgetMessageAcquisition::set_position)
    .def("get_quarternion", &GadgetMessageAcquisition::get_quarternion)
    .def("set_quarternion", &GadgetMessageAcquisition::set_quarternion)
    .def_readwrite("idx", &GadgetMessageAcquisition::idx)
    .def_readwrite("min_idx", &GadgetMessageAcquisition::min_idx)
    .def_readwrite("max_idx", &GadgetMessageAcquisition::max_idx)

    ;

  class_<GadgetMessageImage>("GadgetMessageImage")
    .def("get_matrix_size", &GadgetMessageImage::get_matrix_size)
    .def("set_matrix_size", &GadgetMessageImage::set_matrix_size)
    .def_readwrite("channels", &GadgetMessageImage::channels)
    .def("get_position", &GadgetMessageImage::get_position)
    .def("set_position", &GadgetMessageImage::set_position)
    .def("get_quarternion", &GadgetMessageImage::get_quarternion)
    .def("set_quarternion", &GadgetMessageImage::set_quarternion)
    .def_readwrite("data_idx_min", &GadgetMessageImage::data_idx_min)
    .def_readwrite("data_idx_max", &GadgetMessageImage::data_idx_max)
    .def_readwrite("data_idx_current", &GadgetMessageImage::data_idx_current)
    .def_readwrite("time_stamp", &GadgetMessageImage::time_stamp)

    ;

  class_<GadgetReference>("GadgetReference")
    .def("return_acquisition", &GadgetReference::return_data<GadgetMessageAcquisition>)
    .def("return_image", &GadgetReference::return_data<GadgetMessageImage>)

    ;

}
