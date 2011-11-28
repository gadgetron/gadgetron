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
    .def_readwrite("table_position", &GadgetMessageAcquisition::table_position)
    .def_readwrite("idx", &GadgetMessageAcquisition::idx)
    .def_readwrite("min_idx", &GadgetMessageAcquisition::min_idx)
    .def_readwrite("max_idx", &GadgetMessageAcquisition::max_idx)

    ;

  class_<GadgetMessageImage>("GadgetMessageImage")
	.def_readwrite("flags", &GadgetMessageImage::flags)
    .def("get_matrix_size", &GadgetMessageImage::get_matrix_size)
    .def("set_matrix_size", &GadgetMessageImage::set_matrix_size)
    .def_readwrite("channels", &GadgetMessageImage::channels)
    .def("get_position", &GadgetMessageImage::get_position)
    .def("set_position", &GadgetMessageImage::set_position)
    .def("get_quarternion", &GadgetMessageImage::get_quarternion)
    .def("set_quarternion", &GadgetMessageImage::set_quarternion)
    .def_readwrite("table_position", &GadgetMessageImage::table_position)
    .def_readwrite("data_idx_min", &GadgetMessageImage::data_idx_min)
    .def_readwrite("data_idx_max", &GadgetMessageImage::data_idx_max)
    .def_readwrite("data_idx_current", &GadgetMessageImage::data_idx_current)
    .def_readwrite("time_stamp", &GadgetMessageImage::time_stamp)
    .def_readwrite("image_format", &GadgetMessageImage::image_format)
    .def_readwrite("image_type", &GadgetMessageImage::image_type)
    .def_readwrite("image_index", &GadgetMessageImage::image_index)
    .def_readwrite("image_series_index", &GadgetMessageImage::image_series_index)

    ;

  class_<GadgetReference>("GadgetReference")
    .def("return_acquisition", &GadgetReference::return_data<GadgetMessageAcquisition>)
    .def("return_image", &GadgetReference::return_data<GadgetMessageImage>)

    ;

  enum_<GadgetImageFormats>("GadgetImageFormats")
       .value("GADGET_IMAGE_COMPLEX_FLOAT", GADGET_IMAGE_COMPLEX_FLOAT)
       .value("GADGET_IMAGE_REAL_FLOAT", GADGET_IMAGE_REAL_FLOAT)
       .value("GADGET_IMAGE_REAL_UNSIGNED_SHORT", GADGET_IMAGE_REAL_UNSIGNED_SHORT)
       ;


  enum_<GadgetImageTypes>("GadgetImageTypes")
		  .value("GADGET_IMAGE_MAGNITUDE",GADGET_IMAGE_MAGNITUDE)
		  .value("GADGET_IMAGE_PHASE", GADGET_IMAGE_PHASE)
		  .value("GADGET_IMAGE_REAL",GADGET_IMAGE_REAL)
		  .value("GADGET_IMAGE_IMAG",GADGET_IMAGE_IMAG)
		  ;

  enum_<GadgetMessageID>("GadgetMessageID")
		  .value("GADGET_MESSAGE_EXT_ID_MIN",GADGET_MESSAGE_EXT_ID_MIN)
		  .value("GADGET_MESSAGE_ACQUISITION",GADGET_MESSAGE_ACQUISITION)
		  .value("GADGET_MESSAGE_NEW_MEASUREMENT",GADGET_MESSAGE_NEW_MEASUREMENT)
		  .value("GADGET_MESSAGE_END_OF_SCAN",GADGET_MESSAGE_END_OF_SCAN)
		  .value("GADGET_MESSAGE_IMAGE_CPLX_FLOAT",GADGET_MESSAGE_IMAGE_CPLX_FLOAT)
		  .value("GADGET_MESSAGE_IMAGE_REAL_FLOAT",GADGET_MESSAGE_IMAGE_REAL_FLOAT)
		  .value("GADGET_MESSAGE_IMAGE_REAL_USHORT",GADGET_MESSAGE_IMAGE_REAL_USHORT)
		  .value("GADGET_MESSAGE_EMPTY",GADGET_MESSAGE_EMPTY)
		  .value("GADGET_MESSAGE_EXT_ID_MAX",GADGET_MESSAGE_EXT_ID_MAX)
		  ;
}
