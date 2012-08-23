#include <boost/python.hpp>
#include <numpy/arrayobject.h>

#include "../core/GadgetMRIHeaders.h"
#include "GadgetReference.h"
#include "ismrmrd.h"

using namespace boost::python;

void set_physiology_time_stamp(ISMRMRD::AcquisitionHeader &h, unsigned short i, uint32_t v)
{
	if (i < 3) {
		h.physiology_time_stamp[i] = v;
	}
}

uint32_t get_physiology_time_stamp(ISMRMRD::AcquisitionHeader &h, unsigned short i)
{
	if (i < 3) {
		return h.physiology_time_stamp[i];
	}
	return 0;
}

void set_channel_mask(ISMRMRD::AcquisitionHeader &h, unsigned short i, uint64_t v)
{
	if (i < 16) {
		h.channel_mask[i] = v;
	}
}

uint64_t get_channel_mask(ISMRMRD::AcquisitionHeader &h, unsigned short i)
{
	if (i < 16) {
		return h.channel_mask[i];
	}
	return 0;
}

void set_position(ISMRMRD::AcquisitionHeader &h, unsigned short i, float v)
{
	if (i < 3) {
		h.position[i] = v;
	}
}

float get_position(ISMRMRD::AcquisitionHeader &h, unsigned short i)
{
	if (i < 3) {
		return h.position[i];
	}
	return 0.0f;
}

void set_quaternion(ISMRMRD::AcquisitionHeader &h, unsigned short i, float v)
{
	if (i < 4) {
		h.quaternion[i] = v;
	}
}

float get_quaternion(ISMRMRD::AcquisitionHeader &h, unsigned short i)
{
	if (i < 4) {
		return h.quaternion[i];
	}
	return 0.0f;
}

void set_patient_table_position(ISMRMRD::AcquisitionHeader &h, unsigned short i, float v)
{
	if (i < 3) {
		h.patient_table_position[i] = v;
	}
}

float get_patient_table_position(ISMRMRD::AcquisitionHeader &h, unsigned short i)
{
	if (i < 3) {
		return h.patient_table_position[i];
	}
	return 0.0f;
}

void set_user_int(ISMRMRD::AcquisitionHeader &h, unsigned short i, int32_t v)
{
	if (i < 8) {
		h.user_int[i] = v;
	}
}

int32_t get_user_int(ISMRMRD::AcquisitionHeader &h, unsigned short i)
{
	if (i < 8) {
		return h.user_int[i];
	}
	return 0;
}

void set_user_float(ISMRMRD::AcquisitionHeader &h, unsigned short i, float v)
{
	if (i < 8) {
		h.user_float[i] = v;
	}
}

float get_user_float(ISMRMRD::AcquisitionHeader &h, unsigned short i)
{
	if (i < 8) {
		return h.user_float[i];
	}
	return 0.0f;
}

void set_encoding_user(ISMRMRD::EncodingCounters&e, unsigned short i, uint16_t v)
{
	if (i < 8) {
		e.user[i] = v;
	}
}

uint16_t get_encoding_user(ISMRMRD::EncodingCounters&e, unsigned short i)
{
	if (i < 8) {
		return e.user[i];
	}
	return 0;
}

BOOST_PYTHON_MODULE(GadgetronPythonMRI)
{

	//import_array();
	boost::python::numeric::array::set_module_and_type("numpy", "ndarray");


	class_<ISMRMRD::EncodingCounters>("EncodingCounters")
			.def_readwrite("kspace_encode_step_1",   &ISMRMRD::EncodingCounters::kspace_encode_step_1)
			.def_readwrite("kspace_encode_step_2",   &ISMRMRD::EncodingCounters::kspace_encode_step_2)
			.def_readwrite("average",                &ISMRMRD::EncodingCounters::average)
			.def_readwrite("slice",                  &ISMRMRD::EncodingCounters::slice)
			.def_readwrite("contrast",               &ISMRMRD::EncodingCounters::contrast)
			.def_readwrite("phase",                  &ISMRMRD::EncodingCounters::phase)
			.def_readwrite("repetition",             &ISMRMRD::EncodingCounters::repetition)
			.def_readwrite("segment",                &ISMRMRD::EncodingCounters::segment)
			;

	def("set_physiology_time_stamp", set_physiology_time_stamp);
	def("get_physiology_time_stamp", get_physiology_time_stamp);
	def("set_channel_mask", set_channel_mask);
	def("get_channel_mask", get_channel_mask);
	def("set_position",set_position);
	def("get_position",get_position);
	def("set_quaternion",set_quaternion);
	def("get_quaternion",get_quaternion);
	def("set_patient_table_position", set_patient_table_position);
	def("get_patient_table_position", get_patient_table_position);
	def("set_user_int", set_user_int);
	def("get_user_int", get_user_int);
	def("set_user_float", set_user_float);
	def("get_user_float", get_user_float);
	def("set_encoding_user", set_encoding_user);
	def("get_encoding_user", get_encoding_user);

	class_<ISMRMRD::AcquisitionHeader>("AcquisitionHeader")
			.def_readwrite("version",                &ISMRMRD::AcquisitionHeader::version)
			.def_readwrite("flags",                  &ISMRMRD::AcquisitionHeader::flags)
			.def_readwrite("measurement_uid",        &ISMRMRD::AcquisitionHeader::measurement_uid)
			.def_readwrite("scan_counter",           &ISMRMRD::AcquisitionHeader::scan_counter)
			.def_readwrite("acquisition_time_stamp", &ISMRMRD::AcquisitionHeader::acquisition_time_stamp)
			.def_readwrite("number_of_samples",      &ISMRMRD::AcquisitionHeader::number_of_samples)
			.def_readwrite("available_channels",     &ISMRMRD::AcquisitionHeader::available_channels)
			.def_readwrite("active_channels",        &ISMRMRD::AcquisitionHeader::active_channels)
			.def_readwrite("discard_pre",            &ISMRMRD::AcquisitionHeader::discard_pre)
			.def_readwrite("discard_post",           &ISMRMRD::AcquisitionHeader::discard_post)
			.def_readwrite("centre_sample",          &ISMRMRD::AcquisitionHeader::center_sample)
			.def_readwrite("encoding_space_ref",     &ISMRMRD::AcquisitionHeader::encoding_space_ref)
			.def_readwrite("trajectory_dimensions",  &ISMRMRD::AcquisitionHeader::trajectory_dimensions)
			.def_readwrite("sample_time_us",         &ISMRMRD::AcquisitionHeader::sample_time_us)
			.def_readwrite("idx",                    &ISMRMRD::AcquisitionHeader::idx)
			;


	class_<GadgetMessageImage>("GadgetMessageImage")
			.def_readwrite("flags", &GadgetMessageImage::flags)
			.def("get_matrix_size", &GadgetMessageImage::get_matrix_size)
			.def("set_matrix_size", &GadgetMessageImage::set_matrix_size)
			.def_readwrite("channels", &GadgetMessageImage::channels)
			.def("get_position", &GadgetMessageImage::get_position)
			.def("set_position", &GadgetMessageImage::set_position)
			.def("get_quaternion", &GadgetMessageImage::get_quaternion)
			.def("set_quaternion", &GadgetMessageImage::set_quaternion)
			.def_readwrite("table_position", &GadgetMessageImage::table_position)
			.def_readwrite("slice", &GadgetMessageImage::slice)
			.def_readwrite("contrast", &GadgetMessageImage::contrast)
			.def_readwrite("set", &GadgetMessageImage::set)
			.def_readwrite("phase", &GadgetMessageImage::phase)
			.def_readwrite("average", &GadgetMessageImage::average)
			.def_readwrite("repetition", &GadgetMessageImage::repetition)
			.def_readwrite("time_stamp", &GadgetMessageImage::time_stamp)
			.def_readwrite("pmu_time_stamp", &GadgetMessageImage::pmu_time_stamp)
			.def_readwrite("image_format", &GadgetMessageImage::image_format)
			.def_readwrite("image_type", &GadgetMessageImage::image_type)
			.def_readwrite("image_index", &GadgetMessageImage::image_index)
			.def_readwrite("image_series_index", &GadgetMessageImage::image_series_index)
			;

	class_<GadgetReference>("GadgetReference")
    		.def("return_acquisition", &GadgetReference::return_data<ISMRMRD::AcquisitionHeader>)
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
