#include "GadgetReference.h"
#include <boost/python.hpp>
//#include <numpy/arrayobject.h>

#include "../core/GadgetMRIHeaders.h"

#include "ismrmrd.h"

using namespace boost::python;

void acq_set_physiology_time_stamp(ISMRMRD::AcquisitionHeader &h, unsigned short i, uint32_t v)
{
	if (i < 3) {
		h.physiology_time_stamp[i] = v;
	}
}

uint32_t acq_get_physiology_time_stamp(ISMRMRD::AcquisitionHeader &h, unsigned short i)
{
	if (i < 3) {
		return h.physiology_time_stamp[i];
	}
	return 0;
}

void acq_set_channel_mask(ISMRMRD::AcquisitionHeader &h, unsigned short i, uint64_t v)
{
	if (i < 16) {
		h.channel_mask[i] = v;
	}
}

uint64_t acq_get_channel_mask(ISMRMRD::AcquisitionHeader &h, unsigned short i)
{
	if (i < 16) {
		return h.channel_mask[i];
	}
	return 0;
}

void acq_set_position(ISMRMRD::AcquisitionHeader &h, unsigned short i, float v)
{
	if (i < 3) {
		h.position[i] = v;
	}
}

float acq_get_position(ISMRMRD::AcquisitionHeader &h, unsigned short i)
{
	if (i < 3) {
		return h.position[i];
	}
	return 0.0f;
}

void acq_set_read_dir(ISMRMRD::AcquisitionHeader &h, unsigned short i, float v)
{
	if (i < 3) {
		h.read_dir[i] = v;
	}
}

float acq_get_read_dir(ISMRMRD::AcquisitionHeader &h, unsigned short i)
{
	if (i < 3) {
		return h.read_dir[i];
	}
	return 0.0f;
}

void acq_set_phase_dir(ISMRMRD::AcquisitionHeader &h, unsigned short i, float v)
{
	if (i < 3) {
		h.phase_dir[i] = v;
	}
}

float acq_get_phase_dir(ISMRMRD::AcquisitionHeader &h, unsigned short i)
{
	if (i < 3) {
		return h.phase_dir[i];
	}
	return 0.0f;
}

void acq_set_slice_dir(ISMRMRD::AcquisitionHeader &h, unsigned short i, float v)
{
	if (i < 3) {
		h.slice_dir[i] = v;
	}
}

float acq_get_slice_dir(ISMRMRD::AcquisitionHeader &h, unsigned short i)
{
	if (i < 3) {
		return h.slice_dir[i];
	}
	return 0.0f;
}

void acq_set_patient_table_position(ISMRMRD::AcquisitionHeader &h, unsigned short i, float v)
{
	if (i < 3) {
		h.patient_table_position[i] = v;
	}
}

float acq_get_patient_table_position(ISMRMRD::AcquisitionHeader &h, unsigned short i)
{
	if (i < 3) {
		return h.patient_table_position[i];
	}
	return 0.0f;
}

void acq_set_user_int(ISMRMRD::AcquisitionHeader &h, unsigned short i, int32_t v)
{
	if (i < 8) {
		h.user_int[i] = v;
	}
}

int32_t acq_get_user_int(ISMRMRD::AcquisitionHeader &h, unsigned short i)
{
	if (i < 8) {
		return h.user_int[i];
	}
	return 0;
}

void acq_set_user_float(ISMRMRD::AcquisitionHeader &h, unsigned short i, float v)
{
	if (i < 8) {
		h.user_float[i] = v;
	}
}

float acq_get_user_float(ISMRMRD::AcquisitionHeader &h, unsigned short i)
{
	if (i < 8) {
		return h.user_float[i];
	}
	return 0.0f;
}

void acq_set_encoding_user(ISMRMRD::EncodingCounters&e, unsigned short i, uint16_t v)
{
	if (i < 8) {
		e.user[i] = v;
	}
}

uint16_t acq_get_encoding_user(ISMRMRD::EncodingCounters&e, unsigned short i)
{
	if (i < 8) {
		return e.user[i];
	}
	return 0;
}

void img_set_matrix_size(ISMRMRD::ImageHeader &h, unsigned short i, uint16_t v)
{
	if (i < 3) {
		h.matrix_size[i] = v;
	}
}

uint16_t img_get_matrix_size(ISMRMRD::ImageHeader &h, unsigned short i)
{
	if (i < 3) {
		return h.matrix_size[i];
	}
	return 0;
}

void img_set_physiology_time_stamp(ISMRMRD::ImageHeader &h, unsigned short i, uint32_t v)
{
	if (i < 3) {
		h.physiology_time_stamp[i] = v;
	}
}

uint32_t img_get_physiology_time_stamp(ISMRMRD::ImageHeader &h, unsigned short i)
{
	if (i < 3) {
		return h.physiology_time_stamp[i];
	}
	return 0;
}

void img_set_position(ISMRMRD::ImageHeader &h, unsigned short i, float v)
{
	if (i < 3) {
		h.position[i] = v;
	}
}

float img_get_position(ISMRMRD::ImageHeader &h, unsigned short i)
{
	if (i < 3) {
		return h.position[i];
	}
	return 0.0f;
}

void img_set_read_dir(ISMRMRD::ImageHeader &h, unsigned short i, float v)
{
	if (i < 3) {
		h.read_dir[i] = v;
	}
}

float img_get_read_dir(ISMRMRD::ImageHeader &h, unsigned short i)
{
	if (i < 3) {
		return h.read_dir[i];
	}
	return 0.0f;
}

void img_set_phase_dir(ISMRMRD::ImageHeader &h, unsigned short i, float v)
{
	if (i < 3) {
		h.phase_dir[i] = v;
	}
}

float img_get_phase_dir(ISMRMRD::ImageHeader &h, unsigned short i)
{
	if (i < 3) {
		return h.phase_dir[i];
	}
	return 0.0f;
}

void img_set_slice_dir(ISMRMRD::ImageHeader &h, unsigned short i, float v)
{
	if (i < 3) {
		h.slice_dir[i] = v;
	}
}

float img_get_slice_dir(ISMRMRD::ImageHeader &h, unsigned short i)
{
	if (i < 3) {
		return h.slice_dir[i];
	}
	return 0.0f;
}

void img_set_patient_table_position(ISMRMRD::ImageHeader &h, unsigned short i, float v)
{
	if (i < 3) {
		h.patient_table_position[i] = v;
	}
}

float img_get_patient_table_position(ISMRMRD::ImageHeader &h, unsigned short i)
{
	if (i < 3) {
		return h.patient_table_position[i];
	}
	return 0.0f;
}


void img_set_user_int(ISMRMRD::ImageHeader &h, unsigned short i, int32_t v)
{
	if (i < 8) {
		h.user_int[i] = v;
	}
}

int32_t img_get_user_int(ISMRMRD::ImageHeader &h, unsigned short i)
{
	if (i < 8) {
		return h.user_int[i];
	}
	return 0;
}

void img_set_user_float(ISMRMRD::ImageHeader &h, unsigned short i, float v)
{
	if (i < 8) {
		h.user_float[i] = v;
	}
}

float img_get_user_float(ISMRMRD::ImageHeader &h, unsigned short i)
{
	if (i < 8) {
		return h.user_float[i];
	}
	return 0.0f;
}

void img_set_field_of_view(ISMRMRD::ImageHeader &h, unsigned short i, float v)
{
	if (i < 3) {
		h.field_of_view[i] = v;
	}
}

float img_get_field_of_view(ISMRMRD::ImageHeader &h, unsigned short i)
{
	if (i < 3) {
		return h.field_of_view[i];
	}
	return 0.0f;
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

	def("acq_set_physiology_time_stamp", acq_set_physiology_time_stamp);
	def("acq_get_physiology_time_stamp", acq_get_physiology_time_stamp);
	def("acq_set_channel_mask", acq_set_channel_mask);
	def("acq_get_channel_mask", acq_get_channel_mask);
	def("acq_set_position",acq_set_position);
	def("acq_get_position",acq_get_position);
	def("acq_set_read_dir",acq_set_read_dir);
	def("acq_get_read_dir",acq_get_read_dir);
	def("acq_set_phase_dir",acq_set_phase_dir);
	def("acq_get_phase_dir",acq_get_phase_dir);
	def("acq_set_slice_dir",acq_set_slice_dir);
	def("acq_get_slice_dir",acq_get_slice_dir);
	def("acq_set_patient_table_position", acq_set_patient_table_position);
	def("acq_get_patient_table_position", acq_get_patient_table_position);
	def("acq_set_user_int", acq_set_user_int);
	def("acq_get_user_int", acq_get_user_int);
	def("acq_set_user_float", acq_set_user_float);
	def("acq_get_user_float", acq_get_user_float);
	def("acq_set_encoding_user", acq_set_encoding_user);
	def("acq_get_encoding_user", acq_get_encoding_user);

	def("img_set_physiology_time_stamp", img_set_physiology_time_stamp);
	def("img_get_physiology_time_stamp", img_get_physiology_time_stamp);
	def("img_set_position",img_set_position);
	def("img_get_position",img_get_position);
	def("img_set_read_dir",img_set_read_dir);
	def("img_get_read_dir",img_get_read_dir);
	def("img_set_phase_dir",img_set_phase_dir);
	def("img_get_phase_dir",img_get_phase_dir);
	def("img_set_slice_dir",img_set_slice_dir);
	def("img_get_slice_dir",img_get_slice_dir);
	def("img_set_patient_table_position", img_set_patient_table_position);
	def("img_get_patient_table_position", img_get_patient_table_position);
	def("img_set_user_int", img_set_user_int);
	def("img_get_user_int", img_get_user_int);
	def("img_set_user_float", img_set_user_float);
	def("img_get_user_float", img_get_user_float);
	def("img_get_field_of_view", img_get_field_of_view);
	def("img_set_field_of_view", img_set_field_of_view);
	def("img_get_matrix_size", img_get_matrix_size);
	def("img_set_matrix_size", img_set_matrix_size);

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


	class_<ISMRMRD::ImageHeader>("ImageHeader")
			.def_readwrite("flags", &ISMRMRD::ImageHeader::flags)
			.def_readwrite("channels", &ISMRMRD::ImageHeader::channels)
			.def_readwrite("slice", &ISMRMRD::ImageHeader::slice)
			.def_readwrite("contrast", &ISMRMRD::ImageHeader::contrast)
			.def_readwrite("set", &ISMRMRD::ImageHeader::set)
			.def_readwrite("phase", &ISMRMRD::ImageHeader::phase)
			.def_readwrite("average", &ISMRMRD::ImageHeader::average)
			.def_readwrite("repetition", &ISMRMRD::ImageHeader::repetition)
			.def_readwrite("acquisition_time_stamp", &ISMRMRD::ImageHeader::acquisition_time_stamp)
			.def_readwrite("image_data_type", &ISMRMRD::ImageHeader::image_data_type)
			.def_readwrite("image_type", &ISMRMRD::ImageHeader::image_type)
			.def_readwrite("image_index", &ISMRMRD::ImageHeader::image_index)
			.def_readwrite("image_series_index", &ISMRMRD::ImageHeader::image_series_index)
			;

	class_<Gadgetron::GadgetReference>("GadgetReference")
    		.def("return_acquisition", &Gadgetron::GadgetReference::return_data<ISMRMRD::AcquisitionHeader>)
    		.def("return_image", &Gadgetron::GadgetReference::return_data<ISMRMRD::ImageHeader>)

    		;

	enum_<ISMRMRD::ImageDataType>("ImageDataType")
    		   .value("DATA_COMPLEX_FLOAT", ISMRMRD::DATA_COMPLEX_FLOAT)
    		   .value("DATA_FLOAT", ISMRMRD::DATA_FLOAT)
    		   .value("DATA_UNSIGNED_SHORT", ISMRMRD::DATA_UNSIGNED_SHORT)
    		   ;


	enum_<ISMRMRD::ImageType>("ImageType")
				  .value("TYPE_MAGNITUDE",ISMRMRD::TYPE_MAGNITUDE)
				  .value("TYPE_PHASE", ISMRMRD::TYPE_PHASE)
				  .value("TYPE_REAL",ISMRMRD::TYPE_REAL)
				  .value("TYPE_IMAG",ISMRMRD::TYPE_IMAG)
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
