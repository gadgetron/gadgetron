#include <octave/oct.h>
#include <octave/ov-struct.h>

#include "OctaveCommunicator.h"
#include "ismrmrd.h"
#include "hoNDArray.h"

using namespace Gadgetron;

DEFUN_DLD (GadgetronReturnIsmrmrdImage, args, nargout,
	   "GadgetronReturnIsmrmrdImage return Image to the Gadgetron")
{
	 int nargin = args.length ();

	  octave_value retval;

	  if (nargin != 3) {
	    print_usage();
	  } else {
	    std::string id(args(0).string_value());
	    Octave_map h(args(1).map_value());
	    FloatComplexNDArray d(args(2).complex_array_value());

	    GadgetContainerMessage<ISMRMRD::ImageHeader>* m1 =
	    		new GadgetContainerMessage<ISMRMRD::ImageHeader>();

	    ISMRMRD::ImageHeader* head = m1->getObjectPtr();

	    head->version = octave_value(h.contents("version")(0)).uint16_scalar_value();
	    head->flags = octave_value(h.contents("flags")(0)).uint64_scalar_value();
	    head->measurement_uid = octave_value(h.contents("measurement_uid")(0)).uint32_scalar_value();
	    head->matrix_size[0] = octave_value(h.contents("matrix_size")(0)).uint16_array_value()(0);
	    head->matrix_size[1] = octave_value(h.contents("matrix_size")(0)).uint16_array_value()(1);
	    head->matrix_size[2] = octave_value(h.contents("matrix_size")(0)).uint16_array_value()(2);
		head->field_of_view[0] = octave_value(h.contents("field_of_view")(0)).float_array_value()(0);
	    head->field_of_view[1] = octave_value(h.contents("field_of_view")(0)).float_array_value()(1);
	    head->field_of_view[2] = octave_value(h.contents("field_of_view")(0)).float_array_value()(2);
		head->channels = octave_value(h.contents("channels")(0)).uint16_scalar_value();
	    head->position[0] = octave_value(h.contents("position")(0)).float_array_value()(0);
	    head->position[1] = octave_value(h.contents("position")(0)).float_array_value()(1);
	    head->position[2] = octave_value(h.contents("position")(0)).float_array_value()(2);
	    head->read_dir[0] = octave_value(h.contents("read_dir")(0)).float_array_value()(0);
	    head->read_dir[1] = octave_value(h.contents("read_dir")(0)).float_array_value()(1);
	    head->read_dir[2] = octave_value(h.contents("read_dir")(0)).float_array_value()(2);
	    head->phase_dir[0] = octave_value(h.contents("phase_dir")(0)).float_array_value()(0);
	    head->phase_dir[1] = octave_value(h.contents("phase_dir")(0)).float_array_value()(1);
	    head->phase_dir[2] = octave_value(h.contents("phase_dir")(0)).float_array_value()(2);
	    head->slice_dir[0] = octave_value(h.contents("slice_dir")(0)).float_array_value()(0);
	    head->slice_dir[1] = octave_value(h.contents("slice_dir")(0)).float_array_value()(1);
	    head->slice_dir[2] = octave_value(h.contents("slice_dir")(0)).float_array_value()(2);
	    head->patient_table_position[0] = octave_value(h.contents("patient_table_position")(0)).float_array_value()(0);
	    head->patient_table_position[1] = octave_value(h.contents("patient_table_position")(0)).float_array_value()(1);
	    head->patient_table_position[2] = octave_value(h.contents("patient_table_position")(0)).float_array_value()(2);
	    head->average = octave_value(h.contents("average")(0)).uint16_scalar_value();
	    head->slice = octave_value(h.contents("slice")(0)).uint16_scalar_value();
	    head->contrast = octave_value(h.contents("contrast")(0)).uint16_scalar_value();
	    head->phase = octave_value(h.contents("phase")(0)).uint16_scalar_value();
	    head->repetition = octave_value(h.contents("repetition")(0)).uint16_scalar_value();
	    head->set = octave_value(h.contents("set")(0)).uint16_scalar_value();
	    head->acquisition_time_stamp = octave_value(h.contents("acquisition_time_stamp")(0)).uint32_scalar_value();
	    head->physiology_time_stamp[0] = octave_value(h.contents("physiology_time_stamp")(0)).uint32_array_value()(0);
	    head->physiology_time_stamp[1] = octave_value(h.contents("physiology_time_stamp")(0)).uint32_array_value()(1);
	    head->physiology_time_stamp[2] = octave_value(h.contents("physiology_time_stamp")(0)).uint32_array_value()(2);
	    head->image_data_type          = octave_value(h.contents("image_data_type")(0)).uint16_scalar_value();
	    head->image_data_type          = octave_value(h.contents("image_data_type")(0)).uint16_scalar_value();
	    head->image_data_type          = octave_value(h.contents("image_data_type")(0)).uint16_scalar_value();
	    head->image_data_type          = octave_value(h.contents("image_data_type")(0)).uint16_scalar_value();
	    head->user_int[0]              = octave_value(h.contents("user_int")(0)).int32_array_value()(0);
	    head->user_int[1]              = octave_value(h.contents("user_int")(0)).int32_array_value()(1);
	    head->user_int[2]              = octave_value(h.contents("user_int")(0)).int32_array_value()(2);
	    head->user_int[3]              = octave_value(h.contents("user_int")(0)).int32_array_value()(3);
	    head->user_int[4]              = octave_value(h.contents("user_int")(0)).int32_array_value()(4);
	    head->user_int[5]              = octave_value(h.contents("user_int")(0)).int32_array_value()(5);
	    head->user_int[6]              = octave_value(h.contents("user_int")(0)).int32_array_value()(6);
	    head->user_int[7]              = octave_value(h.contents("user_int")(0)).int32_array_value()(7);
	    head->user_float[0]            = octave_value(h.contents("user_float")(0)).int32_array_value()(0);
	    head->user_float[1]            = octave_value(h.contents("user_float")(0)).int32_array_value()(1);
	    head->user_float[2]            = octave_value(h.contents("user_float")(0)).int32_array_value()(2);
	    head->user_float[3]            = octave_value(h.contents("user_float")(0)).int32_array_value()(3);
	    head->user_float[4]            = octave_value(h.contents("user_float")(0)).int32_array_value()(4);
	    head->user_float[5]            = octave_value(h.contents("user_float")(0)).int32_array_value()(5);
	    head->user_float[6]            = octave_value(h.contents("user_float")(0)).int32_array_value()(6);
	    head->user_float[7]            = octave_value(h.contents("user_float")(0)).int32_array_value()(7);

	    GadgetContainerMessage< hoNDArray<std::complex<float> > >* m2 =
	    		new GadgetContainerMessage< hoNDArray<std::complex<float> > >();

	    std::vector<unsigned int> dims;
	    for (unsigned int i = 0; i < d.dims().length(); i++) {
	    	dims.push_back(d.dims()(i));
	    }

	    try {
	        m2->getObjectPtr()->create(&dims);
	    } catch (...) {
	    	GADGET_DEBUG1("Failed to allocate return array\n");
	    	m1->release();
	    }

	    memcpy(m2->getObjectPtr()->get_data_ptr(), &d(0), sizeof(float)*2*d.nelem());

	    m1->cont(m2);

	    OctaveCommunicator::instance()->message_gadget(id, m1);
	  }

	  return octave_value_list ();
}
