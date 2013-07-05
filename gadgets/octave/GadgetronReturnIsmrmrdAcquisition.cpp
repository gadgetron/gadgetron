#include <octave/oct.h>
#include <octave/ov-struct.h>

#include "GadgetContainerMessage.h"
#include "ismrmrd.h"
#include "hoNDArray.h"
#include "OctaveCommunicator.h"

using namespace Gadgetron;

DEFUN_DLD (GadgetronReturnIsmrmrdAcquisition, args, nargout,
	   "GadgetronReturnIsmrmrdAcquisition Returns Acquisition to the Gadgetron")
{
  int nargin = args.length ();

  octave_value retval;
     
  if (nargin != 3) {
    print_usage(); 
  } else {
    std::string id(args(0).string_value());
    Octave_map h(args(1).map_value());
    FloatComplexNDArray d(args(2).complex_array_value());

    GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1 =
    		new GadgetContainerMessage<ISMRMRD::AcquisitionHeader>();

    ISMRMRD::AcquisitionHeader* head = m1->getObjectPtr();

    head->version = octave_value(h.contents("version")(0)).uint16_scalar_value();
    head->flags = octave_value(h.contents("flags")(0)).uint64_scalar_value();
    head->measurement_uid = octave_value(h.contents("measurement_uid")(0)).uint32_scalar_value();
    head->scan_counter = octave_value(h.contents("scan_counter")(0)).uint32_scalar_value();
    head->acquisition_time_stamp = octave_value(h.contents("acquisition_time_stamp")(0)).uint32_scalar_value();
    head->measurement_uid = octave_value(h.contents("measurement_uid")(0)).uint32_scalar_value();
    head->physiology_time_stamp[0] = octave_value(h.contents("physiology_time_stamp")(0)).uint32_array_value()(0);
    head->physiology_time_stamp[1] = octave_value(h.contents("physiology_time_stamp")(0)).uint32_array_value()(1);
    head->physiology_time_stamp[2] = octave_value(h.contents("physiology_time_stamp")(0)).uint32_array_value()(2);
    head->number_of_samples = octave_value(h.contents("number_of_samples")(0)).uint16_scalar_value();
    head->available_channels = octave_value(h.contents("available_channels")(0)).uint16_scalar_value();
    head->active_channels = octave_value(h.contents("active_channels")(0)).uint16_scalar_value();
    head->channel_mask[0] = octave_value(h.contents("channel_mask")(0)).uint64_array_value()(0);
    head->channel_mask[1] = octave_value(h.contents("channel_mask")(0)).uint64_array_value()(1);
    head->channel_mask[2] = octave_value(h.contents("channel_mask")(0)).uint64_array_value()(2);
    head->channel_mask[3] = octave_value(h.contents("channel_mask")(0)).uint64_array_value()(3);
    head->channel_mask[4] = octave_value(h.contents("channel_mask")(0)).uint64_array_value()(4);
    head->channel_mask[5] = octave_value(h.contents("channel_mask")(0)).uint64_array_value()(5);
    head->channel_mask[6] = octave_value(h.contents("channel_mask")(0)).uint64_array_value()(6);
    head->channel_mask[7] = octave_value(h.contents("channel_mask")(0)).uint64_array_value()(7);
    head->channel_mask[8] = octave_value(h.contents("channel_mask")(0)).uint64_array_value()(8);
    head->channel_mask[9] = octave_value(h.contents("channel_mask")(0)).uint64_array_value()(9);
    head->channel_mask[10] = octave_value(h.contents("channel_mask")(0)).uint64_array_value()(10);
    head->channel_mask[11] = octave_value(h.contents("channel_mask")(0)).uint64_array_value()(11);
    head->channel_mask[12] = octave_value(h.contents("channel_mask")(0)).uint64_array_value()(12);
    head->channel_mask[13] = octave_value(h.contents("channel_mask")(0)).uint64_array_value()(13);
    head->channel_mask[14] = octave_value(h.contents("channel_mask")(0)).uint64_array_value()(14);
    head->channel_mask[15] = octave_value(h.contents("channel_mask")(0)).uint64_array_value()(15);
    head->discard_pre = octave_value(h.contents("discard_pre")(0)).uint16_scalar_value();
    head->discard_post = octave_value(h.contents("discard_post")(0)).uint16_scalar_value();
    head->center_sample = octave_value(h.contents("center_sample")(0)).uint16_scalar_value();
    head->encoding_space_ref = octave_value(h.contents("encoding_space_ref")(0)).uint16_scalar_value();
    head->trajectory_dimensions = octave_value(h.contents("trajectory_dimensions")(0)).uint16_scalar_value();
    head->sample_time_us = octave_value(h.contents("sample_time_us")(0)).float_scalar_value();
    head->position[0] = octave_value(h.contents("position")(0)).float_array_value()(0);
    head->position[1] = octave_value(h.contents("position")(0)).float_array_value()(1);
    head->position[2] = octave_value(h.contents("position")(0)).float_array_value()(2);
    head->read_dir[0] = octave_value(h.contents("read_dir")(0)).float_array_value()(0);
    head->read_dir[1] = octave_value(h.contents("read_dir")(0)).float_array_value()(1);
    head->read_dir[2] = octave_value(h.contents("read_dir")(0)).float_array_value()(2);
    head->phase_dir[0] = octave_value(h.contents("phase_dir")(0)).float_array_value()(0);
    head->phase_dir[1] = octave_value(h.contents("phase_dir")(0)).float_array_value()(1);
    head->phase_dir[2] = octave_value(h.contents("phase_dir")(0)).float_array_value()(2);
    head->slice_dir[0] = octave_value(h.contents("read_dir")(0)).float_array_value()(0);
    head->slice_dir[1] = octave_value(h.contents("read_dir")(0)).float_array_value()(1);
    head->slice_dir[2] = octave_value(h.contents("read_dir")(0)).float_array_value()(2);
    head->patient_table_position[0] = octave_value(h.contents("patient_table_position")(0)).float_array_value()(0);
    head->patient_table_position[1] = octave_value(h.contents("patient_table_position")(0)).float_array_value()(1);
    head->patient_table_position[2] = octave_value(h.contents("patient_table_position")(0)).float_array_value()(2);
    head->idx.kspace_encode_step_1 = octave_value(octave_value(h.contents("idx")(0)).map_value().contents("kspace_encode_step_1")(0)).uint16_scalar_value();
    head->idx.kspace_encode_step_2 = octave_value(octave_value(h.contents("idx")(0)).map_value().contents("kspace_encode_step_2")(0)).uint16_scalar_value();
    head->idx.average              = octave_value(octave_value(h.contents("idx")(0)).map_value().contents("average")(0)).uint16_scalar_value();
    head->idx.slice                = octave_value(octave_value(h.contents("idx")(0)).map_value().contents("slice")(0)).uint16_scalar_value();
    head->idx.contrast             = octave_value(octave_value(h.contents("idx")(0)).map_value().contents("contrast")(0)).uint16_scalar_value();
    head->idx.phase                = octave_value(octave_value(h.contents("idx")(0)).map_value().contents("phase")(0)).uint16_scalar_value();
    head->idx.repetition           = octave_value(octave_value(h.contents("idx")(0)).map_value().contents("repetition")(0)).uint16_scalar_value();
    head->idx.set                  = octave_value(octave_value(h.contents("idx")(0)).map_value().contents("set")(0)).uint16_scalar_value();
    head->idx.segment              = octave_value(octave_value(h.contents("idx")(0)).map_value().contents("segment")(0)).uint16_scalar_value();
    head->idx.user[0]              = octave_value(octave_value(h.contents("idx")(0)).map_value().contents("user")(0)).uint16_array_value()(0);
    head->idx.user[1]              = octave_value(octave_value(h.contents("idx")(0)).map_value().contents("user")(0)).uint16_array_value()(1);
    head->idx.user[2]              = octave_value(octave_value(h.contents("idx")(0)).map_value().contents("user")(0)).uint16_array_value()(2);
    head->idx.user[3]              = octave_value(octave_value(h.contents("idx")(0)).map_value().contents("user")(0)).uint16_array_value()(3);
    head->idx.user[4]              = octave_value(octave_value(h.contents("idx")(0)).map_value().contents("user")(0)).uint16_array_value()(4);
    head->idx.user[5]              = octave_value(octave_value(h.contents("idx")(0)).map_value().contents("user")(0)).uint16_array_value()(5);
    head->idx.user[6]              = octave_value(octave_value(h.contents("idx")(0)).map_value().contents("user")(0)).uint16_array_value()(6);
    head->idx.user[7]              = octave_value(octave_value(h.contents("idx")(0)).map_value().contents("user")(0)).uint16_array_value()(7);
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
