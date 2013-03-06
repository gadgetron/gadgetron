#include "OctaveGadget.h"


 int AcquisitionOctaveGadget::process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1,
	      GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2)
 {

   //We want to avoid a deadlock for the Python GIL if this python call results in an output that the GadgetReference will not be able to get rid of.
   //This is kind of a nasty busy wait, maybe we should add an event handler to the NotificationStrategy of the Q or something, but for now, this will do it.
   while (this->next()->msg_queue()->is_full()) {
     //GADGET_DEBUG2("Gadget (%s) sleeping while downstream Gadget (%s) does some work\n", this->module()->name(), this->next()->module()->name());
     ACE_Time_Value tv(0,10000); //Sleep for 10ms while the downstream Gadget does some work
      ACE_OS::sleep(tv);
   }


	Octave_map m;
	ISMRMRD::AcquisitionHeader h = *m1->getObjectPtr();

	m.assign("version",                    h.version);
	m.assign("flags",                      h.flags);
	m.assign("measurement_uid",            h.measurement_uid);
	m.assign("scan_counter",               h.scan_counter);
	m.assign("acquisition_time_stamp",     h.acquisition_time_stamp);

	dim_vector d(1); d(0) = 3;
	uint32NDArray phys_time(d);
	memcpy(&phys_time(0),h.physiology_time_stamp,sizeof(uint32_t)*3);
	m.assign("physiology_time_stamp",            octave_value(phys_time));


	m.assign("number_of_samples",                h.number_of_samples);
	m.assign("available_channels",               h.available_channels);
	m.assign("active_channels",                  h.active_channels);

	d(0) = 16;
	uint64NDArray channel_mask(d);
	memcpy(&channel_mask(0),h.channel_mask,sizeof(uint64_t)*16);
	m.assign("channel_mask",                     octave_value(channel_mask));

	m.assign("discard_pre",                      h.discard_pre);
	m.assign("discard_post",                     h.discard_post);
	m.assign("center_sample",                    h.center_sample);
	m.assign("encoding_space_ref",               h.encoding_space_ref);
	m.assign("trajectory_dimensions",            h.trajectory_dimensions);
	m.assign("sample_time_us",                   h.sample_time_us);

	d(0) = 3;
	FloatNDArray position(d);
	memcpy(&position(0),h.position,sizeof(float)*3);
	m.assign("position",                         octave_value(position));

	d(0) = 3;
	FloatNDArray read_dir(d);
	memcpy(&read_dir(0),h.read_dir,sizeof(float)*3);
	m.assign("read_dir",            octave_value(read_dir));

	d(0) = 3;
	FloatNDArray phase_dir(d);
	memcpy(&phase_dir(0),h.phase_dir,sizeof(float)*3);
	m.assign("phase_dir",            octave_value(phase_dir));

	d(0) = 3;
	FloatNDArray slice_dir(d);
	memcpy(&slice_dir(0),h.slice_dir,sizeof(float)*3);
	m.assign("slice_dir",            octave_value(slice_dir));

	d(0) = 3;
	FloatNDArray patient_table_position(d);
	memcpy(&patient_table_position(0),h.patient_table_position,sizeof(float)*3);
	m.assign("patient_table_position",         octave_value(patient_table_position));

	Octave_map idx;

	idx.assign("kspace_encode_step_1",       h.idx.kspace_encode_step_1);
	idx.assign("kspace_encode_step_2",       h.idx.kspace_encode_step_2);
	idx.assign("average",                    h.idx.average);
	idx.assign("slice",                      h.idx.slice);
	idx.assign("contrast",                   h.idx.contrast);
	idx.assign("phase",                      h.idx.phase);
	idx.assign("repetition",                 h.idx.phase);
	idx.assign("phase",                      h.idx.repetition);
	idx.assign("set",                        h.idx.set);
	idx.assign("segment",                    h.idx.segment);

	d(0) = 8;
	uint16NDArray user(d);
	memcpy(&user(0),h.idx.user,sizeof(uint16_t)*8);
	idx.assign("user",                    octave_value(user));
	m.assign("idx",                         octave_value(idx));

	d(0) = 8;
	int32NDArray user_int(d);
	memcpy(&user_int(0),h.user_int,sizeof(int32_t)*8);
	m.assign("user_int",                         octave_value(user_int));

	d(0) = 8;
	FloatNDArray user_float(d);
	memcpy(&user_float(0),h.user_float,sizeof(float)*8);
	m.assign("user_float",                         octave_value(user_float));

	//Make a copy of the data for sending to Octave.
    dim_vector dims;
    for (unsigned int i =0; i < m2->getObjectPtr()->get_number_of_dimensions(); i++) {
    	dims(i) = m2->getObjectPtr()->get_size(i);
    }
    FloatComplexNDArray data(dims);
    memcpy(data.fortran_vec(),m2->getObjectPtr()->get_data_ptr(),sizeof(float)*2*data.nelem());

    octave_value_list in;
    in(0) = m;
    in(1) = data;

    octave_value_list out = OctaveCommunicator::instance()->octave_feval (datafunc_->c_str(), in, 2);

    //We are now done with the data
    m1->release();

    return GADGET_OK;
 }

 int ImageOctaveGadget::process(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1,
	      GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2)
 {

   //We want to avoid a deadlock for the Python GIL if this python call results in an output that the GadgetReference will not be able to get rid of.
   //This is kind of a nasty busy wait, maybe we should add an event handler to the NotificationStrategy of the Q or something, but for now, this will do it.
   while (this->next()->msg_queue()->is_full()) {
     //GADGET_DEBUG2("Gadget (%s) sleeping while downstream Gadget (%s) does some work\n", this->module()->name(), this->next()->module()->name());
     ACE_Time_Value tv(0,10000); //Sleep for 10ms while the downstream Gadget does some work
     ACE_OS::sleep(tv);
   }


	Octave_map m;
	ISMRMRD::ImageHeader h = *m1->getObjectPtr();

	m.assign("version",                    h.version);
	m.assign("flags",                      h.flags);
	m.assign("measurement_uid",            h.measurement_uid);

	dim_vector d(1);
	d(0) = 3;
	uint16NDArray matrix_size(d);
	memcpy(&matrix_size(0),h.matrix_size,sizeof(uint16_t)*3);
	m.assign("matrix_size",            octave_value(matrix_size));

	d(0) = 3;
	FloatNDArray field_of_view(d);
	memcpy(&field_of_view(0),h.field_of_view,sizeof(float)*3);
	m.assign("field_of_view",            octave_value(field_of_view));

	m.assign("channels",                    h.channels);

	d(0) = 3;
	FloatNDArray position(d);
	memcpy(&position(0),h.position,sizeof(float)*3);
	m.assign("position",            octave_value(position));

	d(0) = 3;
	FloatNDArray read_dir(d);
	memcpy(&read_dir(0),h.read_dir,sizeof(float)*3);
	m.assign("read_dir",            octave_value(read_dir));

	d(0) = 3;
	FloatNDArray phase_dir(d);
	memcpy(&phase_dir(0),h.phase_dir,sizeof(float)*3);
	m.assign("phase_dir",            octave_value(phase_dir));

	d(0) = 3;
	FloatNDArray slice_dir(d);
	memcpy(&slice_dir(0),h.slice_dir,sizeof(float)*3);
	m.assign("slice_dir",            octave_value(slice_dir));

	d(0) = 3;
	FloatNDArray patient_table_position(d);
	memcpy(&patient_table_position(0),h.patient_table_position,sizeof(float)*3);
	m.assign("patient_table_position",            octave_value(patient_table_position));

	m.assign("average",                    h.average);
	m.assign("slice",                    h.slice);
	m.assign("contrast",                    h.contrast);
	m.assign("phase",                    h.phase);
	m.assign("repetition",                    h.repetition);
	m.assign("set",                    h.set);

	d(0) = 3;
	uint32NDArray physiology_time_stamp(d);
	memcpy(&physiology_time_stamp(0),h.physiology_time_stamp,sizeof(uint32_t)*3);
	m.assign("physiology_time_stamp",            octave_value(physiology_time_stamp));

	m.assign("image_data_type",                    h.image_data_type);
	m.assign("image_type",                    h.image_type);

	m.assign("image_index",                    h.image_index);
	m.assign("image_series_index",            h.image_series_index);

	d(0) = 8;
	int32NDArray user_int(d);
	memcpy(&user_int(0),h.user_int,sizeof(int32_t)*8);
	m.assign("user_int",                         octave_value(user_int));

	d(0) = 8;
	FloatNDArray user_float(d);
	memcpy(&user_float(0),h.user_float,sizeof(float)*8);
	m.assign("user_float",                         octave_value(user_float));

    dim_vector dims;
    for (unsigned int i =0; i < m2->getObjectPtr()->get_number_of_dimensions(); i++) {
    	dims(i) = m2->getObjectPtr()->get_size(i);
    }

    FloatComplexNDArray data(dims);
    memcpy(&data(0),m2->getObjectPtr()->get_data_ptr(),sizeof(float)*2*data.nelem());

    octave_value_list in;
    in(0) = m; //octave_value (this->next()->module()->name());
    in(1) = data;

    octave_value_list out = OctaveCommunicator::instance()->octave_feval (datafunc_->c_str(), in, 2);

    m1->release();

    return GADGET_OK;
 }


GADGET_FACTORY_DECLARE(AcquisitionOctaveGadget)
GADGET_FACTORY_DECLARE(ImageOctaveGadget)
