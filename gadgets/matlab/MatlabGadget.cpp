#include "MatlabGadget.h"

#define NUMBER_OF_ACQ_FIELDS    24
#define NUMBER_OF_IDX_FIELDS    10
#define NUMBER_OF_IMG_FIELDS    25

const char* ismrmrd_acq_field_names[] = {
    "version",
    "flags",
    "measurement_uid",
    "scan_counter",
    "acquisition_time_stamp",
    "physiology_time_stamp",
    "number_of_samples",
    "available_channels",
    "active_channels",
    "channel_mask",
    "discard_pre",
    "discard_post",
    "center_sample",
    "encoding_space_ref",
    "trajectory_dimensions",
    "sample_time_us",
    "position",
    "read_dir",
    "phase_dir",
    "slice_dir",
    "patient_table_position",
    "idx",
    "user_int",
    "user_float"
};

const char* ismrmrd_idx_field_names[] = {
    "kspace_encode_step_1",
    "kspace_encode_step_2",
    "average",
    "slice",
    "contrast",
    "phase",
    "repetition",
    "set",
    "segment",
    "user"
};

const char* ismrmrd_img_field_names[] = {
    "version",
    "flags",
    "measurement_uid",
    "matrix_size",
    "field_of_view",
    "channels",
    "position",
    "read_dir",
    "phase_dir",
    "slice_dir",
    "patient_table_position",
    "average",
    "slice",
    "contrast",
    "phase",
    "repetition",
    "set",
    "acquisition_time_stamp",
    "physiology_time_stamp",
    "image_data_type",
    "image_type",
    "image_index",
    "image_series_index",
    "user_int",
    "user_float"
};

int AcquisitionMatlabGadget::process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1,
        GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2)
{
    ISMRMRD::AcquisitionHeader *acq = m1->getObjectPtr();

    // Create struct array for storing a single ISMRMRD Encoding Counters
    mwSize idx_dims[2] = {1, 1};
    mxArray *idx = mxCreateStructArray(2, idx_dims, NUMBER_OF_IDX_FIELDS, ismrmrd_idx_field_names);
    mxArray *idx_values[NUMBER_OF_IDX_FIELDS];

    idx_values[0] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
    memcpy(mxGetData(idx_values[0]), &acq->idx.kspace_encode_step_1, mxGetElementSize(idx_values[0]));

    idx_values[1] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
    memcpy(mxGetData(idx_values[1]), &acq->idx.kspace_encode_step_2, mxGetElementSize(idx_values[1]));

    idx_values[2] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
    memcpy(mxGetData(idx_values[2]), &acq->idx.average, mxGetElementSize(idx_values[2]));

    idx_values[3] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
    memcpy(mxGetData(idx_values[3]), &acq->idx.slice, mxGetElementSize(idx_values[3]));

    idx_values[4] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
    memcpy(mxGetData(idx_values[4]), &acq->idx.contrast, mxGetElementSize(idx_values[4]));

    idx_values[5] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
    memcpy(mxGetData(idx_values[5]), &acq->idx.phase, mxGetElementSize(idx_values[5]));

    idx_values[6] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
    memcpy(mxGetData(idx_values[6]), &acq->idx.repetition, mxGetElementSize(idx_values[6]));

    idx_values[7] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
    memcpy(mxGetData(idx_values[7]), &acq->idx.set, mxGetElementSize(idx_values[7]));

    idx_values[8] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
    memcpy(mxGetData(idx_values[8]), &acq->idx.segment, mxGetElementSize(idx_values[8]));

    idx_values[9] = mxCreateNumericMatrix(ISMRMRD_USER_INTS, 1, mxUINT16_CLASS, mxREAL);
    memcpy(mxGetData(idx_values[9]), &acq->idx.user, ISMRMRD_USER_INTS * mxGetElementSize(idx_values[9]));

    // Set the fields in the ISMRMRD Acquisition Header struct mxArray
    for (mwIndex i = 0; i < NUMBER_OF_IDX_FIELDS; i++) {
        mxSetFieldByNumber(idx, 0, i, idx_values[i]);
    }


    // Create struct array for storing a single ISMRMRD Acquisition Header
    mwSize acq_dims[2] = {1, 1};
    mxArray *acqhdr = mxCreateStructArray(2, acq_dims, NUMBER_OF_ACQ_FIELDS, ismrmrd_acq_field_names);

    mxArray *acq_values[NUMBER_OF_ACQ_FIELDS];

    acq_values[0] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
    memcpy(mxGetData(acq_values[0]), &acq->version, mxGetElementSize(acq_values[0]));

    acq_values[1] = mxCreateNumericMatrix(1, 1, mxUINT64_CLASS, mxREAL);
    memcpy(mxGetData(acq_values[1]), &acq->flags, mxGetElementSize(acq_values[1]));

    acq_values[2] = mxCreateNumericMatrix(1, 1, mxUINT32_CLASS, mxREAL);
    memcpy(mxGetData(acq_values[2]), &acq->measurement_uid, mxGetElementSize(acq_values[2]));

    acq_values[3] = mxCreateNumericMatrix(1, 1, mxUINT32_CLASS, mxREAL);
    memcpy(mxGetData(acq_values[3]), &acq->scan_counter, mxGetElementSize(acq_values[3]));

    acq_values[4] = mxCreateNumericMatrix(1, 1, mxUINT32_CLASS, mxREAL);
    memcpy(mxGetData(acq_values[4]), &acq->acquisition_time_stamp, mxGetElementSize(acq_values[4]));

    acq_values[5] = mxCreateNumericMatrix(ISMRMRD_PHYS_STAMPS, 1, mxUINT32_CLASS, mxREAL);
    memcpy(mxGetData(acq_values[5]), &acq->physiology_time_stamp, ISMRMRD_PHYS_STAMPS * mxGetElementSize(acq_values[5]));

    acq_values[6] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
    memcpy(mxGetData(acq_values[6]), &acq->number_of_samples, mxGetElementSize(acq_values[6]));

    acq_values[7] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
    memcpy(mxGetData(acq_values[7]), &acq->available_channels, mxGetElementSize(acq_values[7]));

    acq_values[8] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
    memcpy(mxGetData(acq_values[8]), &acq->active_channels, mxGetElementSize(acq_values[8]));

    acq_values[9] = mxCreateNumericMatrix(ISMRMRD_CHANNEL_MASKS, 1, mxUINT64_CLASS, mxREAL);
    memcpy(mxGetData(acq_values[9]), &acq->channel_mask, ISMRMRD_CHANNEL_MASKS * mxGetElementSize(acq_values[9]));

    acq_values[10] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
    memcpy(mxGetData(acq_values[10]), &acq->discard_pre, mxGetElementSize(acq_values[10]));

    acq_values[11] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
    memcpy(mxGetData(acq_values[11]), &acq->discard_post, mxGetElementSize(acq_values[11]));

    acq_values[12] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
    memcpy(mxGetData(acq_values[12]), &acq->center_sample, mxGetElementSize(acq_values[12]));

    acq_values[13] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
    memcpy(mxGetData(acq_values[13]), &acq->encoding_space_ref, mxGetElementSize(acq_values[13]));

    acq_values[14] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
    memcpy(mxGetData(acq_values[14]), &acq->trajectory_dimensions, mxGetElementSize(acq_values[14]));

    acq_values[15] = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
    memcpy(mxGetData(acq_values[15]), &acq->sample_time_us, mxGetElementSize(acq_values[15]));

    acq_values[16] = mxCreateNumericMatrix(ISMRMRD_POSITION_LENGTH, 1, mxSINGLE_CLASS, mxREAL);
    memcpy(mxGetData(acq_values[16]), &acq->position, ISMRMRD_POSITION_LENGTH * mxGetElementSize(acq_values[16]));

    acq_values[17] = mxCreateNumericMatrix(ISMRMRD_DIRECTION_LENGTH, 1, mxSINGLE_CLASS, mxREAL);
    memcpy(mxGetData(acq_values[17]), &acq->read_dir, ISMRMRD_DIRECTION_LENGTH * mxGetElementSize(acq_values[17]));

    acq_values[18] = mxCreateNumericMatrix(ISMRMRD_DIRECTION_LENGTH, 1, mxSINGLE_CLASS, mxREAL);
    memcpy(mxGetData(acq_values[18]), &acq->phase_dir, ISMRMRD_DIRECTION_LENGTH * mxGetElementSize(acq_values[18]));

    acq_values[19] = mxCreateNumericMatrix(ISMRMRD_DIRECTION_LENGTH, 1, mxSINGLE_CLASS, mxREAL);
    memcpy(mxGetData(acq_values[19]), &acq->slice_dir, ISMRMRD_DIRECTION_LENGTH * mxGetElementSize(acq_values[19]));

    acq_values[20] = mxCreateNumericMatrix(ISMRMRD_POSITION_LENGTH, 1, mxSINGLE_CLASS, mxREAL);
    memcpy(mxGetData(acq_values[20]), &acq->patient_table_position, ISMRMRD_POSITION_LENGTH * mxGetElementSize(acq_values[20]));

    acq_values[21] = idx;

    acq_values[22] = mxCreateNumericMatrix(ISMRMRD_USER_INTS, 1, mxINT32_CLASS, mxREAL);
    memcpy(mxGetData(acq_values[22]), &acq->user_int, ISMRMRD_USER_INTS * mxGetElementSize(acq_values[22]));

    acq_values[23] = mxCreateNumericMatrix(ISMRMRD_USER_FLOATS, 1, mxSINGLE_CLASS, mxREAL);
    memcpy(mxGetData(acq_values[23]), &acq->user_float, ISMRMRD_USER_FLOATS * mxGetElementSize(acq_values[23]));

    // Set the fields in the ISMRMRD Acquisition Header struct mxArray
    mwIndex i;
    for (mwIndex i = 0; i < NUMBER_OF_ACQ_FIELDS; i++) {
        mxSetFieldByNumber(acqhdr, 0, i, acq_values[i]);
    }

    // Copy the data
    std::complex<float> *raw_data = m2->getObjectPtr()->get_data_ptr();
    if (!raw_data) {
        GADGET_DEBUG1("Broken raw_data pointer\n");
        return GADGET_FAIL;
    }

    unsigned long num_elements = m2->getObjectPtr()->get_number_of_elements();

    float *real_data = (float *)mxCalloc(num_elements, sizeof(float));
    if (!real_data) {
        GADGET_DEBUG1("Failed to allocate double* for real_data\n");
        return GADGET_FAIL;
    }
    float *imag_data = (float *)mxCalloc(num_elements, sizeof(float));
    if (!imag_data) {
        GADGET_DEBUG1("Failed to allocate double* for imag_data\n");
        return GADGET_FAIL;
    }

    for (int i = 0; i < num_elements; i++) {
        //std::cout << i << ": " << raw_data[i].real() << ", " << raw_data[i].imag() << endl;
        real_data[i] = raw_data[i].real();
        imag_data[i] = raw_data[i].imag();
    }

    mxArray *acqdata = mxCreateNumericMatrix(0, 0, mxSINGLE_CLASS, mxCOMPLEX);
    mxSetData(acqdata, real_data);
    mxSetImagData(acqdata, imag_data);
    mxSetN(acqdata, m1->getObjectPtr()->number_of_samples);
    mxSetM(acqdata, m1->getObjectPtr()->active_channels);

    engPutVariable(engine_, "acqhdr", acqhdr);
    engPutVariable(engine_, "acqdata", acqdata);

    // Prepare a buffer for collecting Matlab's output
    //char buffer[2049] = "\0";
    //engOutputBuffer(engine_, buffer, 2048);

    // instantiate a Matlab ismrmrd.AcquisitionHeader from our struct
    engEvalString(engine_, "acqhdr = ismrmrd.AcquisitionHeader(acqhdr);");
    engEvalString(engine_, "[res_acqhdr, res_data] = matgadget.process(acqhdr, acqdata);");

    // Convert object back to struct
    engEvalString(engine_, "res_acqhdr = struct(res_acqhdr);");
    engEvalString(engine_, "res_acqhdr.idx = struct(res_acqhdr.idx);");

    mxArray *res_hdr = engGetVariable(engine_, "res_acqhdr");
    if (res_hdr == NULL) {
        GADGET_DEBUG1("Failed to get header back from Matlab\n");
        return GADGET_FAIL;
    }

    // grab the modified AcquisitionHeader and convert it back to C++
    GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m3 =
            new GadgetContainerMessage<ISMRMRD::AcquisitionHeader>();
    ISMRMRD::AcquisitionHeader *hdr_new = m3->getObjectPtr();
    mxArray *tmp;

    tmp = mxGetFieldByNumber(res_hdr, 0, 0);
    memcpy(&hdr_new->version, mxGetData(tmp), mxGetElementSize(tmp));

    tmp = mxGetFieldByNumber(res_hdr, 0, 1);
    memcpy(&hdr_new->flags, mxGetData(tmp), mxGetElementSize(tmp));

    tmp = mxGetFieldByNumber(res_hdr, 0, 2);
    memcpy(&hdr_new->measurement_uid, mxGetData(tmp), mxGetElementSize(tmp));

    tmp = mxGetFieldByNumber(res_hdr, 0, 3);
    memcpy(&hdr_new->scan_counter, mxGetData(tmp), mxGetElementSize(tmp));

    tmp = mxGetFieldByNumber(res_hdr, 0, 4);
    memcpy(&hdr_new->acquisition_time_stamp, mxGetData(tmp), mxGetElementSize(tmp));

    tmp = mxGetFieldByNumber(res_hdr, 0, 5);
    memcpy(&hdr_new->physiology_time_stamp, mxGetData(tmp), ISMRMRD_PHYS_STAMPS * mxGetElementSize(tmp));

    tmp = mxGetFieldByNumber(res_hdr, 0, 6);
    memcpy(&hdr_new->number_of_samples, mxGetData(tmp), mxGetElementSize(tmp));

    tmp = mxGetFieldByNumber(res_hdr, 0, 7);
    memcpy(&hdr_new->available_channels, mxGetData(tmp), mxGetElementSize(tmp));

    tmp = mxGetFieldByNumber(res_hdr, 0, 8);
    memcpy(&hdr_new->active_channels, mxGetData(tmp), mxGetElementSize(tmp));

    tmp = mxGetFieldByNumber(res_hdr, 0, 9);
    memcpy(&hdr_new->channel_mask, mxGetData(tmp), ISMRMRD_CHANNEL_MASKS * mxGetElementSize(tmp));

    tmp = mxGetFieldByNumber(res_hdr, 0, 10);
    memcpy(&hdr_new->discard_pre, mxGetData(tmp), mxGetElementSize(tmp));

    tmp = mxGetFieldByNumber(res_hdr, 0, 11);
    memcpy(&hdr_new->discard_post, mxGetData(tmp), mxGetElementSize(tmp));

    tmp = mxGetFieldByNumber(res_hdr, 0, 12);
    memcpy(&hdr_new->center_sample, mxGetData(tmp), mxGetElementSize(tmp));

    tmp = mxGetFieldByNumber(res_hdr, 0, 13);
    memcpy(&hdr_new->encoding_space_ref, mxGetData(tmp), mxGetElementSize(tmp));

    tmp = mxGetFieldByNumber(res_hdr, 0, 14);
    memcpy(&hdr_new->trajectory_dimensions, mxGetData(tmp), mxGetElementSize(tmp));

    tmp = mxGetFieldByNumber(res_hdr, 0, 15);
    memcpy(&hdr_new->sample_time_us, mxGetData(tmp), ISMRMRD_POSITION_LENGTH * mxGetElementSize(tmp));

    tmp = mxGetFieldByNumber(res_hdr, 0, 16);
    memcpy(&hdr_new->position, mxGetData(tmp), ISMRMRD_DIRECTION_LENGTH * mxGetElementSize(tmp));

    tmp = mxGetFieldByNumber(res_hdr, 0, 17);
    memcpy(&hdr_new->read_dir, mxGetData(tmp), ISMRMRD_DIRECTION_LENGTH * mxGetElementSize(tmp));

    tmp = mxGetFieldByNumber(res_hdr, 0, 18);
    memcpy(&hdr_new->phase_dir, mxGetData(tmp), ISMRMRD_DIRECTION_LENGTH * mxGetElementSize(tmp));

    tmp = mxGetFieldByNumber(res_hdr, 0, 19);
    memcpy(&hdr_new->slice_dir, mxGetData(tmp), ISMRMRD_DIRECTION_LENGTH * mxGetElementSize(tmp));

    tmp = mxGetFieldByNumber(res_hdr, 0, 20);
    memcpy(&hdr_new->patient_table_position, mxGetData(tmp), ISMRMRD_POSITION_LENGTH * mxGetElementSize(tmp));

    // grab and populate the EncodingCounters
    mxArray *res_idx = mxGetFieldByNumber(res_hdr, 0, 21);
    tmp = mxGetFieldByNumber(res_idx, 0, 0);
    memcpy(&hdr_new->idx.kspace_encode_step_1, mxGetData(tmp), mxGetElementSize(tmp));

    tmp = mxGetFieldByNumber(res_idx, 0, 1);
    memcpy(&hdr_new->idx.kspace_encode_step_2, mxGetData(tmp), mxGetElementSize(tmp));

    tmp = mxGetFieldByNumber(res_idx, 0, 2);
    memcpy(&hdr_new->idx.average, mxGetData(tmp), mxGetElementSize(tmp));

    tmp = mxGetFieldByNumber(res_idx, 0, 3);
    memcpy(&hdr_new->idx.slice, mxGetData(tmp), mxGetElementSize(tmp));

    tmp = mxGetFieldByNumber(res_idx, 0, 4);
    memcpy(&hdr_new->idx.contrast, mxGetData(tmp), mxGetElementSize(tmp));

    tmp = mxGetFieldByNumber(res_idx, 0, 5);
    memcpy(&hdr_new->idx.phase, mxGetData(tmp), mxGetElementSize(tmp));

    tmp = mxGetFieldByNumber(res_idx, 0, 6);
    memcpy(&hdr_new->idx.repetition, mxGetData(tmp), mxGetElementSize(tmp));

    tmp = mxGetFieldByNumber(res_idx, 0, 7);
    memcpy(&hdr_new->idx.set, mxGetData(tmp), mxGetElementSize(tmp));

    tmp = mxGetFieldByNumber(res_idx, 0, 8);
    memcpy(&hdr_new->idx.segment, mxGetData(tmp), mxGetElementSize(tmp));

    tmp = mxGetFieldByNumber(res_idx, 0, 9);
    memcpy(&hdr_new->idx.user, mxGetData(tmp), ISMRMRD_USER_INTS * mxGetElementSize(tmp));

    // back to AcquisitionHeader
    tmp = mxGetFieldByNumber(res_hdr, 0, 22);
    memcpy(&hdr_new->user_int, mxGetData(tmp), ISMRMRD_USER_INTS * mxGetElementSize(tmp));

    tmp = mxGetFieldByNumber(res_hdr, 0, 23);
    memcpy(&hdr_new->user_float, mxGetData(tmp), ISMRMRD_USER_FLOATS * mxGetElementSize(tmp));


    mxArray *res_data = engGetVariable(engine_, "res_data");
    if (res_data == NULL) {
        GADGET_DEBUG1("Failed to get data back from Matlab\n");
        return GADGET_FAIL;
    }

    size_t number_of_samples = mxGetN(res_data);
    size_t active_channels = mxGetM(res_data);

    GadgetContainerMessage<hoNDArray< std::complex<float> > >* m4 =
            new GadgetContainerMessage< hoNDArray< std::complex<float> > >();

    m3->cont(m4);
    std::vector<unsigned int> dims;
    dims.push_back(number_of_samples);
    dims.push_back(active_channels);
    if (!m4->getObjectPtr()->create(&dims)) {
        GADGET_DEBUG1("Failed to create new hoNDArray\n");
        return GADGET_FAIL;
    }

    real_data = (float *)mxGetData(res_data);
    imag_data = (float *)mxGetImagData(res_data);
    for (int i = 0; i < number_of_samples*active_channels; i++) {
        m4->getObjectPtr()->get_data_ptr()[i] = std::complex<float>(real_data[i],imag_data[i]);
    }

    //printf("%s", buffer);
    mxDestroyArray(res_hdr);
    mxDestroyArray(res_data);
    mxDestroyArray(acqhdr);
    mxDestroyArray(acqdata);

    return this->next()->putq(m3);
}


int ImageMatlabGadget::process(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1,
        GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2)
{
    ISMRMRD::ImageHeader *img = m1->getObjectPtr();

    mwSize img_dims[2] = {1, 1};
    mxArray *imghdr = mxCreateStructArray(2, img_dims, NUMBER_OF_IMG_FIELDS, ismrmrd_img_field_names);

    mxArray *img_values[NUMBER_OF_IMG_FIELDS];

    img_values[0] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
    memcpy(mxGetData(img_values[0]), &img->version, mxGetElementSize(img_values[0]));

    img_values[1] = mxCreateNumericMatrix(1, 1, mxUINT64_CLASS, mxREAL);
    memcpy(mxGetData(img_values[1]), &img->flags, mxGetElementSize(img_values[1]));

    img_values[2] = mxCreateNumericMatrix(1, 1, mxUINT32_CLASS, mxREAL);
    memcpy(mxGetData(img_values[2]), &img->measurement_uid, mxGetElementSize(img_values[2]));

    img_values[3] = mxCreateNumericMatrix(3, 1, mxUINT16_CLASS, mxREAL);
    memcpy(mxGetData(img_values[3]), &img->matrix_size, 3 * mxGetElementSize(img_values[3]));

    img_values[4] = mxCreateNumericMatrix(3, 1, mxUINT16_CLASS, mxREAL);
    memcpy(mxGetData(img_values[4]), &img->field_of_view, 3 * mxGetElementSize(img_values[4]));

    img_values[5] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
    memcpy(mxGetData(img_values[5]), &img->channels, mxGetElementSize(img_values[5]));

    img_values[6] = mxCreateNumericMatrix(ISMRMRD_POSITION_LENGTH, 1, mxSINGLE_CLASS, mxREAL);
    memcpy(mxGetData(img_values[6]), &img->position, ISMRMRD_POSITION_LENGTH * mxGetElementSize(img_values[6]));

    img_values[7] = mxCreateNumericMatrix(ISMRMRD_DIRECTION_LENGTH, 1, mxSINGLE_CLASS, mxREAL);
    memcpy(mxGetData(img_values[7]), &img->read_dir, ISMRMRD_DIRECTION_LENGTH * mxGetElementSize(img_values[7]));

    img_values[8] = mxCreateNumericMatrix(ISMRMRD_DIRECTION_LENGTH, 1, mxSINGLE_CLASS, mxREAL);
    memcpy(mxGetData(img_values[8]), &img->phase_dir, ISMRMRD_DIRECTION_LENGTH * mxGetElementSize(img_values[8]));

    img_values[9] = mxCreateNumericMatrix(ISMRMRD_DIRECTION_LENGTH, 1, mxSINGLE_CLASS, mxREAL);
    memcpy(mxGetData(img_values[9]), &img->slice_dir, ISMRMRD_DIRECTION_LENGTH * mxGetElementSize(img_values[9]));

    img_values[10] = mxCreateNumericMatrix(ISMRMRD_POSITION_LENGTH, 1, mxSINGLE_CLASS, mxREAL);
    memcpy(mxGetData(img_values[10]), &img->patient_table_position, ISMRMRD_POSITION_LENGTH * mxGetElementSize(img_values[10]));

    img_values[11] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
    memcpy(mxGetData(img_values[11]), &img->average, mxGetElementSize(img_values[11]));

    img_values[12] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
    memcpy(mxGetData(img_values[12]), &img->slice, mxGetElementSize(img_values[12]));

    img_values[13] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
    memcpy(mxGetData(img_values[13]), &img->contrast, mxGetElementSize(img_values[13]));

    img_values[14] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
    memcpy(mxGetData(img_values[14]), &img->phase, mxGetElementSize(img_values[14]));

    img_values[15] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
    memcpy(mxGetData(img_values[15]), &img->repetition, mxGetElementSize(img_values[15]));

    img_values[16] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
    memcpy(mxGetData(img_values[16]), &img->set, mxGetElementSize(img_values[16]));

    img_values[17] = mxCreateNumericMatrix(1, 1, mxUINT32_CLASS, mxREAL);
    memcpy(mxGetData(img_values[17]), &img->acquisition_time_stamp, mxGetElementSize(img_values[17]));

    img_values[18] = mxCreateNumericMatrix(ISMRMRD_PHYS_STAMPS, 1, mxUINT32_CLASS, mxREAL);
    memcpy(mxGetData(img_values[18]), &img->physiology_time_stamp, ISMRMRD_PHYS_STAMPS * mxGetElementSize(img_values[18]));

    img_values[19] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
    memcpy(mxGetData(img_values[19]), &img->image_data_type, mxGetElementSize(img_values[19]));

    img_values[20] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
    memcpy(mxGetData(img_values[20]), &img->image_type, mxGetElementSize(img_values[20]));

    img_values[21] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
    memcpy(mxGetData(img_values[21]), &img->image_index, mxGetElementSize(img_values[21]));

    img_values[22] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
    memcpy(mxGetData(img_values[22]), &img->image_series_index, mxGetElementSize(img_values[22]));

    img_values[23] = mxCreateNumericMatrix(ISMRMRD_USER_INTS, 1, mxINT32_CLASS, mxREAL);
    memcpy(mxGetData(img_values[23]), &img->user_int, ISMRMRD_USER_INTS * mxGetElementSize(img_values[23]));

    img_values[24] = mxCreateNumericMatrix(ISMRMRD_USER_FLOATS, 1, mxSINGLE_CLASS, mxREAL);
    memcpy(mxGetData(img_values[24]), &img->user_float, ISMRMRD_USER_FLOATS * mxGetElementSize(img_values[24]));

    // Set the fields in the ISMRMRD Image Header struct mxArray
    for (mwIndex i = 0; i < NUMBER_OF_IMG_FIELDS; i++) {
        mxSetFieldByNumber(imghdr, 0, i, img_values[i]);
    }

    std::complex<float> *raw_data = m2->getObjectPtr()->get_data_ptr();
    if (!raw_data) {
        GADGET_DEBUG1("Broken raw_data pointer\n");
        return GADGET_FAIL;
    }

    unsigned long num_elements = m2->getObjectPtr()->get_number_of_elements();
    printf("%ld\n", num_elements);
    fflush(stdout);

    double *real_data = (double *)mxCalloc(num_elements, sizeof(double));
    if (!real_data) {
        GADGET_DEBUG1("Failed to allocate double* for real_data\n");
        return GADGET_FAIL;
    }
    double *imag_data = (double *)mxCalloc(num_elements, sizeof(double));
    if (!imag_data) {
        GADGET_DEBUG1("Failed to allocate double* for imag_data\n");
        return GADGET_FAIL;
    }

    for (int i = 0; i < num_elements; i++) {
        //std::cout << i << ": " << raw_data[i].real() << ", " << raw_data[i].imag() << endl;
        real_data[i] = raw_data[i].real();
        imag_data[i] = raw_data[i].imag();
    }

    mxArray *imgdata = mxCreateDoubleMatrix(0, 0, mxCOMPLEX);
    mxSetPr(imgdata, real_data);
    mxSetPi(imgdata, imag_data);
    mxSetN(imgdata, 1);
    mxSetM(imgdata, 2);


    engPutVariable(engine_, "imghdr", imghdr);
    engPutVariable(engine_, "imgdata", imgdata);

    // Prepare a buffer for collecting Matlab's output
    //char buffer[2049] = "\0";
    //engOutputBuffer(engine_, buffer, 2048);

    engEvalString(engine_, "d = imgdata");

    // instantiate a Matlab ismrmrd.AcquisitionHeader from our struct
    engEvalString(engine_, "h = ismrmrd.ImageHeader(imghdr);");

    engEvalString(engine_, "s = struct(h)");

    mxArray *res = engGetVariable(engine_, "s");
    if (res == NULL) {
        GADGET_DEBUG1("Failed to get struct back from ImageHeader in Matlab\n");
        return GADGET_FAIL;
    }

    mxDestroyArray(imghdr);
    mxDestroyArray(imgdata);
    mxDestroyArray(res);

    m1->cont(m2);
    return this->next()->putq(m1);
}


GADGET_FACTORY_DECLARE(AcquisitionMatlabGadget)
GADGET_FACTORY_DECLARE(ImageMatlabGadget)
