#include "MatlabGadget.h"

#define NUMBER_OF_FIELDS    6 //23

const char *ismrmrd_acq_field_names[] = {
    "version",
    "flags",
    "measurement_uid",
    "scan_counter",
    "acquisition_time_stamp",
    "physiology_time_stamp",
    /*
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
 //   "EncodingCounters   idx",
    "user_int",
    "user_float"
    */
};

int AcquisitionMatlabGadget::process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1,
        GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2)
{
    if (this->next()->putq(m1) < 0) {
        return GADGET_FAIL;
    }

    ISMRMRD::AcquisitionHeader *acq = m1->getObjectPtr();
    hoNDArray< std::complex<float> > *data = m2->getObjectPtr();

    std::cout << "Preparing Matlab struct array for Acquisition header" << std::endl;

    // Create struct array for storing a single ISMRMRD Acquisition Header
    mwSize dims[2] = {1, 1};
    mxArray *acqhdr = mxCreateStructArray(2, dims, NUMBER_OF_FIELDS, ismrmrd_acq_field_names);

    std::cout << "Created struct array" << std::endl;

    mxArray *acqhdr_values[NUMBER_OF_FIELDS];

    acqhdr_values[0] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
    memcpy(mxGetData(acqhdr_values[0]), &acq->version, mxGetElementSize(acqhdr_values[0]));

    acqhdr_values[1] = mxCreateNumericMatrix(1, 1, mxUINT64_CLASS, mxREAL);
    memcpy(mxGetData(acqhdr_values[1]), &acq->flags, mxGetElementSize(acqhdr_values[1]));

    acqhdr_values[2] = mxCreateNumericMatrix(1, 1, mxUINT32_CLASS, mxREAL);
    memcpy(mxGetData(acqhdr_values[2]), &acq->measurement_uid, mxGetElementSize(acqhdr_values[2]));

    acqhdr_values[3] = mxCreateNumericMatrix(1, 1, mxUINT32_CLASS, mxREAL);
    memcpy(mxGetData(acqhdr_values[3]), &acq->scan_counter, mxGetElementSize(acqhdr_values[3]));

    acqhdr_values[4] = mxCreateNumericMatrix(1, 1, mxUINT32_CLASS, mxREAL);
    memcpy(mxGetData(acqhdr_values[4]), &acq->acquisition_time_stamp, mxGetElementSize(acqhdr_values[4]));

    acqhdr_values[5] = mxCreateNumericMatrix(3, 1, mxUINT32_CLASS, mxREAL);
    memcpy(mxGetData(acqhdr_values[5]), acq->physiology_time_stamp, 3 * mxGetElementSize(acqhdr_values[5]));

    std::cout << "Set 6 array values" << std::endl;
    fflush(stdout);
/*
    acqhdr_values[6] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
    mxSetData(acqhdr_values[6], &macq->number_of_samples);

    acqhdr_values[7] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
    mxSetData(acqhdr_values[7], &macq->available_channels);

    acqhdr_values[8] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
    mxSetData(acqhdr_values[8], &macq->active_channels);

    acqhdr_values[9] = mxCreateNumericMatrix(ISMRMRD_CHANNEL_MASKS, 1, mxUINT64_CLASS, mxREAL);
    mxSetData(acqhdr_values[9], macq->channel_mask);

    acqhdr_values[10] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
    mxSetData(acqhdr_values[10], &macq->discard_pre);

    acqhdr_values[11] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
    mxSetData(acqhdr_values[11], &macq->discard_post);

    acqhdr_values[12] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
    mxSetData(acqhdr_values[12], &macq->center_sample);

    acqhdr_values[13] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
    mxSetData(acqhdr_values[13], &macq->encoding_space_ref);

    acqhdr_values[14] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
    mxSetData(acqhdr_values[14], &macq->trajectory_dimensions);

    acqhdr_values[15] = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
    mxSetData(acqhdr_values[15], &macq->sample_time_us);

    acqhdr_values[16] = mxCreateNumericMatrix(ISMRMRD_POSITION_LENGTH, 1, mxSINGLE_CLASS, mxREAL);
    mxSetData(acqhdr_values[16], macq->position);

    acqhdr_values[17] = mxCreateNumericMatrix(ISMRMRD_DIRECTION_LENGTH, 1, mxSINGLE_CLASS, mxREAL);
    mxSetData(acqhdr_values[17], macq->read_dir);

    acqhdr_values[18] = mxCreateNumericMatrix(ISMRMRD_DIRECTION_LENGTH, 1, mxSINGLE_CLASS, mxREAL);
    mxSetData(acqhdr_values[18], macq->phase_dir);

    acqhdr_values[19] = mxCreateNumericMatrix(ISMRMRD_DIRECTION_LENGTH, 1, mxSINGLE_CLASS, mxREAL);
    mxSetData(acqhdr_values[19], macq->slice_dir);

    acqhdr_values[20] = mxCreateNumericMatrix(ISMRMRD_POSITION_LENGTH, 1, mxSINGLE_CLASS, mxREAL);
    mxSetData(acqhdr_values[20], macq->patient_table_position);

    //acqhdr_values[21] = mxCreateNumericMatrix(1, 1, mxUINT32_CLASS, mxREAL);
    acqhdr_values[21] = mxCreateNumericMatrix(ISMRMRD_USER_INTS, 1, mxINT32_CLASS, mxREAL);
    mxSetData(acqhdr_values[21], macq->user_int);

    acqhdr_values[22] = mxCreateNumericMatrix(ISMRMRD_USER_FLOATS, 1, mxSINGLE_CLASS, mxREAL);
    mxSetData(acqhdr_values[22], macq->user_float);
*/
    // Set the fields in the ISMRMRD Acquisition Header struct mxArray
    mwIndex i;
    for (i = 0; i < NUMBER_OF_FIELDS; i++) {
        mxSetFieldByNumber(acqhdr, 0, i, acqhdr_values[i]);
    }

    engPutVariable(engine_, "acqhdr", acqhdr);

    // Prepare a buffer for collecting Matlab's output
    char buffer[2049] = "\0";
    engOutputBuffer(engine_, buffer, 2048);

    // instantiate a Matlab ismrmrd.AcquisitionHeader from our struct
    engEvalString(engine_, "h = ismrmrd.AcquisitionHeader(acqhdr);");

    engEvalString(engine_, "h.scan_counter = 42;");

    // cast EncodingCounters to a struct
    engEvalString(engine_, "idx = struct(h.idx);");

    engEvalString(engine_, "s = struct(h)");

    // overwrite the EncodingCounter class with the EncodingCounters struct
    engEvalString(engine_, "s.idx = idx;");

    mxArray *res = engGetVariable(engine_, "s");
    if (res == NULL) {
        GADGET_DEBUG1("Failed to get struct back from AcquisitionHeader in Matlab\n");
        return GADGET_FAIL;
    }

    printf("%s", buffer);
    printf("s is class %s\n", mxGetClassName(res));

    mxArray *scan_counter_array = mxGetField(res, 0, "scan_counter");
    uint32_t scan_counter = ((uint32_t *)mxGetData(scan_counter_array))[0];
    printf("Scan counter: %u\n", scan_counter);
    mxDestroyArray(scan_counter_array);

    mxDestroyArray(acqhdr);
    mxDestroyArray(res);

    return GADGET_OK;
}


int ImageMatlabGadget::process(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1,
        GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2)
{

    if (this->next()->putq(m1) < 0) {
        return GADGET_FAIL;
    }

    ISMRMRD::ImageHeader *img = m1->getObjectPtr();
    hoNDArray< std::complex<float> > *data = m2->getObjectPtr();

    return GADGET_OK;
}


GADGET_FACTORY_DECLARE(AcquisitionMatlabGadget)
GADGET_FACTORY_DECLARE(ImageMatlabGadget)
