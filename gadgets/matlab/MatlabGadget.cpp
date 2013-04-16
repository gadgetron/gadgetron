#include "MatlabGadget.h"

namespace Gadgetron{

int AcquisitionMatlabGadget::process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1,
        GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2)
{
    // Initialize a string for matlab commands
    std::string cmd;

    ISMRMRD::AcquisitionHeader *acq = m1->getObjectPtr();

    mwSize acq_hdr_dims[2] = {1, sizeof(ISMRMRD::AcquisitionHeader)};
    mxArray *acq_hdr_bytes = mxCreateNumericArray(2, acq_hdr_dims, mxUINT8_CLASS, mxREAL);
    memcpy(mxGetData(acq_hdr_bytes), acq, sizeof(ISMRMRD::AcquisitionHeader));

    // Copy the data
    std::complex<float> *raw_data = m2->getObjectPtr()->get_data_ptr();
    if (!raw_data) {
        GADGET_DEBUG1("Broken raw_data pointer\n");
        return GADGET_FAIL;
    }

    unsigned long num_elements = m2->getObjectPtr()->get_number_of_elements();

    float *real_data = (float *)mxCalloc(num_elements, sizeof(float));
    if (!real_data) {
        GADGET_DEBUG1("Failed to allocate float* for real_data\n");
        return GADGET_FAIL;
    }
    float *imag_data = (float *)mxCalloc(num_elements, sizeof(float));
    if (!imag_data) {
        GADGET_DEBUG1("Failed to allocate float* for imag_data\n");
        return GADGET_FAIL;
    }

    for (int i = 0; i < num_elements; i++) {
        //std::cout << i << ": " << raw_data[i].real() << ", " << raw_data[i].imag() << endl;
        real_data[i] = raw_data[i].real();
        imag_data[i] = raw_data[i].imag();
    }

    mxArray *acq_data = mxCreateNumericMatrix(0, 0, mxSINGLE_CLASS, mxCOMPLEX);
    mxSetData(acq_data, real_data);
    mxSetImagData(acq_data, imag_data);
    mxSetM(acq_data, m1->getObjectPtr()->number_of_samples);
    mxSetN(acq_data, m1->getObjectPtr()->active_channels);

    // Logic:
    // send AcquisitionHeader as a byte array
    // send AcquisitionData as a complex float array
    // Call the run_process function in the BaseGadget
    // Empty the gadget's queue.
    // This puts a copy of the queue on the workspace.
    // The queue is a structure array and we read it back
    // TODO put this in a readme file somewhere useful
    engPutVariable(engine_, "hdr_bytes", acq_hdr_bytes);
    engPutVariable(engine_, "data", acq_data);
    cmd = "Q = matgadget.run_process(1, hdr_bytes, data); matgadget.emptyQ();";
    send_matlab_command(cmd);

    // Get the size of the gadget's queue
    mxArray *Q = engGetVariable(engine_, "Q");
    if (Q == NULL) {
        GADGET_DEBUG1("Failed to get the Queue from matgadget\n");
        return GADGET_FAIL;
    }
    size_t qlen = mxGetNumberOfElements(Q);
    //GADGET_DEBUG2("Queue size: %ld", qlen);

    // Loop over the elements of the Q, reading one entry at a time
    // to get a structure with type, headerbytes, and data
    mwIndex idx;
    for (idx = 0; idx < qlen; idx++) {
        mxArray *res_type = mxGetField(Q, idx, "type");
        mxArray *res_hdr  = mxGetField(Q, idx, "bytes");
        mxArray *res_data = mxGetField(Q, idx, "data");

        // determine the type of the object on the quue (i.e. acquisition or image)
        int tp = *((int *)mxGetData(res_type));
        switch (tp) {
        case 1:     // AcquisitionHeader
        {
            // grab the modified AcquisitionHeader and convert it back to C++
            GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m3 =
                    new GadgetContainerMessage<ISMRMRD::AcquisitionHeader>();
            ISMRMRD::AcquisitionHeader *hdr_new = m3->getObjectPtr();
            memcpy(hdr_new, mxGetData(res_hdr), sizeof(ISMRMRD::AcquisitionHeader));

            size_t number_of_samples = mxGetM(res_data);
            size_t active_channels = mxGetN(res_data);

            GadgetContainerMessage<hoNDArray< std::complex<float> > >* m4 =
                    new GadgetContainerMessage< hoNDArray< std::complex<float> > >();

            m3->cont(m4);
            std::vector<unsigned int> dims;
            dims.push_back(number_of_samples);
            dims.push_back(active_channels);
            try {
                m4->getObjectPtr()->create(&dims);
            } catch (bad_alloc& err) {
                GADGET_DEBUG1("Failed to create new hoNDArray\n");
                return GADGET_FAIL;
            }

            float *real_data = (float *)mxGetData(res_data);
            float *imag_data = (float *)mxGetImagData(res_data);
            for (int i = 0; i < number_of_samples*active_channels; i++) {
                m4->getObjectPtr()->get_data_ptr()[i] = std::complex<float>(real_data[i],imag_data[i]);
            }

            if (this->next()->putq(m3) < 0) {
                GADGET_DEBUG1("Failed to put Acquisition message on queue\n");
                return GADGET_FAIL;
            }

            break;
        }
        case 2:     // ImageHeader
        {
            // grab the modified AcquisitionHeader and convert it back to C++
            GadgetContainerMessage<ISMRMRD::ImageHeader>* m3 =
                    new GadgetContainerMessage<ISMRMRD::ImageHeader>();
            ISMRMRD::ImageHeader *hdr_new = m3->getObjectPtr();
            memcpy(hdr_new, mxGetData(res_hdr), sizeof(ISMRMRD::ImageHeader));

            GadgetContainerMessage<hoNDArray< std::complex<float> > >* m4 =
                    new GadgetContainerMessage< hoNDArray< std::complex<float> > >();

            m3->cont(m4);
            std::vector<unsigned int> dims;
            dims.push_back(hdr_new->matrix_size[0]);
            dims.push_back(hdr_new->matrix_size[1]);
            dims.push_back(hdr_new->matrix_size[2]);
            dims.push_back(hdr_new->channels);
            try {
                m4->getObjectPtr()->create(&dims);
            } catch (bad_alloc& err) {
                GADGET_DEBUG1("Failed to create new hoNDArray\n");
                return GADGET_FAIL;
            }

            float *real_data = (float *)mxGetData(res_data);
            float *imag_data = (float *)mxGetImagData(res_data);
            for (int i = 0; i < m4->getObjectPtr()->get_number_of_elements(); i++) {
                m4->getObjectPtr()->get_data_ptr()[i] = std::complex<float>(real_data[i],imag_data[i]);
            }

            if (this->next()->putq(m3) < 0) {
                GADGET_DEBUG1("Failed to put Image message on queue\n");
                return GADGET_FAIL;
            }

            break;
        }
        default:
            GADGET_DEBUG1("Matlab gadget returned undefined header type\n");
            return GADGET_FAIL;
        }
    }

    // Match all mxCreate___s with mxDestroy___s
    mxDestroyArray(acq_hdr_bytes);
    mxDestroyArray(acq_data);

    return GADGET_OK;
}


// TODO: The ImageMatlabGadget is not currently templated
//      It only works for images of type std::complex<float>
int ImageMatlabGadget::process(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1,
        GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2)
{
    // Initialize a string for matlab commands
    std::string cmd;

    ISMRMRD::ImageHeader *img = m1->getObjectPtr();

    // Create a mxArray of bytes for the ISMRMRD::ImageHeader
    mwSize img_hdr_dims[2] = {1, sizeof(ISMRMRD::ImageHeader)};
    mxArray *img_hdr_bytes = mxCreateNumericArray(2, img_hdr_dims, mxUINT8_CLASS, mxREAL);
    memcpy(mxGetData(img_hdr_bytes), img, sizeof(ISMRMRD::ImageHeader));

    // Create a mxArray for the Image data
    std::complex<float> *raw_data = m2->getObjectPtr()->get_data_ptr();
    if (!raw_data) {
        GADGET_DEBUG1("Broken raw_data pointer\n");
        return GADGET_FAIL;
    }

    if (img->matrix_size[0] == 0) img->matrix_size[0] = 1;
    if (img->matrix_size[1] == 0) img->matrix_size[1] = 1;
    if (img->matrix_size[2] == 0) img->matrix_size[2] = 1;
    if (img->channels == 0) img->channels = 1;

    mwSize ndim = 4;
    mwSize dims[4] = {img->matrix_size[0], img->matrix_size[1], img->matrix_size[2], img->channels};
    mxArray *img_data = mxCreateNumericArray(ndim, dims, mxSINGLE_CLASS, mxCOMPLEX);

    float *real_data = (float *)mxGetData(img_data);
    float *imag_data = (float *)mxGetImagData(img_data);
    unsigned long num_elements = m2->getObjectPtr()->get_number_of_elements();
    for (int i = 0; i < num_elements; i++) {
        real_data[i] = raw_data[i].real();
        imag_data[i] = raw_data[i].imag();
    }

    engPutVariable(engine_, "hdr_bytes", img_hdr_bytes);
    engPutVariable(engine_, "data", img_data);
    cmd = "Q = matgadget.run_process(2, hdr_bytes, data); matgadget.emptyQ();";
    send_matlab_command(cmd);

    // Get the size of the gadget's queue
    mxArray *Q = engGetVariable(engine_, "Q");
    if (Q == NULL) {
        GADGET_DEBUG1("Failed to get the Queue from matgadget\n");
        return GADGET_FAIL;
    }
    size_t qlen = mxGetNumberOfElements(Q);

    // Loop over the elements of the Q, reading one entry at a time
    // to get a structure with type, headerbytes, and data
    mwIndex idx;
    for (idx = 0; idx < qlen; idx++) {
        mxArray *res_type = mxGetField(Q, idx, "type");
        mxArray *res_hdr  = mxGetField(Q, idx, "bytes");
        mxArray *res_data = mxGetField(Q, idx, "data");

        // determine the type of the object on the queue (i.e. acquisition or image)
        // although, since this is an Image gadget, it better be an image
        int tp = *((int *)mxGetData(res_type));
        switch (tp) {
        case 2:     // ImageHeader
        {
            // grab the modified AcquisitionHeader and convert it back to C++
            GadgetContainerMessage<ISMRMRD::ImageHeader>* m3 =
                    new GadgetContainerMessage<ISMRMRD::ImageHeader>();
            ISMRMRD::ImageHeader *hdr_new = m3->getObjectPtr();
            memcpy(hdr_new, mxGetData(res_hdr), sizeof(ISMRMRD::ImageHeader));

            GadgetContainerMessage<hoNDArray< std::complex<float> > >* m4 =
                    new GadgetContainerMessage< hoNDArray< std::complex<float> > >();

            m3->cont(m4);
            std::vector<unsigned int> dims;
            dims.push_back(hdr_new->matrix_size[0]);
            dims.push_back(hdr_new->matrix_size[1]);
            dims.push_back(hdr_new->matrix_size[2]);
            dims.push_back(hdr_new->channels);
            try {
                m4->getObjectPtr()->create(&dims);
            } catch (bad_alloc& err) {
                GADGET_DEBUG1("Failed to create new hoNDArray\n");
                return GADGET_FAIL;
            }

            float *real_data = (float *)mxGetData(res_data);
            float *imag_data = (float *)mxGetImagData(res_data);
            for (int i = 0; i < m4->getObjectPtr()->get_number_of_elements(); i++) {
                m4->getObjectPtr()->get_data_ptr()[i] = std::complex<float>(real_data[i],imag_data[i]);
            }

            if (this->next()->putq(m3) < 0) {
                GADGET_DEBUG1("Failed to put Image message on queue\n");
                return GADGET_FAIL;
            }

            break;
        }
        default:
            GADGET_DEBUG1("Matlab gadget returned undefined header type\n");
            return GADGET_FAIL;
        }
    }

    // Match all mxCreate___s with mxDestroy___s
    mxDestroyArray(img_hdr_bytes);
    mxDestroyArray(img_data);

    return GADGET_OK;
}


GADGET_FACTORY_DECLARE(AcquisitionMatlabGadget)
GADGET_FACTORY_DECLARE(ImageMatlabGadget)
}
