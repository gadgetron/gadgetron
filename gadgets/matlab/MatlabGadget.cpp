#include "MatlabGadget.h"

int AcquisitionMatlabGadget::process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1,
        GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2)
{
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
    mxSetN(acq_data, m1->getObjectPtr()->number_of_samples);
    mxSetM(acq_data, m1->getObjectPtr()->active_channels);

    // Prepare a buffer for collecting Matlab's output
    char buffer[2049] = "\0";
    engOutputBuffer(engine_, buffer, 2048);

    engPutVariable(engine_, "acq_hdr_bytes", acq_hdr_bytes);
    engPutVariable(engine_, "acqData", acq_data);
    engEvalString(engine_, "matgadget.putIn(1, acq_hdr_bytes, acqData);");
    engEvalString(engine_, "matgadget.process();");

    printf("%s\n", buffer);

    engEvalString(engine_, "smatgadget = struct(matgadget);");
    mxArray *matgadget = engGetVariable(engine_, "smatgadget");
    if (!matgadget) {
        GADGET_DEBUG1("Failed to get matgadget\n");
        return GADGET_FAIL;
    }

    mxArray *Q = mxGetField(matgadget, 0, "Q");
    if (!Q) {
        GADGET_DEBUG1("Failed to get Q from matgadget\n");
        return GADGET_FAIL;
    }
    if (!mxIsCell(Q)) {
        GADGET_DEBUG1("Q is NOT A CELL");
        return GADGET_FAIL;
    }
    mwSize num_cells = mxGetNumberOfElements(Q);
    mwIndex idx;
    printf("# cells: %ld\n", num_cells);

    for (idx = 0; idx < num_cells; idx++) {
        mxArray *entry = mxGetCell(Q, idx);
        mxArray *res_type = mxGetField(entry, 0, "type");
        mxArray *res_hdr = mxGetField(entry, 0, "head");
        mxArray *res_data = mxGetField(entry, 0, "data");

        int tp = *((int *)mxGetData(res_type));
        printf("idx: %ld, header type: %d\n", idx, tp);
        switch (tp) {
        case 1:     // AcquisitionHeader
        {
            // grab the modified AcquisitionHeader and convert it back to C++
            GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m3 =
                    new GadgetContainerMessage<ISMRMRD::AcquisitionHeader>();
            ISMRMRD::AcquisitionHeader *hdr_new = m3->getObjectPtr();
            memcpy(hdr_new, mxGetData(res_hdr), sizeof(ISMRMRD::AcquisitionHeader));

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
            break;
        }
        default:
            GADGET_DEBUG1("Matlab gadget returned undefined header type\n");
            return GADGET_FAIL;
        }
    }

    // Clear the gadget's queue
    engEvalString(engine_, "matgadget.Q = {};");

    printf("%s", buffer);

    //mxDestroyArray(acq_hdr_bytes);
    //mxDestroyArray(acq_data);


    return GADGET_OK;
}


int ImageMatlabGadget::process(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1,
        GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2)
{
    ISMRMRD::ImageHeader *img = m1->getObjectPtr();

    mwSize img_hdr_dims[2] = {1, sizeof(ISMRMRD::ImageHeader)};
    mxArray *img_hdr_bytes = mxCreateNumericArray(2, img_hdr_dims, mxUINT8_CLASS, mxREAL);
    memcpy(mxGetData(img_hdr_bytes), img, sizeof(ISMRMRD::ImageHeader));

    std::complex<float> *raw_data = m2->getObjectPtr()->get_data_ptr();
    if (!raw_data) {
        GADGET_DEBUG1("Broken raw_data pointer\n");
        return GADGET_FAIL;
    }

    unsigned long num_elements = m2->getObjectPtr()->get_number_of_elements();
    printf("%ld\n", num_elements);
    fflush(stdout);

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

    mxArray *img_data = mxCreateNumericMatrix(0, 0, mxSINGLE_CLASS, mxCOMPLEX);
    mxSetData(img_data, real_data);
    mxSetImagData(img_data, imag_data);
    mxSetN(img_data, 1);
    mxSetM(img_data, 2);


    engPutVariable(engine_, "img_hdr_bytes", img_hdr_bytes);
    engPutVariable(engine_, "img_data", img_data);

    // Prepare a buffer for collecting Matlab's output
    //char buffer[2049] = "\0";
    //engOutputBuffer(engine_, buffer, 2048);

    //engEvalString(engine_, "d = img_data");

    // instantiate a Matlab ismrmrd.AcquisitionHeader from our struct
    //engEvalString(engine_, "h = ismrmrd.ImageHeader(imghdr);");

    //engEvalString(engine_, "s = struct(h)");

    /*
    mxArray *res = engGetVariable(engine_, "s");
    if (res == NULL) {
        GADGET_DEBUG1("Failed to get struct back from ImageHeader in Matlab\n");
        return GADGET_FAIL;
    }
    */

    //mxDestroyArray(res);
    mxDestroyArray(img_hdr_bytes);
    mxDestroyArray(img_data);

    m1->cont(m2);
    return this->next()->putq(m1);
}


GADGET_FACTORY_DECLARE(AcquisitionMatlabGadget)
GADGET_FACTORY_DECLARE(ImageMatlabGadget)
