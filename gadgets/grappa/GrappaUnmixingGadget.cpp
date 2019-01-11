#include "readers/GadgetIsmrmrdReader.h"
#include "GrappaUnmixingGadget.h"
#include "hoNDFFT.h"

namespace Gadgetron {

    GrappaUnmixingGadget::GrappaUnmixingGadget() {}
    GrappaUnmixingGadget::~GrappaUnmixingGadget() {}

    int GrappaUnmixingGadget::process(GadgetContainerMessage<GrappaUnmixingJob> *unmixing_job_message,
                                      GadgetContainerMessage<ISMRMRD::ImageHeader> *image_header_message,
                                      GadgetContainerMessage<hoNDArray<std::complex<float> > > *image_data_message) {
        auto *image_data = new GadgetContainerMessage<hoNDArray<std::complex<float> > >();

        std::vector<size_t> combined_dims(3, 0);
        combined_dims[0] = image_header_message->getObjectPtr()->matrix_size[0];
        combined_dims[1] = image_header_message->getObjectPtr()->matrix_size[1];
        combined_dims[2] = image_header_message->getObjectPtr()->matrix_size[2];

        if (image_header_message->getObjectPtr()->channels > 1) {
            combined_dims.push_back(image_header_message->getObjectPtr()->channels);
        }

        try { image_data->getObjectPtr()->create(&combined_dims); }
        catch (std::runtime_error &err) { GEXCEPTION(err, "Unable to create combined image array\n");
            return GADGET_FAIL;
        }

        unmixing_job_message->cont(0);
        image_header_message->cont(image_data);

        hoNDFFT<float>::instance()->ifft3c(*image_data_message->getObjectPtr());

        if (!unmixing_job_message->getObjectPtr()->weights_) {
            GDEBUG("Weights are a NULL\n");
            return GADGET_FAIL;
        }

        float scale_factor = 1.0;
        int appl_result = unmixing_job_message->getObjectPtr()->weights_->apply(image_data_message->getObjectPtr(), image_data->getObjectPtr(), scale_factor);
        if (appl_result < 0) {
            GDEBUG("Failed to apply GRAPPA weights: error code %d\n", appl_result);
            return GADGET_FAIL;
        }

        unmixing_job_message->release();
        image_data_message->release();

        if (this->next()->putq(image_header_message) < 0) {
            return GADGET_FAIL;
        }

        return GADGET_OK;
    }

    GADGET_FACTORY_DECLARE(GrappaUnmixingGadget)
}
