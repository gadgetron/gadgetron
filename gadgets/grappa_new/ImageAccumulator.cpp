#include "ImageAccumulator.h"

#include <chrono>
#include <ismrmrd/ismrmrd.h>

#include "Node.h"

#include "log.h"

namespace {
    using namespace Gadgetron;
    using namespace Gadgetron::Core;
/*
    void copy_acquisition_data_to_buffer(
            const Acquisition &acquisition,
            hoNDArray<std::complex<float>> &buffer
    ) {
        auto &header = std::get<ISMRMRD::AcquisitionHeader>(acquisition);
        auto &data   = std::get<hoNDArray<std::complex<float>>>(acquisition);

        auto current_slice      = header.idx.slice;
//        auto current_line       = header.idx.kspace_encode_step_1 + line_offset_;
        auto current_partition  = header.idx.kspace_encode_step_2;
        auto samples            = header.number_of_samples;

        std::complex<float> *b = image_data_[slice]->getObjectPtr()->get_data_ptr();
        std::complex<float> *d = m2->getObjectPtr()->get_data_ptr();

        size_t offset = 0;
        //Copy the data for all the channels
        for (int c = 0; c < acq_head->active_channels; c++) {
            offset =
                    c * image_dimensions_[0] * image_dimensions_[1] * image_dimensions_[2] +
                    partition * image_dimensions_[0] * image_dimensions_[1] +
                    line * image_dimensions_[0];

            memcpy(b + offset, d + c * samples, sizeof(std::complex<float>) * samples);
        }
    }
    */
}

namespace Gadgetron::Grappa {

    ImageAccumulator::ImageAccumulator(
            const Core::Context &context,
            const std::unordered_map<std::string, std::string> &props
    ) : TypedGadgetNode<Acquisition>(props) {}

    void ImageAccumulator::process(TypedInputChannel<Acquisition> &in, OutputChannel &out) {

        GINFO_STREAM("Hello, I'm the ImageAccumulator process function. I'm running!");

        std::vector<hoNDArray<std::complex<float>>> image_buffers;
        std::vector<std::chrono::milliseconds> time_stamps;

        for (auto acquisition : in) {
            auto &header = std::get<ISMRMRD::AcquisitionHeader>(acquisition);

            auto slice = header.idx.slice;

            GINFO_STREAM("SLICE! " << slice);

            /*
            copy_acquisition_data_to_buffer(acquisition, image_buffers[slice]);

            if (header.isFlagSet(ISMRMRD::ISMRMRD_ACQ_LAST_IN_SLICE)) {
                push_recon_job(acquisition, image_buffers[slice], out);
            }
             */
        }

        GINFO_STREAM("ImageAccumulator acquisition loop over - channel closed.")
    }

    GADGETRON_GADGET_EXPORT(ImageAccumulator);
}
