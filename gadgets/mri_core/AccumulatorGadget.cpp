#include "AccumulatorGadget.h"

namespace {
int addPrePostZeros(size_t centre_column, size_t samples) {
    // 1 : pre zeros
    // 2 : post zeros
    // 0 : no zeros
    if (2 * centre_column == samples) {
        return 0;
    }
    if (2 * centre_column < samples) {
        return 1;
    }
    if (2 * centre_column > samples) {
        return 2;
    }
    return 0;
}
} // namespace

namespace Gadgetron {

AccumulatorGadget::AccumulatorGadget(const Core::Context& context, const Core::GadgetProperties& props)
    : Core::ChannelGadget<Core::Acquisition>(context, props) {
    buffer_ = 0;
    image_counter_ = 0;
    image_series_ = 0;

    auto h = (context.header);
    if (h.encoding.size() != 1) {
        GDEBUG("Number of encoding spaces: %d\n", h.encoding.size());
        GDEBUG("This simple AccumulatorGadget only supports one encoding space\n");
        // TODO: How to throw Gadget failures?
    }

    ISMRMRD::EncodingSpace e_space = h.encoding[0].encodedSpace;
    ISMRMRD::EncodingSpace r_space = h.encoding[0].reconSpace;
    ISMRMRD::EncodingLimits e_limits = h.encoding[0].encodingLimits;

    GDEBUG("Matrix size: %d, %d, %d\n", r_space.matrixSize.x, e_space.matrixSize.y, e_space.matrixSize.z);
    dimensions_.push_back(r_space.matrixSize.x);
    dimensions_.push_back(e_space.matrixSize.y);
    dimensions_.push_back(e_space.matrixSize.z);

    field_of_view_.push_back(r_space.fieldOfView_mm.x);
    field_of_view_.push_back(e_space.fieldOfView_mm.y);
    field_of_view_.push_back(e_space.fieldOfView_mm.z);
    GDEBUG("FOV: %f, %f, %f\n", r_space.fieldOfView_mm.x, e_space.fieldOfView_mm.y, e_space.fieldOfView_mm.z);

    slices_ = e_limits.slice ? e_limits.slice->maximum + 1 : 1;
}

AccumulatorGadget::~AccumulatorGadget() {
    if (buffer_)
        delete buffer_;
}

void AccumulatorGadget::process(Core::InputChannel<Core::Acquisition>& in, Core::OutputChannel& out) {
    for (auto [header, acq, traj] : in) {
        if (!buffer_) {
            dimensions_.push_back(header.active_channels);
            dimensions_.push_back(slices_);

            if (!(buffer_ = new hoNDArray<std::complex<float>>())) {
                GDEBUG("Failed create buffer\n");
              // TODO: How to throw Gadget failures?
            }

            try {
                buffer_->create(dimensions_);
            } 
            catch (std::runtime_error& err) {
                GEXCEPTION(err, "Failed allocate buffer array\n");
                // TODO: How to throw Gadget failures?
            }
            image_series_ = image_series;
        }

        std::complex<float>* b = buffer_->get_data_ptr();
        std::complex<float>* d = acq.get_data_ptr();

        int samples = header.number_of_samples;
        int line = header.idx.kspace_encode_step_1;
        int partition = header.idx.kspace_encode_step_2;
        int slice = header.idx.slice;

        if (samples > static_cast<int>(dimensions_[0])) {
            GDEBUG("Wrong number of samples received\n");
            // TODO: How to throw Gadget failures
        }

        size_t offset = 0;
        // Copy the data for all the channels
        for (int c = 0; c < header.active_channels; c++) {
            offset = slice * dimensions_[0] * dimensions_[1] * dimensions_[2] * dimensions_[3] +
                     c * dimensions_[0] * dimensions_[1] * dimensions_[2] +
                     partition * dimensions_[0] * dimensions_[1] + line * dimensions_[0] + (dimensions_[0] >> 1) -
                     header.center_sample;
            memcpy(b + offset, d + c * samples, sizeof(std::complex<float>) * samples);
        }

        bool is_last_scan_in_slice = header.isFlagSet(ISMRMRD::ISMRMRD_ACQ_LAST_IN_SLICE);

        if (is_last_scan_in_slice) {
            ISMRMRD::ImageHeader newHeader = ISMRMRD::ImageHeader();
            newHeader.clearAllFlags();
            hoNDArray<std::complex<float>> newData = hoNDArray<std::complex<float>>();

            std::vector<size_t> img_dims(4);
            img_dims[0] = dimensions_[0];
            img_dims[1] = dimensions_[1];
            img_dims[2] = dimensions_[2];
            img_dims[3] = dimensions_[3];

            try {
                newData.create(img_dims);
            } catch (std::runtime_error& err) {
                GEXCEPTION(err, "Unable to allocate new image array\n");
                // TODO: How to throw Gadget failures
            }

            size_t data_length = dimensions_[0] * dimensions_[1] * dimensions_[2] * dimensions_[3];

            offset = slice * data_length;

            memcpy(newData.get_data_ptr(), b + offset, sizeof(std::complex<float>) * data_length);

            newHeader.matrix_size[0] = (uint16_t)img_dims[0];
            newHeader.matrix_size[1] = (uint16_t)img_dims[1];
            newHeader.matrix_size[2] = (uint16_t)img_dims[2];
            newHeader.field_of_view[0] = field_of_view_[0];
            newHeader.field_of_view[1] = field_of_view_[1];
            newHeader.field_of_view[2] = field_of_view_[2];
            newHeader.channels = (uint16_t)img_dims[3];
            newHeader.slice = header.idx.slice;

            memcpy(newHeader.position, header.position, sizeof(float) * 3);

            memcpy(newHeader.read_dir, header.read_dir, sizeof(float) * 3);

            memcpy(newHeader.phase_dir, header.phase_dir, sizeof(float) * 3);

            memcpy(newHeader.slice_dir, header.slice_dir, sizeof(float) * 3);

            memcpy(newHeader.patient_table_position, header.patient_table_position, sizeof(float) * 3);

            newHeader.data_type = ISMRMRD::ISMRMRD_CXFLOAT;
            newHeader.image_index = (uint16_t)(++image_counter_);
            newHeader.image_series_index = (uint16_t)image_series_;
            auto newMetaContainer = std::optional<ISMRMRD::MetaContainer>(); // TODO: Should this be empty? Seems like it should be populated somehow
            out.push(Core::Image<std::complex<float>>(newHeader, std::move(newData), newMetaContainer));
        }
    }
}
GADGETRON_GADGET_EXPORT(AccumulatorGadget)
} // namespace Gadgetron
