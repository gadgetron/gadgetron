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
    : Core::ChannelGadget<mrd::Acquisition>(context, props) {
    buffer_ = 0;
    image_counter_ = 0;
    image_series_ = 0;

    auto h = (context.header);
    if (h.encoding.size() != 1) {
        GDEBUG("Number of encoding spaces: %d\n", h.encoding.size());
        GDEBUG("This simple AccumulatorGadget only supports one encoding space\n");
        // TODO: How to throw Gadget failures?
    }

    mrd::EncodingSpaceType e_space = h.encoding[0].encoded_space;
    mrd::EncodingSpaceType r_space = h.encoding[0].recon_space;
    mrd::EncodingLimitsType e_limits = h.encoding[0].encoding_limits;

    GDEBUG("Matrix size: %d, %d, %d\n", r_space.matrix_size.x, e_space.matrix_size.y, e_space.matrix_size.z);
    dimensions_.push_back(r_space.matrix_size.x);
    dimensions_.push_back(e_space.matrix_size.y);
    dimensions_.push_back(e_space.matrix_size.z);

    field_of_view_[0] = r_space.field_of_view_mm.x;
    field_of_view_[1] = e_space.field_of_view_mm.y;
    field_of_view_[2] = e_space.field_of_view_mm.z;
    GDEBUG("FOV: %f, %f, %f\n", r_space.field_of_view_mm.x, e_space.field_of_view_mm.y, e_space.field_of_view_mm.z);

    slices_ = e_limits.slice ? e_limits.slice->maximum + 1 : 1;
}

AccumulatorGadget::~AccumulatorGadget() {
    if (buffer_)
        delete buffer_;
}

void AccumulatorGadget::process(Core::InputChannel<mrd::Acquisition>& in, Core::OutputChannel& out) {
    mrd::Acquisition ref_acq;
    for (auto acq : in) {
        if (!buffer_) {
            dimensions_.push_back(acq.Coils());
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

        bool is_first_scan_in_slice = acq.head.flags.HasFlags(mrd::AcquisitionFlags::kFirstInSlice);
        if (is_first_scan_in_slice) {
            ref_acq = acq;
        }

        std::complex<float>* buffer_raw = buffer_->get_data_ptr();
        std::complex<float>* data_raw = acq.data.data();

        int samples = acq.Samples();
        int line = acq.head.idx.kspace_encode_step_1.value_or(0);
        int partition = acq.head.idx.kspace_encode_step_2.value_or(0);
        int slice = acq.head.idx.slice.value_or(0);
        int center_sample = acq.head.center_sample.value_or(samples / 2);

        if (samples > dimensions_[0]) {
            GDEBUG("Wrong number of samples received\n");
            // TODO: How to throw Gadget failures
        }

        size_t offset = 0;
        // Copy the data for all the channels
        for (int chan = 0; chan < acq.Coils(); chan++) {
            offset = slice * dimensions_[0] * dimensions_[1] * dimensions_[2] * dimensions_[3] +
                     chan * dimensions_[0] * dimensions_[1] * dimensions_[2] +
                     partition * dimensions_[0] * dimensions_[1] + line * dimensions_[0] + (dimensions_[0] >> 1) -
                     center_sample;
            memcpy(buffer_raw + offset, data_raw + chan * samples, sizeof(std::complex<float>) * samples);
        }

        bool is_last_scan_in_slice = acq.head.flags.HasFlags(mrd::AcquisitionFlags::kLastInSlice);
        if (is_last_scan_in_slice) {
            mrd::Image<std::complex<float>> img;

            std::vector<size_t> img_dims(4);
            img_dims[0] = dimensions_[0];
            img_dims[1] = dimensions_[1];
            img_dims[2] = dimensions_[2];
            img_dims[3] = dimensions_[3];

            img.data.create(img_dims);

            size_t data_length = img.data.size();
            offset = slice * data_length;
            memcpy(img.data.data(), buffer_raw + offset, sizeof(std::complex<float>) * data_length);
            img.head.measurement_uid = acq.head.measurement_uid;
            img.head.field_of_view[0] = field_of_view_[0];
            img.head.field_of_view[1] = field_of_view_[1];
            img.head.field_of_view[2] = field_of_view_[2];
            img.head.position = acq.head.position;
            img.head.col_dir = acq.head.read_dir;
            img.head.line_dir = acq.head.phase_dir;
            img.head.slice_dir = acq.head.slice_dir;
            img.head.patient_table_position = acq.head.patient_table_position;
            img.head.average = acq.head.idx.average;
            img.head.slice = acq.head.idx.slice;
            img.head.contrast = acq.head.idx.contrast;
            img.head.phase = acq.head.idx.phase;
            img.head.repetition = acq.head.idx.repetition;
            img.head.set = acq.head.idx.set;
            img.head.acquisition_time_stamp = acq.head.acquisition_time_stamp;
            img.head.physiology_time_stamp = acq.head.physiology_time_stamp;
            img.head.image_type = mrd::ImageType::kComplex;
            img.head.image_index = ++image_counter_;
            img.head.image_series_index = image_series_;

            out.push(img);
        }
    }
}
GADGETRON_GADGET_EXPORT(AccumulatorGadget)
} // namespace Gadgetron