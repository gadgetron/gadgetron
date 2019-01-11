#include "GrappaCalibrationBuffer.h"
#include "readers/GadgetIsmrmrdReader.h"

namespace Gadgetron {

    GrappaCalibrationBuffer::GrappaCalibrationBuffer(std::vector<size_t> dimensions,
                                                     boost::shared_ptr<GrappaWeights<float> > weights,
                                                     GrappaWeightsCalculator<float> *weights_calculator)
            : weights_(weights), weights_calculator_(weights_calculator), buffer_counter_(dimensions[1]),
              biggest_gap_current_(0), acceleration_factor_(0), last_line_(0), weights_invalid_(true) {
        dimensions_ = dimensions;
        try {
            buffer_.create(&dimensions_);
            buffer_.fill(std::complex<float>(0.0, 0.0));
        } catch (std::runtime_error &err) { GEXCEPTION(err, "Unable to allocate memory for GRAPPA buffer");
        }

    }

    int GrappaCalibrationBuffer::add_data(ISMRMRD::AcquisitionHeader *acq_header, hoNDArray<std::complex<float> > *acq_data,
                                          unsigned short line_offset, unsigned short partition_offset) {
        if (!buffer_.get_data_ptr()) {
            GDEBUG("Buffer not allocated, cannot add data");
            return GADGET_FAIL;
        }

        unsigned int samples = acq_header->number_of_samples;
        unsigned int line = acq_header->idx.kspace_encode_step_1 + line_offset;
        unsigned int partition = acq_header->idx.kspace_encode_step_2 + partition_offset;
        unsigned int slice = acq_header->idx.slice; //We should probably check this // KLK: How? Why?

        if (samples != dimensions_[0]) {
            GDEBUG("Wrong number of samples received\n");
            return GADGET_FAIL;
        }

        std::complex<float> *b = buffer_.get_data_ptr();
        std::complex<float> *d = acq_data->get_data_ptr();

        size_t offset = 0;
        //Copy the data for all the channels
        for (int c = 0; c < acq_header->active_channels; c++) {
            offset =
                    c * dimensions_[0] * dimensions_[1] * dimensions_[2] +
                    partition * dimensions_[0] * dimensions_[1] +
                    line * dimensions_[0];
            memcpy(b + offset, d + c * samples, sizeof(std::complex<float>) * samples);
        }

        int buf_update = buffer_counter_.update_line(line, acq_header->position,
                                                     acq_header->read_dir, acq_header->phase_dir, acq_header->slice_dir);

        if (buf_update < 0) {
            GDEBUG("Unable to update buffer counter for line %d\n", line);
            return GADGET_FAIL;
        }

        //Let's figure out if we should start a weight calculation job

        //This means that the orientation changed
        if (buf_update == 1) {
            weights_invalid_ = true;
        }

        bool is_first_scan_in_slice = acq_header->isFlagSet(ISMRMRD::ISMRMRD_ACQ_FIRST_IN_SLICE);

        //Depending on the sequence used, we could get into trouble if the sequence switches slice acquisition scheme before finishing a slice.
        bool acquiring_sequentially = line > last_line_;

        if (is_first_scan_in_slice) {
            biggest_gap_current_ = 0;
        } else if (acquiring_sequentially) {
            unsigned int gap = std::abs(static_cast<int>(last_line_) - static_cast<int>(line));
            if (gap != biggest_gap_current_) biggest_gap_current_ = gap;
        } else {
            biggest_gap_current_ = 0;
        }
        last_line_ = line;


        bool is_last_scan_in_slice = acq_header->isFlagSet(ISMRMRD::ISMRMRD_ACQ_LAST_IN_SLICE);

        if (is_last_scan_in_slice && acquiring_sequentially) {
            unsigned int min_ky, max_ky;

            if (biggest_gap_current_ != acceleration_factor_) {
                acceleration_factor_ = biggest_gap_current_;
                weights_invalid_ = true;
            }

            if (buffer_counter_.get_region_of_support(min_ky, max_ky) < 0) {
                GDEBUG("Unable to query min_ky, max_ky\n");
                return GADGET_FAIL;
            }

            //If there is nothing on the queue, we might as well recalculate
            if (weights_calculator_->msg_queue()->message_count() < 1) {
                weights_invalid_ = true;
            }

            if (weights_invalid_ && ((max_ky - min_ky) > acceleration_factor_)) {
                std::vector<std::pair<unsigned int, unsigned int> > sampled_region;
                sampled_region.push_back(std::pair<unsigned int, unsigned int>(0, samples - 1));
                sampled_region.push_back(std::pair<unsigned int, unsigned int>(min_ky, max_ky));

                std::vector<unsigned int> uncombined_channel_weights;

                if (!weights_calculator_) {
                    GDEBUG("Weights calculator not defined\n");
                    return GADGET_FAIL;
                }

                weights_calculator_->add_job(&buffer_,
                                             sampled_region,
                                             acceleration_factor_,
                                             weights_,
                                             uncombined_channel_weights,
                                             true);

                weights_invalid_ = false;
            }
        }
        return GADGET_OK;
    }
}
