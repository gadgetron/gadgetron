#include "GrappaGadget.h"
#include "GrappaUnmixingGadget.h"
#include "ismrmrd/xml.h"

#include <ace/OS_NS_stdlib.h>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>


namespace Gadgetron {

    GrappaGadget::GrappaGadget()
            : image_counter_(0), image_series_(0), first_call_(true), target_coils_(0), use_gpu_(true) {
    }

    GrappaGadget::~GrappaGadget() {
        for (unsigned int i = 0; i < buffers_.size(); i++) {
            if (buffers_[i]) delete buffers_[i];
            buffers_[i] = 0;


            if (image_data_[i]) {
                image_data_[i]->release();
                image_data_[i] = 0;
            }
        }
    }

    int GrappaGadget::close(unsigned long flags) {
        int ret = Gadget::close(flags);
        GDEBUG("Shutting down GRAPPA Gadget\n");

        if (weights_calculator_.close(flags) < 0) {
            GDEBUG("Failed to close down weights calculator\n");
            return GADGET_FAIL;
        }

        return ret;
    }

    int GrappaGadget::process_config(ACE_Message_Block *mb) {
        ISMRMRD::IsmrmrdHeader h;
        ISMRMRD::deserialize(mb->rd_ptr(), h);

        if (h.encoding.size() != 1) {
            GDEBUG("Number of encoding spaces: %d\n", h.encoding.size());
            GDEBUG("This Gadget only supports one encoding space\n");
            return GADGET_FAIL;
        }

        ISMRMRD::EncodingSpace e_space = h.encoding[0].encodedSpace;
        ISMRMRD::EncodingSpace r_space = h.encoding[0].reconSpace;
        ISMRMRD::EncodingLimits e_limits = h.encoding[0].encodingLimits;

        unsigned int slices = e_limits.slice ? e_limits.slice->maximum + 1 : 1;
        dimensions_.push_back(e_space.matrixSize.x);
        dimensions_.push_back(e_space.matrixSize.y);
        dimensions_.push_back(e_space.matrixSize.z);
        dimensions_.push_back((h.acquisitionSystemInformation && h.acquisitionSystemInformation->receiverChannels)
                              ? *(h.acquisitionSystemInformation->receiverChannels) : 1);
        dimensions_.push_back(slices);

        GDEBUG("Dimensions %d, %d, %d, %d, %d\n", dimensions_[0], dimensions_[1], dimensions_[2], dimensions_[3],
               dimensions_[4]);

        image_dimensions_.push_back(r_space.matrixSize.x);
        image_dimensions_.push_back(r_space.matrixSize.y);
        image_dimensions_.push_back(r_space.matrixSize.z);
        image_dimensions_.push_back(dimensions_[3]);

        fov_.push_back(r_space.fieldOfView_mm.x);
        fov_.push_back(r_space.fieldOfView_mm.y);
        fov_.push_back(r_space.fieldOfView_mm.z);

        line_offset_ = (dimensions_[1] / 2) - e_limits.kspace_encoding_step_1->center;

        if (h.userParameters) {
            for (size_t i = 0; i < h.userParameters->userParameterString.size(); i++) {
                std::string name = h.userParameters->userParameterString[i].name;
                std::string value = h.userParameters->userParameterString[i].value;
                if (name.substr(0, 5) == std::string("COIL_")) {
                    int coil_num = std::atoi(name.substr(5, name.size() - 5).c_str());
                    channel_map_[value] = coil_num;
                }
            }
        }

        return GADGET_OK;
    }


    int GrappaGadget::initial_setup() {


        weights_ = std::vector<boost::shared_ptr<GrappaWeights<float>>>(dimensions_[4]);

        buffers_ = std::vector<GrappaCalibrationBuffer *>(dimensions_[4], 0);
        time_stamps_ = std::vector<ACE_UINT32>(dimensions_[4], 0);

        //Let's figure out the number of target coils
        target_coils_ = target_coils.value();
        if ((target_coils_ <= 0) || (target_coils_ > dimensions_[3])) {
            target_coils_ = dimensions_[3];
        }

        GDEBUG("Running GRAPPA recon with %d source channels and %d target channels\n", dimensions_[3], target_coils_);

        weights_calculator_.set_number_of_target_coils(target_coils_);

        bool use_gpu_ = use_gpu.value();
        GDEBUG_STREAM("use_gpu_ is " << use_gpu_);

        weights_calculator_.set_use_gpu(use_gpu_);

        if (device_channels.value()) {
            GDEBUG("We got the number of device channels from other gadget: %d\n", device_channels.value());
            for (int i = 0; i < device_channels.value(); i++) {
                weights_calculator_.add_uncombined_channel((unsigned int) i);
            }
        } else {
            //Let's figure out if we have channels that are supposed to be uncombined
            std::string uncomb_str = uncombined_channels.value();
            std::vector<std::string> uncomb;
            boost::split(uncomb, uncomb_str, boost::is_any_of(","));
            for (unsigned int i = 0; i < uncomb.size(); i++) {
                std::string ch = boost::algorithm::trim_copy(uncomb[i]);
                if (ch.size() > 0) {
                    unsigned int channel_id = static_cast<unsigned int>(ACE_OS::atoi(ch.c_str()));
                    weights_calculator_.add_uncombined_channel(channel_id);
                }
            }

            uncomb_str = uncombined_channels_by_name.value();
            if (uncomb_str.size()) {
                GDEBUG("uncomb_str: %s\n", uncomb_str.c_str());
                boost::split(uncomb, uncomb_str, boost::is_any_of(","));
                for (unsigned int i = 0; i < uncomb.size(); i++) {
                    std::string ch = boost::algorithm::trim_copy(uncomb[i]);
                    map_type_::iterator it = channel_map_.find(ch);
                    if (it != channel_map_.end()) {
                        unsigned int channel_id = static_cast<unsigned int>(it->second);
                        GDEBUG("Device channel: %s (%d)\n", uncomb[i].c_str(), channel_id);
                        weights_calculator_.add_uncombined_channel(channel_id);
                    }
                }
            }
        }


        for (unsigned int i = 0; i < buffers_.size(); i++) {
            weights_[i] = boost::shared_ptr<GrappaWeights<float> >(new GrappaWeights<float>());
            buffers_[i] = new GrappaCalibrationBuffer(image_dimensions_,
                                                      weights_[i],
                                                      &weights_calculator_);
        }


        if (weights_calculator_.open() < 0) {
            GDEBUG("Failed to open GrappaWeightsCalculator\n");
            return GADGET_FAIL;
        }

        image_data_ = std::vector<GadgetContainerMessage<hoNDArray<std::complex<float> > > *>(dimensions_[4], 0);
        for (unsigned int i = 0; i < image_data_.size(); i++) {
            if (create_image_buffer(i) != GADGET_OK) {
                GDEBUG("Unable to create image buffers");
                return GADGET_FAIL;
            }
        }

        image_series_ = image_series.value();

        return GADGET_OK;
    }


    int GrappaGadget::
    process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader> *m1,
            GadgetContainerMessage<hoNDArray<std::complex<float> > > *m2) {
        bool is_noise = m1->getObjectPtr()->isFlagSet(ISMRMRD::ISMRMRD_ACQ_IS_NOISE_MEASUREMENT);

        //We should not be receiving noise here
        if (is_noise) {
            m1->release();
            return GADGET_OK;
        }

        if (first_call_) {
            if (m1->getObjectPtr()->active_channels != dimensions_[3]) {
                GDEBUG("Detected coil number change. Maybe due to upstream channel reduction\n");
                dimensions_[3] = m1->getObjectPtr()->active_channels;
            }

            if (initial_setup() != GADGET_OK) {
                GDEBUG("Initial Setup Failed\n");
                m1->release();
                return GADGET_FAIL;
            }
            first_call_ = false;
        }

        ISMRMRD::AcquisitionHeader *acq_head = m1->getObjectPtr();

        unsigned int samples = acq_head->number_of_samples;
        unsigned int line = acq_head->idx.kspace_encode_step_1 + line_offset_;
        unsigned int partition = acq_head->idx.kspace_encode_step_2;
        unsigned int slice = acq_head->idx.slice;

        if (samples != image_dimensions_[0]) {
            GDEBUG("GrappaGadget: wrong number of samples received %d, expected %d\n", samples, image_dimensions_[0]);
            return GADGET_FAIL;
        }

        if (slice >= image_data_.size()) {
            GDEBUG("Invalid slice number received\n");
            return GADGET_FAIL;
        }

        if (!image_data_[0]) {
            if (create_image_buffer(slice) != GADGET_OK) {
                GDEBUG("Failed to allocate new slice buffer\n");
                return GADGET_FAIL;
            }
        }

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


        bool is_last_scan_in_slice = acq_head->isFlagSet(ISMRMRD::ISMRMRD_ACQ_LAST_IN_SLICE);
        bool is_first_scan_in_slice = acq_head->isFlagSet(ISMRMRD::ISMRMRD_ACQ_FIRST_IN_SLICE);

        if (is_first_scan_in_slice) {
            time_stamps_[slice] = acq_head->acquisition_time_stamp;
        }

        if (is_last_scan_in_slice) {

            GadgetContainerMessage<GrappaUnmixingJob> *unmixing_job_container =
                    new GadgetContainerMessage<GrappaUnmixingJob>();

            GadgetContainerMessage<ISMRMRD::ImageHeader> *image_header_container =
                    new GadgetContainerMessage<ISMRMRD::ImageHeader>();

            image_header_container->getObjectPtr()->matrix_size[0] = image_dimensions_[0];
            image_header_container->getObjectPtr()->matrix_size[1] = image_dimensions_[1];
            image_header_container->getObjectPtr()->matrix_size[2] = image_dimensions_[2];

            image_header_container->getObjectPtr()->field_of_view[0] = fov_[0];
            image_header_container->getObjectPtr()->field_of_view[1] = fov_[1];
            image_header_container->getObjectPtr()->field_of_view[2] = fov_[2];

            image_header_container->getObjectPtr()->channels = 1 + weights_calculator_.get_number_of_uncombined_channels();
            image_header_container->getObjectPtr()->slice = m1->getObjectPtr()->idx.slice;
            image_header_container->getObjectPtr()->acquisition_time_stamp = time_stamps_[slice];

            memcpy(image_header_container->getObjectPtr()->position, m1->getObjectPtr()->position,
                   sizeof(float) * 3);

            memcpy(image_header_container->getObjectPtr()->read_dir, m1->getObjectPtr()->read_dir,
                   sizeof(float) * 3);

            memcpy(image_header_container->getObjectPtr()->phase_dir, m1->getObjectPtr()->phase_dir,
                   sizeof(float) * 3);

            memcpy(image_header_container->getObjectPtr()->slice_dir, m1->getObjectPtr()->slice_dir,
                   sizeof(float) * 3);

            memcpy(image_header_container->getObjectPtr()->patient_table_position, m1->getObjectPtr()->patient_table_position,
                   sizeof(float) * 3);

            image_header_container->getObjectPtr()->image_index = ++image_counter_;
            image_header_container->getObjectPtr()->image_series_index = image_series_;

            unmixing_job_container->getObjectPtr()->weights_ = weights_[slice];
            unmixing_job_container->cont(image_header_container);
            image_header_container->cont(image_data_[slice]);

            image_data_[slice] = 0;
            if (create_image_buffer(slice) != GADGET_OK) {
                GDEBUG("Failed to create image buffer");
                return GADGET_FAIL;
            }

            if (this->next()->putq(unmixing_job_container) < 0) {
                GDEBUG("Failed to pass image on to next Gadget in chain\n");
                return GADGET_FAIL;
            }
        }

        if (buffers_[slice]->add_data(m1->getObjectPtr(), m2->getObjectPtr(), line_offset_) < 0) {
            GDEBUG("Failed to add incoming data to grappa calibration buffer\n");
            return GADGET_FAIL;
        }

        m1->release();
        return GADGET_OK;
    }


    int GrappaGadget::create_image_buffer(unsigned int slice) {
        if (slice >= image_data_.size()) {
            return GADGET_FAIL;
        }

        if (image_data_[slice] != 0) {
            image_data_[slice]->release();
            image_data_[slice] = 0;
        }

        image_data_[slice] = new GadgetContainerMessage<hoNDArray<std::complex<float> > >();
        try { image_data_[slice]->getObjectPtr()->create(&image_dimensions_); }
        catch (std::runtime_error &err) { GEXCEPTION(err, "Unable to create image buffers");
            return GADGET_FAIL;
        }

        std::fill(image_data_[slice]->getObjectPtr()->get_data_ptr(),
                  image_data_[slice]->getObjectPtr()->get_data_ptr() +
                  image_data_[slice]->getObjectPtr()->get_number_of_elements(),
                  std::complex<float>(0.0f, 0.0f));

        return GADGET_OK;

    }

    GADGET_FACTORY_DECLARE(GrappaGadget)
}
