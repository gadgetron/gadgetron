#include "GrappaWeightsCalculator.h"
#include "GadgetContainerMessage.h"
#include "readers/GadgetIsmrmrdReader.h"
#include "hoNDArray_fileio.h"
#include "hoNDArray_reductions.h"
#include "GadgetronTimer.h"

#ifdef USE_CUDA
#include "GPUTimer.h"
#include "cuNDFFT.h"
#include "b1_map.h"
#include "htgrappa.h"
#include <cuComplex.h>
#endif // USE_CUDA

#include "complext.h"

#include "hoNDFFT.h"
#include "hoNDArray_elemwise.h"
#include "mri_core_grappa.h"
#include "mri_core_coil_map_estimation.h"


namespace Gadgetron {

    int GrappaWeightsCalculator::svc() {

        ACE_Message_Block *mb = nullptr;

        while (this->getq(mb) >= 0) {

            if (mb->msg_type() == ACE_Message_Block::MB_HANGUP) {
                GDEBUG("Hanging up in weights calculator\n");

                if (this->putq(mb) == -1) {
                    GERROR("GrappaWeightsCalculator::svc, putq");
                    return -1;
                }

                break;
            }

            auto *mb1 = AsContainerMessage<WeightsDescription>(mb);
            if (!mb1) {
                mb->release();
                return -2;
            }

            auto *mb2 = AsContainerMessage<hoNDArray<std::complex<float>>>(mb1->cont());

            if (!mb2) {
                mb->release();
                return -3;
            }

            WeightsDescription *description = mb1->getObjectPtr();
            hoNDArray<std::complex<float>> *host_data = mb2->getObjectPtr();

            size_t ks = 5;
            size_t power = 3;

            host_data->squeeze();

            size_t RO = host_data->get_size(0);
            size_t E1 = host_data->get_size(1);
            size_t CHA = host_data->get_size(2);

            std::vector<size_t> data_dimensions;
            host_data->get_dimensions(data_dimensions);

            if (uncombined_channels_.size() > 0) {
                data_dimensions.push_back(uncombined_channels_.size() + 1);
            }

            unmixing_.create(data_dimensions);

            // compute the unmixing coefficients
            size_t numUnCombined = uncombined_channels_.size();

            double thres = 0.0005;
            size_t kRO = 5;
            size_t kNE1 = 4;

            if (numUnCombined == 0) {
                hoNDArray<std::complex<float> > acs(RO, E1, target_coils_, host_data->begin());
                hoNDArray<std::complex<float> > target_acs(RO, E1, target_coils_, acs.begin());

                // estimate coil map
                if (!complex_im_.dimensions_equal(&target_acs)) {
                    complex_im_.create(RO, E1, target_coils_);
                }

                hoNDFFT<float>::instance()->ifft2c(target_acs, complex_im_);
                Gadgetron::coil_map_2d_Inati(complex_im_, coil_map_, ks, power);

                // compute unmixing coefficients
                if (description->acceleration_factor == 1) {
                    Gadgetron::conjugate(coil_map_, coil_map_);
                    Gadgetron::clear(unmixing_);
                    memcpy(unmixing_.begin(), coil_map_.begin(), coil_map_.get_number_of_bytes());
                } else {
                    size_t startRO = description->sampled_region[0].first;
                    size_t endRO = description->sampled_region[0].second;

                    size_t startE1 = description->sampled_region[1].first;
                    size_t endE1 = description->sampled_region[1].second;

                    Gadgetron::grappa2d_calib_convolution_kernel(acs, target_acs, description->acceleration_factor,
                                                                 thres, kRO, kNE1, startRO, endRO, startE1, endE1,
                                                                 conv_ker_);

                    Gadgetron::grappa2d_image_domain_kernel(conv_ker_, RO, E1, kIm_);

                    Gadgetron::clear(unmixing_);

                    Gadgetron::grappa2d_unmixing_coeff(kIm_, coil_map_,
                                                       description->acceleration_factor,
                                                       unmixing_, gFactor_);
                }
            } else {
                hoNDArray<std::complex<float> > acs(RO, E1, CHA, host_data->begin());

                // handle the case that all channels are reconed
                size_t target_coils_with_uncombined = target_coils_ + numUnCombined;
                if (target_coils_with_uncombined > CHA) target_coils_with_uncombined = CHA;

                std::list<unsigned int>::iterator it;

                // compute unmixing coefficients
                if (description->acceleration_factor == 1) {
                    // if no acceleration, the input data is used for coil map estimation
                    // no need to differentiate combined and uncombined channels
                    if (!complex_im_.dimensions_equal(&acs)) {
                        complex_im_.create(RO, E1, acs.get_size(2));
                    }

                    hoNDFFT<float>::instance()->ifft2c(acs, complex_im_);

                    Gadgetron::coil_map_2d_Inati(complex_im_, coil_map_, ks, power);
                    Gadgetron::conjugate(coil_map_, coil_map_);
                    Gadgetron::clear(unmixing_);

                    // copy back to unmixing
                    memcpy(unmixing_.begin(), coil_map_.begin(), sizeof(std::complex<float>) * RO * E1 * CHA);

                    // set uncombined channels
                    size_t ind = 1;
                    for (it = uncombined_channels_.begin(); it != uncombined_channels_.end(); it++, ind++) {
                        std::complex<float> *pUnmixing = unmixing_.begin() + ind * RO * E1 * CHA + (*it) * RO * E1;
                        for (size_t p = 0; p < RO * E1; p++) {
                            pUnmixing[p] = 1;
                        }
                    }
                } else {
                    // first, assemble the target_acs
                    // the combined channel comes first and then all uncombined channels
                    std::vector<size_t> dimTarget(3);
                    dimTarget[0] = RO;
                    dimTarget[1] = E1;
                    dimTarget[2] = target_coils_with_uncombined;

                    if (!target_acs_.dimensions_equal(&dimTarget)) {
                        target_acs_.create(RO, E1, target_coils_with_uncombined);
                    }

                    // copy first combined channels and all uncombined channels to target_acs_
                    size_t sCha, ind(0), ind_uncombined(0);

                    // record from which src channel, a target channel is selected
                    std::vector<size_t> srcChaLoc(target_coils_with_uncombined);
                    for (sCha = 0; sCha < CHA; sCha++) {
                        bool uncombined = false;
                        for (it = uncombined_channels_.begin(); it != uncombined_channels_.end(); it++) {
                            if (sCha == *it) {
                                uncombined = true;
                                break;
                            }
                        }

                        if (!uncombined) {
                            if (ind < (target_coils_with_uncombined - numUnCombined)) {
                                memcpy(
                                        target_acs_.begin() + ind * RO * E1,
                                        acs.begin() + sCha * RO * E1,
                                        sizeof(std::complex<float>) * RO * E1
                                );
                                srcChaLoc[ind] = sCha;
                                ind++;
                            }
                        } else {
                            memcpy(
                                    target_acs_.begin() + (target_coils_with_uncombined - numUnCombined + ind_uncombined) * RO * E1,
                                    acs.begin() + sCha * RO * E1,
                                    sizeof(std::complex<float>) * RO * E1
                            );
                            srcChaLoc[target_coils_with_uncombined - numUnCombined + ind_uncombined] = sCha;
                            ind_uncombined++;
                        }
                    }

                    if (!complex_im_.dimensions_equal(&target_acs_)) {
                        complex_im_.create(RO, E1, target_acs_.get_size(2));
                    }

                    hoNDFFT<float>::instance()->ifft2c(target_acs_, complex_im_);

                    Gadgetron::coil_map_2d_Inati(complex_im_, coil_map_, ks, power);

                    Gadgetron::grappa2d_calib_convolution_kernel(acs, target_acs_,
                                                                 description->acceleration_factor,
                                                                 thres, kRO, kNE1, conv_ker_);

                    Gadgetron::grappa2d_image_domain_kernel(conv_ker_, RO, E1, kIm_);

                    // kIm_ stored the unwrapping coefficients as [RO E1 CHA target_coils_with_uncombined]
                    // for the target_coils_with_uncombined dimension, combined channels come first and then uncombined channels

                    Gadgetron::clear(unmixing_);

                    hoNDArray<std::complex<float> > unmixing_all_channels(RO, E1, CHA, unmixing_.begin());
                    Gadgetron::grappa2d_unmixing_coeff(kIm_, coil_map_,
                                                       description->acceleration_factor,
                                                       unmixing_all_channels, gFactor_);

                    // set unmixing coefficients for uncombined channels
                    for (ind = 0, it = uncombined_channels_.begin(); it != uncombined_channels_.end(); it++, ind++) {
                        memcpy(unmixing_.begin() + ind * RO * E1 * CHA,
                               kIm_.begin() + (target_coils_with_uncombined - numUnCombined + ind) * RO * E1 * CHA,
                               sizeof(std::complex<float>) * RO * E1 * CHA);
                    }
                }
            }

            push_update(description, host_data);

            mb->release();
        }

        return 0;
    }

    void GrappaWeightsCalculator::push_update(const WeightsDescription *description,
                                              const hoNDArray<std::complex<float>> *host_data) const {
        boost::shared_ptr<std::vector<size_t> > tmp_dims = host_data->get_dimensions();

        if (uncombined_channels_.size()) {
            tmp_dims->push_back((size_t) (uncombined_channels_.size() + 1));
        }

        auto unmixing_host = std::make_shared<hoNDArray<std::complex<float>>>(tmp_dims);
        copy(unmixing_.begin(), unmixing_.end(), unmixing_host->begin());

        description->destination->update(unmixing_host.get());
    }

    int GrappaWeightsCalculator::close(unsigned long flags) {
        int rval = 0;
        if (flags == 1) {
            ACE_Message_Block *hangup = new ACE_Message_Block();
            hangup->msg_type(ACE_Message_Block::MB_HANGUP);
            if (this->putq(hangup) == -1) {
                hangup->release();
                GERROR("GrappaWeightsCalculator::close, putq");
                return -1;
            }

            rval = this->wait();
        }
        return rval;
    }


    int GrappaWeightsCalculator::add_job(hoNDArray<std::complex<float>> *ref_data,
                                         std::vector<std::pair<unsigned int, unsigned int>> sampled_region,
                                         unsigned int acceleration_factor,
                                         boost::shared_ptr<GrappaWeights<float>> destination,
                                         std::vector<unsigned int> uncombined_channel_weights,
                                         bool include_uncombined_channels_in_combined_weights) {

        auto *weights_description_message = new GadgetContainerMessage<WeightsDescription>();

        weights_description_message->getObjectPtr()->sampled_region = sampled_region;
        weights_description_message->getObjectPtr()->acceleration_factor = acceleration_factor;
        weights_description_message->getObjectPtr()->destination = destination;
        weights_description_message->getObjectPtr()->uncombined_channel_weights = uncombined_channel_weights;
        weights_description_message->getObjectPtr()->include_uncombined_channels_in_combined_weights =
                include_uncombined_channels_in_combined_weights;

        auto *weights_data_message = new GadgetContainerMessage<hoNDArray<std::complex<float>>>();

        weights_description_message->cont(weights_data_message);

        try { weights_data_message->getObjectPtr()->create(ref_data->get_dimensions().get()); }
        catch (std::runtime_error &err) {
            weights_description_message->release();
            return -3;
        }

        memcpy(weights_data_message->getObjectPtr()->get_data_ptr(), ref_data->get_data_ptr(),
               ref_data->get_number_of_elements() * sizeof(std::complex<float>));

        this->putq(weights_description_message);

        return 0;
    }

    int GrappaWeightsCalculator::add_uncombined_channel(unsigned int channel_id) {
        remove_uncombined_channel(channel_id);
        uncombined_channels_.push_back(channel_id);
        return 0;
    }

    int GrappaWeightsCalculator::remove_uncombined_channel(unsigned int channel_id) {
        uncombined_channels_.remove(channel_id);
        return 0;
    }
}
