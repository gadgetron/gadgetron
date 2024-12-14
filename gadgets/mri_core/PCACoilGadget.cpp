/*
* PCACoilGadget.cpp
*
*  Created on: Dec 13, 2011
*      Author: Michael S. Hansen
*/

#include "PCACoilGadget.h"
#include "hoNDArray_elemwise.h"
#include "hoNDArray_fileio.h"
#include "hoNDKLT.h"
#include "hoNDArray_linalg.h"

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>

namespace Gadgetron {

    namespace {
        int get_location(const mrd::Acquisition& acq)
        {
            return acq.head.idx.slice.value_or(0);
        }
    }

    PCACoilGadget::PCACoilGadget(const Core::Context& context, const Core::GadgetProperties& props)
        : Core::ChannelGadget<mrd::Acquisition>(context, props)
    {
        std::vector<std::string> uncomb;
        if (uncombined_channels_by_name.size()) {
            GDEBUG("uncomb_str: %s\n", uncombined_channels_by_name.c_str());
            boost::split(uncomb, uncombined_channels_by_name, boost::is_any_of(","));
            for (unsigned int i = 0; i < uncomb.size(); i++) {
                std::string ch = boost::algorithm::trim_copy(uncomb[i]);
                if (context.header.acquisition_system_information) {
                    for (size_t i = 0; i < context.header.acquisition_system_information->coil_label.size(); i++) {
                        if (ch == context.header.acquisition_system_information->coil_label[i].coil_name) {
                            uncombined_channels_.push_back(i);//This assumes that the channels are sorted in the header
                            break;
                        }
                    }
                }
            }
        }

        /** NOTE:
         *
         * This PCACoilGadget used to have a GADGET_PROPERTY "present_uncombined_channels" that was
         * updated here to `uncombined_channels_.size()`, then later referenced by the NoiseAdjustGadget
         * **only** in the `interventional_mri/grappa_device.xml` chain, which is untested.
         *
         * However, ChannelGadgets don't seem to support a non-const GADGET_PROPERTY. They only support
         * NODE_PROPERTY, which declares a const member variable.
         *
         * Since the grappa_device chains is not tested anywhere, I think this is a "dead" feature.
         *
         */
        // present_uncombined_channels.value((int)uncombined_channels_.size());
        // GDEBUG("Number of uncombined channels (present_uncombined_channels) set to %d\n", uncombined_channels_.size());

#ifdef USE_OMP
        omp_set_num_threads(1);
#endif // USE_OMP

    }

    PCACoilGadget::~PCACoilGadget()
    {
        std::map<int, hoNDKLT<std::complex<float> >* >::iterator it;
        it = pca_coefficients_.begin();
        while (it != pca_coefficients_.end()) {
            if (it->second) {
                delete it->second;
                it->second = 0;
            }
            it++;
        }
    }


    void PCACoilGadget::process(Core::InputChannel<mrd::Acquisition>& input, Core::OutputChannel& output)
    {
        for (auto acq : input) {
            bool is_noise = acq.head.flags.HasFlags(mrd::AcquisitionFlags::kIsNoiseMeasurement);

            //We should not be receiving noise here
            if (is_noise) {
                GERROR("Received noise in PCACoilGadget\n");
                continue;
            }

            int location = get_location(acq);

            // We will always start in buffering mode for a given location.
            buffering_mode_.try_emplace(location, true);
            bool is_buffering = buffering_mode_[location];

            if (is_buffering)
            {
                buffer_[location].push_back(acq);

                bool is_last_scan_in_slice = acq.head.flags.HasFlags(mrd::AcquisitionFlags::kLastInSlice);
                size_t profiles_available = buffer_[location].size();

                //Are we ready for calculating PCA
                if (is_last_scan_in_slice || (profiles_available >= max_buffered_profiles_))
                {
                    calculate_coefficients(location);

                    //Now we should pump all the profiles that we have buffered back through the system
                    for (size_t p = 0; p < profiles_available; p++) {
                        auto& acq = buffer_[location][p];
                        do_pca(acq);
                        output.push(std::move(acq));
                    }

                    //Switch off buffering for this slice
                    buffering_mode_[location] = false;
                    //Remove references in this buffer
                    buffer_[location].clear();
                }
            }
            else {
                // GDEBUG_STREAM("Not buffering location " << location << " anymore");
                do_pca(acq);
                output.push(std::move(acq));
            }
        }
    }

    void PCACoilGadget::calculate_coefficients(int location)
    {
        size_t profiles_available = buffer_[location].size();

        mrd::Acquisition& ref = buffer_[location][0];
        size_t samples_per_profile = ref.Samples();
        size_t channels = ref.Coils();

        GDEBUG("Calculating PCA coefficients with %d profiles for %d coils\n", profiles_available, channels);
        size_t samples_to_use = samples_per_profile > samples_to_use_ ? samples_to_use_ : samples_per_profile;

        //For some sequences there is so little data, we should just use it all.
        if (profiles_available < 16) {
            samples_to_use = samples_per_profile;
        }

        size_t total_samples = samples_to_use*profiles_available;

        std::vector<size_t> dims{total_samples, channels};
        hoNDArray<std::complex<float>> A(dims);

        size_t sample_counter = 0;

        size_t data_offset = 0;
        auto center_sample = ref.head.center_sample.value_or(0);
        if (center_sample >= (samples_to_use >> 1)) {
            data_offset = center_sample - (samples_to_use >> 1);
        }

        //GDEBUG("Data offset = %d\n", data_offset);

        std::vector<size_t> means_dims{channels};
        hoNDArray<std::complex<float> > means(means_dims);
        means.fill(std::complex<float>(0.0f, 0.0f));

        GDEBUG_STREAM("Finished setting up arrays");

        for (size_t p = 0; p < profiles_available; p++) {
            mrd::Acquisition& tmp = buffer_[location][p];

            for (size_t s = 0; s < samples_to_use; s++) {
                for (size_t c = 0; c < channels; c++) {
                    bool uncombined_channel = std::find(uncombined_channels_.begin(), uncombined_channels_.end(), c) != uncombined_channels_.end();
                    if (uncombined_channel) {
                        A(sample_counter, c) = std::complex<float>(0.0, 0.0);
                    } else {
                        A(sample_counter, c) = tmp.data(data_offset + s, c);
                        means(c) += tmp.data(data_offset + s, c);
                    }
                }

                sample_counter++;
                //GDEBUG("Sample counter = %d/%d\n", sample_counter, total_samples);
            }
        }

        GDEBUG_STREAM("Finished calculating A and means");

        //Subtract off mean
        for (size_t c = 0; c < channels; c++) {
            for (size_t s = 0; s < total_samples; s++) {
                A(s, c) -= means(c) / std::complex<float>(total_samples, 0);
            }
        }

        GDEBUG_STREAM("Finished subtracting mean");

        //Collected data for temp matrix, now let's calculate SVD coefficients

        pca_coefficients_[location] = new hoNDKLT < std::complex<float> > ;
        hoNDKLT< std::complex<float> >* VT = pca_coefficients_[location];

        GDEBUG_STREAM("Preparing VT");
        //We will create a new matrix that explicitly preserves the uncombined channels
        if (uncombined_channels_.size())
        {
            std::vector<size_t> untransformed(uncombined_channels_.size());
            for (size_t un = 0; un < uncombined_channels_.size(); un++)
            {
                untransformed[un] = uncombined_channels_[un];
            }

            VT->prepare(A, (size_t)1, untransformed, (size_t)0, false);

        }
        else
        {
            VT->prepare(A, (size_t)1, (size_t)0, false);
        }
        GDEBUG_STREAM("Finished VT->prepare")
    }

    void PCACoilGadget::do_pca(mrd::Acquisition& acq)
    {
        auto location = get_location(acq);

        hoNDArray<std::complex<float>> data_out(acq.data.dimensions());

        if (pca_coefficients_[location] != 0)
        {
            pca_coefficients_[location]->transform(acq.data, data_out, 1);
        }

        acq.data = data_out;
    }

    GADGETRON_GADGET_EXPORT(PCACoilGadget)
}
