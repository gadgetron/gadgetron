/*
* PCACoilGadget.cpp
*
*  Created on: Dec 13, 2011
*      Author: Michael S. Hansen
*/

#include "PCACoilGadget.h"
#include "hoNDArray_elemwise.h"
#include "ismrmrd/xml.h"
#include "hoNDArray_fileio.h"
#include "hoNDKLT.h"
#include "hoNDArray_linalg.h"

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>

namespace Gadgetron {

    PCACoilGadget::PCACoilGadget()
        : max_buffered_profiles_(100)
        , samples_to_use_(16)
    {
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

    int PCACoilGadget::process_config(ACE_Message_Block *mb)
    {
        ISMRMRD::IsmrmrdHeader h;
        ISMRMRD::deserialize(mb->rd_ptr(), h);

        std::string uncomb_str = uncombined_channels_by_name.value();
        std::vector<std::string> uncomb;
        if (uncomb_str.size()) {
            GDEBUG("uncomb_str: %s\n", uncomb_str.c_str());
            boost::split(uncomb, uncomb_str, boost::is_any_of(","));
            for (unsigned int i = 0; i < uncomb.size(); i++) {
                std::string ch = boost::algorithm::trim_copy(uncomb[i]);
                if (h.acquisitionSystemInformation) {
                    for (size_t i = 0; i < h.acquisitionSystemInformation->coilLabel.size(); i++) {
                        if (ch == h.acquisitionSystemInformation->coilLabel[i].coilName) {
                            uncombined_channels_.push_back(i);//This assumes that the channels are sorted in the header
                            break;
                        }
                    }
                }
            }
        }

        present_uncombined_channels.value((int)uncombined_channels_.size());
        GDEBUG("Number of uncombined channels (present_uncombined_channels) set to %d\n", uncombined_channels_.size());

#ifdef USE_OMP
        omp_set_num_threads(1);
#endif // USE_OMP

        return GADGET_OK;
    }

    int PCACoilGadget::process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader> *m1, GadgetContainerMessage<hoNDArray<std::complex<float> > > *m2)
    {
        bool is_noise = m1->getObjectPtr()->isFlagSet(ISMRMRD::ISMRMRD_ACQ_IS_NOISE_MEASUREMENT);

        //We should not be receiving noise here
        if (is_noise) {
            m1->release();
            return GADGET_OK;
        }


        std::map<int, bool>::iterator it;
        int location = m1->getObjectPtr()->idx.slice;
        bool is_last_scan_in_slice = m1->getObjectPtr()->isFlagSet(ISMRMRD::ISMRMRD_ACQ_LAST_IN_SLICE);
        int samples_per_profile = m1->getObjectPtr()->number_of_samples;
        int channels = m1->getObjectPtr()->active_channels;

        it = buffering_mode_.find(location);

        bool is_buffering = true;
        //Do we have an entry for this location
        if (it != buffering_mode_.end()) {
            is_buffering = it->second;
        }
        else {
            //else make an entry. We will always start in buffering mode for a given location.
            buffering_mode_[location] = is_buffering;
        }

        if (is_buffering)
        {
            buffer_[location].push_back(m1);
            int profiles_available = buffer_[location].size();

            //Are we ready for calculating PCA
            if (is_last_scan_in_slice || (profiles_available >= max_buffered_profiles_))
            {

                //GDEBUG("Calculating PCA coefficients with %d profiles for %d coils\n", profiles_available, channels);
                int samples_to_use = samples_per_profile > samples_to_use_ ? samples_to_use_ : samples_per_profile;

                //For some sequences there is so little data, we should just use it all.
                if (profiles_available < 16) {
                    samples_to_use = samples_per_profile;
                }

                int total_samples = samples_to_use*profiles_available;

                std::vector<size_t> dims(2);
                dims[0] = total_samples; dims[1] = channels;

                hoNDArray< std::complex<float> > A;
                try{ A.create(dims); }
                catch (std::runtime_error & err){
                    GDEBUG("Unable to create array for PCA calculation\n");
                    return GADGET_FAIL;
                }

                std::complex<float>* A_ptr = A.get_data_ptr();
                size_t sample_counter = 0;

                size_t data_offset = 0;
                if (m1->getObjectPtr()->center_sample >= (samples_to_use >> 1)) {
                    data_offset = m1->getObjectPtr()->center_sample - (samples_to_use >> 1);
                }

                //GDEBUG("Data offset = %d\n", data_offset);

                hoNDArray<std::complex<float> > means;
                std::vector<size_t> means_dims; means_dims.push_back(channels);

                try{ means.create(means_dims); }
                catch (std::runtime_error& err){
                    GDEBUG("Unable to create temporary stoorage for mean values\n");
                    return GADGET_FAIL;
                }

                means.fill(std::complex<float>(0.0f, 0.0f));

                std::complex<float>* means_ptr = means.get_data_ptr();

                for (size_t p = 0; p < profiles_available; p++) {
                    GadgetContainerMessage<hoNDArray<std::complex<float> > >* m_tmp =
                        AsContainerMessage<hoNDArray< std::complex<float> > >(buffer_[location][p]->cont());

                    if (!m_tmp) {
                        GDEBUG("Fatal error, unable to recover data from data buffer (%d,%d)\n", p, profiles_available);
                        return GADGET_FAIL;
                    }

                    std::complex<float>* d = m_tmp->getObjectPtr()->get_data_ptr();

                    for (unsigned s = 0; s < samples_to_use; s++) {
                        for (size_t c = 0; c < channels; c++) {
                            bool uncombined_channel = std::find(uncombined_channels_.begin(), uncombined_channels_.end(), c) != uncombined_channels_.end();
                            if (uncombined_channel) {
                                A_ptr[sample_counter + c *total_samples] = std::complex<float>(0.0, 0.0);
                            }
                            else {
                                A_ptr[sample_counter + c *total_samples] = d[c*samples_per_profile + data_offset + s];
                                means_ptr[c] += d[c*samples_per_profile + data_offset + s];
                            }
                        }

                        sample_counter++;
                        //GDEBUG("Sample counter = %d/%d\n", sample_counter, total_samples);
                    }
                }

                //Subtract off mean
                for (size_t c = 0; c < channels; c++) {
                    for (size_t s = 0; s < total_samples; s++) {
                        A_ptr[s + c *total_samples] -= means_ptr[c] / std::complex<float>(total_samples, 0);
                    }
                }

                //Collected data for temp matrix, now let's calculate SVD coefficients

                std::vector<size_t> VT_dims;
                VT_dims.push_back(channels);
                VT_dims.push_back(channels);
                pca_coefficients_[location] = new hoNDKLT < std::complex<float> > ;
                hoNDKLT< std::complex<float> >* VT = pca_coefficients_[location];

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

                //Switch off buffering for this slice
                buffering_mode_[location] = false;

                //Now we should pump all the profiles that we have buffered back through the system
                for (size_t p = 0; p < profiles_available; p++) {
                    ACE_Message_Block* mb = buffer_[location][p];
                    if (inherited::process(mb) != GADGET_OK) {
                        GDEBUG("Failed to reprocess buffered data\n");
                        return GADGET_FAIL;
                    }
                }
                //Remove references in this buffer
                buffer_[location].clear();
            }
        }
        else {
            //GDEBUG("Not buffering anymore\n");
            GadgetContainerMessage< hoNDArray< std::complex<float> > >* m3 =
                new GadgetContainerMessage < hoNDArray< std::complex<float> > > ;

            try{ m3->getObjectPtr()->create(m2->getObjectPtr()->dimensions()); }
            catch (std::runtime_error& err){
                GEXCEPTION(err, "Unable to create storage for PCA coils\n");
                m3->release();
                return GADGET_FAIL;
            }

            if (pca_coefficients_[location] != 0)
            {
                pca_coefficients_[location]->transform(*(m2->getObjectPtr()), *(m3->getObjectPtr()), 1);
            }

            m1->cont(m3);

            //In case there are trajectories attached. 
            m3->cont(m2->cont());
            m2->cont(0);

            m2->release();

            if (this->next()->putq(m1) < 0) {
                GDEBUG("Unable to put message on Q");
                return GADGET_FAIL;
            }
        }
        return GADGET_OK;
    }

    GADGET_FACTORY_DECLARE(PCACoilGadget)
}
