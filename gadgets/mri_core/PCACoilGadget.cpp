/*
 * PCACoilGadget.cpp
 *
 *  Created on: Dec 13, 2011
 *      Author: Michael S. Hansen
 */

#include "PCACoilGadget.h"
#include "hoNDArray_fileio.h"
#include "matrix_decomposition.h"

PCACoilGadget::PCACoilGadget()
 : max_buffered_profiles_(100)
 , samples_to_use_(16)
{

}

PCACoilGadget::~PCACoilGadget()
{

}


int PCACoilGadget::process_config(ACE_Message_Block *mb)
{
	return GADGET_OK;
}



int PCACoilGadget::process(GadgetContainerMessage<GadgetMessageAcquisition> *m1, GadgetContainerMessage<hoNDArray<std::complex<float> > > *m2)
{
	std::map<int, bool>::iterator it;

	int location = m1->getObjectPtr()->idx.slice;
	bool is_last_scan_in_slice = (m1->getObjectPtr()->flags & GADGET_FLAG_LAST_ACQ_IN_SLICE);
	int samples_per_profile = m1->getObjectPtr()->samples;
	int channels = m1->getObjectPtr()->channels;

	it = buffering_mode_.find(location);

	bool is_buffering = true;
	//Do we have an entry for this location
	if (it != buffering_mode_.end()) {
		is_buffering = it->second;
	} else {
		//else make an entry. We will always start in buffering mode for a given location.
		buffering_mode_[location] = is_buffering;
	}

	if (is_buffering) {
		buffer_[location].push_back(m1);

		int profiles_available = buffer_[location].size();

		//Are we ready for calculating PCA
		if (is_last_scan_in_slice || (profiles_available >= max_buffered_profiles_)) {

			GADGET_DEBUG2("Calculating PCA coefficients with %d profiles for %d coils\n", profiles_available, channels);
			int samples_to_use = samples_per_profile > samples_to_use_ ? samples_to_use_ : samples_per_profile;
			int total_samples = samples_to_use*profiles_available;

			std::vector<unsigned int> dims(2);
			dims[0] = total_samples;dims[1] = channels;

			hoNDArray< std::complex<float> > A;
			if (!A.create(&dims)) {
				GADGET_DEBUG1("Unable to create array for PCA calculation\n");
			}

			std::complex<float>* A_ptr = A.get_data_ptr();
			unsigned int sample_counter = 0;
			unsigned int data_offset = (samples_per_profile -samples_to_use)>>1;

			for (unsigned int p = 0; p < profiles_available; p++) {
				GadgetContainerMessage<hoNDArray<std::complex<float> > >* m_tmp =
						AsContainerMessage<hoNDArray< std::complex<float> > >(buffer_[location][p]->cont());

				if (!m_tmp) {
					GADGET_DEBUG2("Fatal error, unable to recover data from data buffer (%d,%d)\n", p, profiles_available);
					return GADGET_FAIL;
				}

				std::complex<float>* d = m_tmp->getObjectPtr()->get_data_ptr();

				for (unsigned s = 0; s < samples_to_use; s++) {
					for (unsigned int c = 0; c < channels; c++) {
						A_ptr[c*total_samples + sample_counter] =
								d[c*samples_per_profile + data_offset + s];
					}
					sample_counter++;
					GADGET_DEBUG2("Sample counter = %d/%d\n", sample_counter, total_samples);
				}
			}

			//Collected data for temp matrix, now let's calculate SVD coefficients

			write_nd_array(&A,"A.cplx");
			hoNDArray_svd< std::complex<float> >(&A, 0 , 0 , 0);
			write_nd_array(&A,"A2.cplx");



			//Switch off buffering for this slice
			buffering_mode_[location] = false;
		}
	}


	return GADGET_OK;//this->next()->putq(m1);
}

GADGET_FACTORY_DECLARE(PCACoilGadget)

