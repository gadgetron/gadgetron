/*
 * PCACoilGadget.cpp
 *
 *  Created on: Dec 13, 2011
 *      Author: Michael S. Hansen
 */

#include "PCACoilGadget.h"
#include "hoNDArray_fileio.h"
#include "matrix_vector_op.h"
#include "matrix_decomposition.h"

PCACoilGadget::PCACoilGadget()
 : max_buffered_profiles_(100)
 , samples_to_use_(16)
{

}

PCACoilGadget::~PCACoilGadget()
{

	std::map<int, hoNDArray<std::complex<float> >* >::iterator it;
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
	return GADGET_OK;
}



int PCACoilGadget::process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader> *m1, GadgetContainerMessage<hoNDArray<std::complex<float> > > *m2)
{
	std::map<int, bool>::iterator it;

	int location = m1->getObjectPtr()->idx.slice;
	bool is_last_scan_in_slice = (ISMRMRD::FlagBit(ISMRMRD::LAST_IN_SLICE).isSet(m1->getObjectPtr()->flags));
	int samples_per_profile = m1->getObjectPtr()->number_of_samples;
	int channels = m1->getObjectPtr()->active_channels;

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

			//GADGET_DEBUG2("Calculating PCA coefficients with %d profiles for %d coils\n", profiles_available, channels);
			int samples_to_use = samples_per_profile > samples_to_use_ ? samples_to_use_ : samples_per_profile;
			int total_samples = samples_to_use*profiles_available;

			std::vector<unsigned int> dims(2);
			dims[0] = channels;dims[1] = total_samples;

			hoNDArray< std::complex<float> > A;
			if (!A.create(&dims)) {
				GADGET_DEBUG1("Unable to create array for PCA calculation\n");
			}

			std::complex<float>* A_ptr = A.get_data_ptr();
			unsigned int sample_counter = 0;
			unsigned int data_offset = (samples_per_profile -samples_to_use)>>1;

			hoNDArray<std::complex<float> > means;
			std::vector<unsigned int> means_dims; means_dims.push_back(channels);

			if (!means.create(&means_dims)) {
				GADGET_DEBUG1("Unable to create temporary stoorage for mean values\n");
				return -1;
			}

			means.clear(std::complex<float>(0.0f,0.0f));

			std::complex<float>* means_ptr = means.get_data_ptr();
			for (unsigned int p = 0; p < profiles_available; p++) {
				GadgetContainerMessage<hoNDArray<std::complex<float> > >* m_tmp =
						AsContainerMessage<hoNDArray< std::complex<float> > >(buffer_[location][p]->cont());

				if (!m_tmp) {
					GADGET_DEBUG2("Fatal error, unable to recover data from data buffer (%d,%d)\n", p, profiles_available);
					return GADGET_FAIL;
				}

				std::complex<float>* d = m_tmp->getObjectPtr()->get_data_ptr();

				for (unsigned s = 0; s < samples_to_use; s++) {
					std::complex<float> mean(0.0,0.0);
					for (unsigned int c = 0; c < channels; c++) {

						//We use the conjugate of the data so that the output VT of the SVD is the actual PCA coefficient matrix\
						A_ptr[c + sample_counter*channels] =
								conj(d[c*samples_per_profile + data_offset + s]);

						means_ptr[c] += conj(d[c*samples_per_profile + data_offset + s]);
					}

					sample_counter++;
					//GADGET_DEBUG2("Sample counter = %d/%d\n", sample_counter, total_samples);
				}
			}

			//Subtract off mean
			for (unsigned int c = 0; c < channels; c++) {
				for (unsigned int s = 0; s < total_samples; s++) {
					A_ptr[c + s*channels] -=  means_ptr[c]/std::complex<float>(total_samples,0);
				}
			}


			//Collected data for temp matrix, now let's calculate SVD coefficients

			//write_nd_array(&A,"A.cplx");

			std::vector<unsigned int> S_dims; S_dims.push_back(channels);
			hoNDArray<float> S;
			if (!S.create(&S_dims)) {
				GADGET_DEBUG1("Failed to create array for singular values\n");
				return GADGET_FAIL;
			}

			std::vector<unsigned int> VT_dims;
			VT_dims.push_back(channels);
			VT_dims.push_back(channels);
			pca_coefficients_[location] = new hoNDArray< std::complex<float> >;
			hoNDArray< std::complex<float> >* VT = pca_coefficients_[location];

			if (!VT->create(&VT_dims)) {
				GADGET_DEBUG1("Failed to create array for VT\n");
				return GADGET_FAIL;
			}

			/*
			std::vector<unsigned int> U_dims;
			U_dims.push_back(channels);
			U_dims.push_back(total_samples);
			hoNDArray< std::complex<float> > U;
			if (!U.create(&U_dims)) {
				GADGET_DEBUG1("Failed to create array for U\n");
				return GADGET_FAIL;
			}
			*/

			//We don't need to calculate U in this case.
			//if (hoNDArray_svd< std::complex<float>, float >(&A, &U , &S , &VT) != 0) {
			if (hoNDArray_svd< std::complex<float>, float >(&A, 0 , &S , VT) != 0) {
				GADGET_DEBUG1("SVD failed\n");
				return GADGET_FAIL;
			}



			//write_nd_array(&S,"S.real");
			//write_nd_array(&U,"U.cplx");
			//write_nd_array(VT,"VT.cplx");

			//Switch off buffering for this slice
			buffering_mode_[location] = false;

			//Now we should pump all the profiles that we have buffered back through the system
			for (unsigned int p = 0; p < profiles_available; p++) {
				ACE_Message_Block* mb = buffer_[location][p];
				if (inherited::process(mb) != GADGET_OK) {
					GADGET_DEBUG1("Failed to reprocess buffered data\n");
					return GADGET_FAIL;
				}
			}
			//Remove references in this buffer
			buffer_[location].clear();


		}
	} else {
		GadgetContainerMessage< hoNDArray< std::complex<float> > >* m3 =
				new GadgetContainerMessage< hoNDArray< std::complex<float> > >;

		if (!m3->getObjectPtr()->create(m2->getObjectPtr()->get_dimensions().get())) {
			GADGET_DEBUG1("Unable to create storage for PCA coils\n");
			m3->release();
			return GADGET_FAIL;
		}

		if (pca_coefficients_[location] != 0) {
			std::complex<float> alpha(1.0,0.0);
			std::complex<float> beta(0.0,0.0);
			if (hoNDArray_gemm( pca_coefficients_[location], m2->getObjectPtr(), alpha,  m3->getObjectPtr(), beta) < 0) {
				GADGET_DEBUG1("Failed to apply PCA coefficients\n");
				return GADGET_FAIL;
			}
		}

		m1->cont(m3);
		m2->release();

		if (this->next()->putq(m1) < 0) {
			GADGET_DEBUG1("Unable to put message on Q");
			return GADGET_FAIL;
		}
	}

	return GADGET_OK;//
}

GADGET_FACTORY_DECLARE(PCACoilGadget)

