/*
 * PCACoilGadget.cpp
 *
 *  Created on: Dec 13, 2011
 *      Author: Michael S. Hansen
 */

#include "PCACoilGadget.h"
#include "GadgetIsmrmrdReadWrite.h"
#include "hoArmadillo.h"
#include "hoNDArray_elemwise.h"

#ifdef HAVE_MKL
#include "mkl.h"
#endif

namespace Gadgetron {

  PCACoilGadget::PCACoilGadget()
    : max_buffered_profiles_(100)
    , samples_to_use_(16)
  {
    // There is a bug in the MKL SVD when running in multi-threaded mode.
    // Set the number of threads to 1 in this gadget.
#ifdef HAVE_MKL
    mkl_set_num_threads(1);
#endif

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
    bool is_last_scan_in_slice = (ISMRMRD::FlagBit(ISMRMRD::ACQ_LAST_IN_SLICE).isSet(m1->getObjectPtr()->flags));
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

	//For some sequences there is so little data, we should just use it all.
	if (profiles_available < 16) {
	  samples_to_use = samples_per_profile;
	}

	int total_samples = samples_to_use*profiles_available;

	std::vector<unsigned int> dims(2);
	dims[0] = channels;dims[1] = total_samples;

	hoNDArray< std::complex<float> > A;
	try{ A.create(&dims); }
	catch (std::runtime_error & err){
	  GADGET_DEBUG1("Unable to create array for PCA calculation\n");
	  return GADGET_FAIL;
	}

	std::complex<float>* A_ptr = A.get_data_ptr();
	unsigned int sample_counter = 0;

	unsigned int data_offset = 0;
	if (m1->getObjectPtr()->center_sample >= (samples_to_use>>1)) {
	  data_offset = m1->getObjectPtr()->center_sample - (samples_to_use>>1);
	}

	//GADGET_DEBUG2("Data offset = %d\n", data_offset);
 
	hoNDArray<std::complex<float> > means;
	std::vector<unsigned int> means_dims; means_dims.push_back(channels);

	try{means.create(&means_dims);}
	catch (std::runtime_error& err){
	  GADGET_DEBUG1("Unable to create temporary stoorage for mean values\n");
	  return GADGET_FAIL;
	}

	means.fill(std::complex<float>(0.0f,0.0f));

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
	    for (unsigned int c = 0; c < channels; c++) {
	      //We use the conjugate of the data so that the output VT of the SVD is the actual PCA coefficient matrix
	      A_ptr[c + sample_counter*channels] = d[c*samples_per_profile + data_offset + s];
	      means_ptr[c] += d[c*samples_per_profile + data_offset + s];
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

	std::vector<unsigned int> VT_dims;
	VT_dims.push_back(channels);
	VT_dims.push_back(channels);
	pca_coefficients_[location] = new hoNDArray< std::complex<float> >;
	hoNDArray< std::complex<float> >* VT = pca_coefficients_[location];

	try {VT->create(&VT_dims);}
	catch (std::runtime_error& err){
	  GADGET_DEBUG_EXCEPTION(err,"Failed to create array for VT\n");
	  return GADGET_FAIL;
	}
	
	arma::cx_fmat Am = as_arma_matrix(&A);
	arma::cx_fmat Vm = as_arma_matrix(VT);
	arma::cx_fmat Um;
	arma::fvec Sv;

	if( !arma::svd_econ(Um,Sv,Vm,Am.st(),'r') ){
	  GADGET_DEBUG1("Failed to compute SVD\n");
	  return GADGET_FAIL;
	}

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
      //GADGET_DEBUG1("Not buffering anymore\n");
      GadgetContainerMessage< hoNDArray< std::complex<float> > >* m3 =
	new GadgetContainerMessage< hoNDArray< std::complex<float> > >;

      try{m3->getObjectPtr()->create(m2->getObjectPtr()->get_dimensions().get()); }
      catch (std::runtime_error& err){
	GADGET_DEBUG_EXCEPTION(err,"Unable to create storage for PCA coils\n");
	m3->release();
	return GADGET_FAIL;
      }

      if (pca_coefficients_[location] != 0) {	
	arma::cx_fmat am3 = as_arma_matrix(m3->getObjectPtr());
	arma::cx_fmat am2 = as_arma_matrix(m2->getObjectPtr());
	arma::cx_fmat aPca = as_arma_matrix(pca_coefficients_[location]);
	am3 = am2*aPca;
      }

      m1->cont(m3);
      
      //In case there are trajectories attached. 
      m3->cont(m2->cont());
      m2->cont(0);

      m2->release();

      if (this->next()->putq(m1) < 0) {
	GADGET_DEBUG1("Unable to put message on Q");
	return GADGET_FAIL;
      }
    }
    return GADGET_OK;
  }

  GADGET_FACTORY_DECLARE(PCACoilGadget)
}
