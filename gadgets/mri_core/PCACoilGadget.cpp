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
#include "ismrmrd/xml.h"
#include "hoNDArray_fileio.h"

#include <ace/OS_NS_stdlib.h>
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
      ISMRMRD::IsmrmrdHeader h;
      ISMRMRD::deserialize(mb->rd_ptr(),h);

      if (h.acquisitionSystemInformation) {
	for (size_t i = 0; i < h.acquisitionSystemInformation->coilLabel.size(); i++) {
	  int coil_num = h.acquisitionSystemInformation->coilLabel[i].coilNumber;
	  channel_map_[h.acquisitionSystemInformation->coilLabel[i].coilName] = coil_num;
	}
      }
      
      
      boost::shared_ptr<std::string> uncomb_str = this->get_string_value("uncombined_channels_by_name");
      std::vector<std::string> uncomb;
      if (uncomb_str->size()) {
	GADGET_DEBUG2("uncomb_str: %s\n",  uncomb_str->c_str());
	boost::split(uncomb, *uncomb_str, boost::is_any_of(","));
	for (unsigned int i = 0; i < uncomb.size(); i++) {
	  std::string ch = boost::algorithm::trim_copy(uncomb[i]);
	  coil_map_type_::iterator it = channel_map_.find(ch);
	  if (it != channel_map_.end()) {
	    unsigned int channel_id = static_cast<unsigned int>(it->second);
	    GADGET_DEBUG2("Device channel: %s (%d)\n",  uncomb[i].c_str(), channel_id);
	    uncombined_channels_.push_back(channel_id);
	  }
	}
      }

      char val[32];
      sprintf(val,"%d",(int)uncombined_channels_.size());
      this->set_parameter("present_uncombined_channels",val);
      GADGET_DEBUG2("Number of uncombined channels (present_uncombined_channels) set to %d\n", uncombined_channels_.size());

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

                std::vector<size_t> dims(2);
                dims[0] = channels;dims[1] = total_samples;

                hoNDArray< std::complex<float> > A;
                try{ A.create(&dims); }
                catch (std::runtime_error & err){
                    GADGET_DEBUG1("Unable to create array for PCA calculation\n");
                    return GADGET_FAIL;
                }

                std::complex<float>* A_ptr = A.get_data_ptr();
                size_t sample_counter = 0;

                size_t data_offset = 0;
                if (m1->getObjectPtr()->center_sample >= (samples_to_use>>1)) {
                    data_offset = m1->getObjectPtr()->center_sample - (samples_to_use>>1);
                }

                //GADGET_DEBUG2("Data offset = %d\n", data_offset);

                hoNDArray<std::complex<float> > means;
                std::vector<size_t> means_dims; means_dims.push_back(channels);

                try{means.create(&means_dims);}
                catch (std::runtime_error& err){
                    GADGET_DEBUG1("Unable to create temporary stoorage for mean values\n");
                    return GADGET_FAIL;
                }

                means.fill(std::complex<float>(0.0f,0.0f));

                std::complex<float>* means_ptr = means.get_data_ptr();
                for (size_t p = 0; p < profiles_available; p++) {
                    GadgetContainerMessage<hoNDArray<std::complex<float> > >* m_tmp =
                        AsContainerMessage<hoNDArray< std::complex<float> > >(buffer_[location][p]->cont());

                    if (!m_tmp) {
                        GADGET_DEBUG2("Fatal error, unable to recover data from data buffer (%d,%d)\n", p, profiles_available);
                        return GADGET_FAIL;
                    }

                    std::complex<float>* d = m_tmp->getObjectPtr()->get_data_ptr();

		      for (unsigned s = 0; s < samples_to_use; s++) {
			for (size_t c = 0; c < channels; c++) {
			  bool uncombined_channel = std::find(uncombined_channels_.begin(),uncombined_channels_.end(), c) != uncombined_channels_.end();
			  //We use the conjugate of the data so that the output VT of the SVD is the actual PCA coefficient matrix
			  if (uncombined_channel) {
			    A_ptr[c + sample_counter*channels] = std::complex<float>(0.0,0.0);
			  } else {
			    A_ptr[c + sample_counter*channels] = d[c*samples_per_profile + data_offset + s];
			    means_ptr[c] += d[c*samples_per_profile + data_offset + s];
			  }
			}
			
			sample_counter++;
			//GADGET_DEBUG2("Sample counter = %d/%d\n", sample_counter, total_samples);
		      }
                }

                //Subtract off mean
                for (size_t c = 0; c < channels; c++) {
                    for (size_t s = 0; s < total_samples; s++) {
                        A_ptr[c + s*channels] -=  means_ptr[c]/std::complex<float>(total_samples,0);
                    }
                }

                //Collected data for temp matrix, now let's calculate SVD coefficients

                std::vector<size_t> VT_dims;
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
		
		//We will create a new matrix that explicitly preserves the uncombined channels
		if (uncombined_channels_.size()) {
		  hoNDArray< std::complex<float> >* VT_new = new hoNDArray< std::complex<float> >;
		  try {VT_new->create(&VT_dims);}
		  catch (std::runtime_error& err){
                    GADGET_DEBUG_EXCEPTION(err,"Failed to create array for VT (new)\n");
                    return GADGET_FAIL;
		  }

		  arma::cx_fmat Vm_new = as_arma_matrix(VT_new);

		  size_t uncomb_count = 0;
		  size_t comb_count = 0;
		  for (size_t c = 0; c < Vm_new.n_cols; c++) {
		    bool uncombined_channel = std::find(uncombined_channels_.begin(),uncombined_channels_.end(), c) != uncombined_channels_.end();
		    if (uncombined_channel) {
		      for (size_t r = 0; r < Vm_new.n_rows; r++) {
			if (r == c) {
			  Vm_new(r,uncomb_count) = 1;
			} else {
			  Vm_new(r,uncomb_count) = 0;
			}
		      }
		      uncomb_count++;
		    } else {
		      for (size_t r = 0; r < Vm_new.n_rows; r++) { 
			bool uncombined_channel_row = std::find(uncombined_channels_.begin(),uncombined_channels_.end(), r) != uncombined_channels_.end();
			if (uncombined_channel_row) {
			  Vm_new(r,comb_count+uncombined_channels_.size()) = 0;
			} else {
			  Vm_new(r,comb_count+uncombined_channels_.size()) = Vm(r,c);
			}
		      }
		      comb_count++;
		    }
		  } 
		  GADGET_DEBUG2("uncomb_count = %d, comb_count = %d\n", uncomb_count, comb_count);

		  //Delete the old one and set the new one
		  delete pca_coefficients_[location];
		  pca_coefficients_[location] = VT_new;
		}


                //Switch off buffering for this slice
                buffering_mode_[location] = false;

                //Now we should pump all the profiles that we have buffered back through the system
                for (size_t p = 0; p < profiles_available; p++) {
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
