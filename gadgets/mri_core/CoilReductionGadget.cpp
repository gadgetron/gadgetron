/*
* CoilReductionGadget.cpp
*
*  Created on: Dec 5, 2011
*      Author: hansenms
*/

#include "CoilReductionGadget.h"
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>
#include "ismrmrd/xml.h"

namespace Gadgetron{

    CoilReductionGadget::CoilReductionGadget() {
    }

    CoilReductionGadget::~CoilReductionGadget() {
    }

    int CoilReductionGadget::process_config(ACE_Message_Block *mb)
    {
      ISMRMRD::IsmrmrdHeader h;
      ISMRMRD::deserialize(mb->rd_ptr(),h);
      
      coils_in_ = h.acquisitionSystemInformation->receiverChannels ? *h.acquisitionSystemInformation->receiverChannels : 128;

      std::string coil_mask_int = coil_mask.value();

      if (coil_mask_int.compare(std::string("")) == 0) {
	if (coils_out.value() <= 0) {
	  GDEBUG("Invalid number of output coils %d\n", coils_out.value());
	  return GADGET_FAIL;
	}
	coil_mask_ = std::vector<unsigned short>(coils_out.value(),1);
      } else {
	std::vector<std::string> chm;
	boost::split(chm, coil_mask_int, boost::is_any_of(" "));
	for (size_t i = 0; i < chm.size(); i++) {
	  std::string ch = boost::algorithm::trim_copy(chm[i]);
	  if (ch.size() > 0) {
	    size_t mv = static_cast<size_t>(std::stoi(ch));
	    //GDEBUG("Coil mask value: %d\n", mv);
	    if (mv > 0) {
	      coil_mask_.push_back(1);
	    } else {
	      coil_mask_.push_back(0);
	    }
	  }
	}
      }
      
      while (coil_mask_.size() < coils_in_) coil_mask_.push_back(0);
      while (coil_mask_.size() > coils_in_) coil_mask_.pop_back();
      
      if (coil_mask_.size() != coils_in_) {
	GDEBUG("Error configuring coils for coil reduction\n");
	return GADGET_FAIL;
      }
      
      coils_out_ = 0;
      for (size_t i = 0; i < coil_mask_.size(); i++) {
	if (coil_mask_[i]) coils_out_++;
      }
      
      GDEBUG("Coil reduction from %d to %d\n", coils_in_, coils_out_);
      
      return GADGET_OK;
    }


    int CoilReductionGadget::process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader> *m1, GadgetContainerMessage<hoNDArray<std::complex<float> > > *m2)
    {
        std::vector<size_t> dims_out(2);
        dims_out[0] = m1->getObjectPtr()->number_of_samples;
        dims_out[1] = coils_out_;

        GadgetContainerMessage< hoNDArray<std::complex<float> > >* m3 =
            new GadgetContainerMessage< hoNDArray<std::complex<float> > >();

        try{ m3->getObjectPtr()->create(dims_out);}
        catch (std::runtime_error &err){
            GEXCEPTION(err,"Unable to create storage for reduced dataset size\n");
            return GADGET_FAIL;
        }

        std::complex<float>* s = m2->getObjectPtr()->get_data_ptr();
        std::complex<float>* d = m3->getObjectPtr()->get_data_ptr();
        size_t samples =  m1->getObjectPtr()->number_of_samples;
        size_t coils_copied = 0;
        for (int c = 0; c < m1->getObjectPtr()->active_channels; c++) {
            if (c > coil_mask_.size()) {
                GDEBUG("Fatal error, too many coils for coil mask\n");
                m3->release();
                return GADGET_FAIL;
            }
            if (coil_mask_[c]) {
                memcpy(d+coils_copied*samples,s+c*samples,sizeof(std::complex<float>)*samples);
                coils_copied++;
            }
        }

        m1->cont(m3);
	
	//In case trajectories are attached
	m3->cont(m2->cont());
	m2->cont(0);

        m2->release();

        m1->getObjectPtr()->active_channels = coils_out_;
	
        if( this->next()->putq(m1) < 0 ){
	  GDEBUG("Failed to put message on queue\n");
	  return GADGET_FAIL;
	}
	
	return GADGET_OK;
    }

    GADGET_FACTORY_DECLARE(CoilReductionGadget)
}
