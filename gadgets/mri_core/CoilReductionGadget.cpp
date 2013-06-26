/*
* CoilReductionGadget.cpp
*
*  Created on: Dec 5, 2011
*      Author: hansenms
*/

#include "CoilReductionGadget.h"
#include "GadgetIsmrmrdReadWrite.h"

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>

namespace Gadgetron{

    CoilReductionGadget::CoilReductionGadget() {
    }

    CoilReductionGadget::~CoilReductionGadget() {
    }

    int CoilReductionGadget::process_config(ACE_Message_Block *mb)
    {
        boost::shared_ptr<ISMRMRD::ismrmrdHeader> cfg = parseIsmrmrdXMLHeader(std::string(mb->rd_ptr()));

        coils_in_ = cfg->acquisitionSystemInformation().get().receiverChannels().get();

        boost::shared_ptr<std::string> coil_mask = this->get_string_value("coil_mask");

        if (coil_mask->compare(std::string("")) == 0) {
            int coils_out = this->get_int_value("coils_out");
            if (coils_out <= 0) {
                GADGET_DEBUG2("Invalid number of output coils %d\n", coils_out);
                return GADGET_FAIL;
            }
            coil_mask_ = std::vector<unsigned short>(coils_out,1);
        } else {
            std::vector<std::string> chm;
            boost::split(chm, *coil_mask, boost::is_any_of(" "));
            for (unsigned int i = 0; i < chm.size(); i++) {
                std::string ch = boost::algorithm::trim_copy(chm[i]);
                if (ch.size() > 0) {
                    unsigned int mv = static_cast<unsigned int>(ACE_OS::atoi(ch.c_str()));
                    //GADGET_DEBUG2("Coil mask value: %d\n", mv);
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
            GADGET_DEBUG1("Error configuring coils for coil reduction\n");
            return GADGET_FAIL;
        }

        coils_out_ = 0;
        for (unsigned int i = 0; i < coil_mask_.size(); i++) {
            if (coil_mask_[i]) coils_out_++;
        }

        GADGET_DEBUG2("Coil reduction from %d to %d\n", coils_in_, coils_out_);

        return GADGET_OK;
    }


    int CoilReductionGadget::process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader> *m1, GadgetContainerMessage<hoNDArray<std::complex<float> > > *m2)
    {
        std::vector<unsigned int> dims_out(2);
        dims_out[0] = m1->getObjectPtr()->number_of_samples;
        dims_out[1] = coils_out_;

        GadgetContainerMessage< hoNDArray<std::complex<float> > >* m3 =
            new GadgetContainerMessage< hoNDArray<std::complex<float> > >();

        try{ m3->getObjectPtr()->create(&dims_out);}
        catch (std::runtime_error &err){
            GADGET_DEBUG_EXCEPTION(err,"Unable to create storage for reduced dataset size\n");
            return GADGET_FAIL;
        }

        std::complex<float>* s = m2->getObjectPtr()->get_data_ptr();
        std::complex<float>* d = m3->getObjectPtr()->get_data_ptr();
        unsigned int samples =  m1->getObjectPtr()->number_of_samples;
        unsigned int coils_copied = 0;
        for (int c = 0; c < m1->getObjectPtr()->active_channels; c++) {
            if (c > coil_mask_.size()) {
                GADGET_DEBUG1("Fatal error, too many coils for coil mask\n");
                m3->release();
                return GADGET_FAIL;
            }
            if (coil_mask_[c]) {
                memcpy(d+coils_copied*samples,s+c*samples,sizeof(std::complex<float>)*samples);
                coils_copied++;
            }
        }

        m2->release();
        m1->cont(m3);
        m1->getObjectPtr()->active_channels = coils_out_;
	
        if( this->next()->putq(m1) < 0 ){
	  GADGET_DEBUG1("Failed to put message on queue\n");
	  return GADGET_FAIL;
	}
	
	return GADGET_OK;
    }

    GADGET_FACTORY_DECLARE(CoilReductionGadget)
}
