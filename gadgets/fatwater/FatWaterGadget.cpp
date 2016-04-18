#include "FatWaterGadget.h"
#include "fatwater.h"
#include <ismrmrd/xml.h>

namespace Gadgetron{

    FatWaterGadget::FatWaterGadget()
    {
        
    }
    
    int FatWaterGadget::process_config(ACE_Message_Block* mb)
    {
        ISMRMRD::IsmrmrdHeader hdr;
        ISMRMRD::deserialize(mb->rd_ptr(), hdr);

        if (hdr.sequenceParameters.is_present()) {
            if (hdr.sequenceParameters->TE.is_present()) {
                for (auto& te: *(hdr.sequenceParameters->TE)) {
                    this->echoTimes_.push_back(te);
                }
            } else {
                GERROR("No echo times found in sequence parameters\n");
                return GADGET_FAIL;
            }
        } else {
            GERROR("Sequence parameters are required to do water fat seperations\n");
            return GADGET_FAIL;
        }

        for (auto& te: echoTimes_) {
            GDEBUG("Echo time: %f\n", te);
        }

        if (hdr.acquisitionSystemInformation.is_present() && hdr.acquisitionSystemInformation->systemFieldStrength_T.is_present()) {
            this->fieldStrength_ = *(hdr.acquisitionSystemInformation->systemFieldStrength_T);
            GDEBUG("Field strength: %f\n", this->fieldStrength_);
        } else {
            GERROR("Field strength not defined. Required for fat-water seperation\n");
            return GADGET_FAIL;
        }
              

        
        return GADGET_OK;
    }
    
    int FatWaterGadget::process( GadgetContainerMessage<IsmrmrdImageArray>* m1)
    {

        GDEBUG("In FW process\n");
    
        //Grab a reference to the buffer containing the imaging data
        IsmrmrdImageArray & imagearr = *m1->getObjectPtr();

        //7D, fixed order [X, Y, Z, CHA, N, S, LOC]
        uint16_t X = imagearr.data_.get_size(0);
        uint16_t Y = imagearr.data_.get_size(1);
        uint16_t Z = imagearr.data_.get_size(2);
        uint16_t CHA = imagearr.data_.get_size(3);
        uint16_t N = imagearr.data_.get_size(4);
        uint16_t S = imagearr.data_.get_size(5);
        uint16_t LOC = imagearr.data_.get_size(6);

        Gadgetron::FatWaterParameters p;
        
        p.echoTimes_ = this->echoTimes_;
        p.fieldStrengthT_ = this->fieldStrength_;

        if (p.echoTimes_.size() < 3) {
            GERROR("At least 3 echo times are required for fw separation\n");
            return GADGET_FAIL;
        }

        if (p.fieldStrengthT_ < 0.5) {
            GERROR("Water fat separation not possible at %f T\n", p.fieldStrengthT_);
            return GADGET_FAIL;
        }
            

        Gadgetron::FatWaterAlgorithm a;
        
        Gadgetron::ChemicalSpecies w("water");
        w.ampFreq_.push_back(std::make_pair(1.0, 0.0));

        Gadgetron::ChemicalSpecies f("fat");
        f.ampFreq_.push_back(std::make_pair(0.087,-3.8));
        f.ampFreq_.push_back(std::make_pair(0.693,-3.4));
        f.ampFreq_.push_back(std::make_pair(0.128,-2.6));
        f.ampFreq_.push_back(std::make_pair(0.004,-1.94));
        f.ampFreq_.push_back(std::make_pair(0.039,-0.39));
        f.ampFreq_.push_back(std::make_pair(0.048,0.60));
        
        a.species_.push_back(w);
        a.species_.push_back(f);

        try {
            hoNDArray< std::complex<float> > wfimages = Gadgetron::fatwater_separation(imagearr.data_, p, a);
        } catch (...) {
            GERROR("Error caught while doing water fat separation\n");
            m1->release();
            return GADGET_FAIL;
        }
        
        //Pass the image down the chain
        if (this->next()->putq(m1) < 0) {
            return GADGET_FAIL;
        }

        return GADGET_OK;  

    }

    GADGET_FACTORY_DECLARE(FatWaterGadget)
}
