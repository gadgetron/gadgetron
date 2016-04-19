#include "FatWaterGadget.h"
#include "fatwater.h"
#include "mri_core_def.h"
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
            //This should return the images
            hoNDArray< std::complex<float> > wfimages = Gadgetron::fatwater_separation(imagearr.data_, p, a);

            //Now let's make an image array to return the f/w images + any additional recon products
            //Right now, we will be using a copy of the original data
            //TODO: Remove this data copy
            //wfimages = imagearr.data_;

            uint16_t n_images = wfimages.get_size(4); 
            uint16_t s_images = wfimages.get_size(5); //S-dimention is the image dimension
            uint16_t loc_images = wfimages.get_size(6);

            if (n_images != N || loc_images != LOC) {
                GERROR("Wrong number of N or LOC images received from fat water seperation\n");
                m1->release();
                return GADGET_FAIL;
            }

            
            auto m2 = new GadgetContainerMessage<IsmrmrdImageArray>();
            m2->getObjectPtr()->data_ = wfimages; //It is a copy, but worth it for simplicty
            m2->getObjectPtr()->headers_.create(n_images, s_images, loc_images);
            for (uint16_t loc = 0; loc < loc_images; loc++) {
                for (uint16_t s = 0; s < s_images; s++) {
                    for (uint16_t n = 0; n < n_images; n++) {
                        m2->getObjectPtr()->headers_[loc*s_images*n_images + s*n_images + n] = m1->getObjectPtr()->headers_[loc*s_images*n_images + 0*n_images + n];

                        m2->getObjectPtr()->headers_[loc*s_images*n_images + s*n_images + n].image_series_index =
                            m2->getObjectPtr()->headers_[loc*s_images*n_images + s*n_images + n].image_series_index + 100;
                        
                        ISMRMRD::MetaContainer meta = m1->getObjectPtr()->meta_[loc*n_images*s_images + n];
                        //TODO: These sepcies and image type specifiers should come from the toolbox
                        if (s < a.species_.size()) {
                            if (a.species_[s].name_ == "water") {
                                meta.set(GADGETRON_DATA_ROLE, GADGETRON_IMAGE_WATER);
                            } else if (a.species_[s].name_ == "fat") {
                                meta.set(GADGETRON_DATA_ROLE, GADGETRON_IMAGE_FAT);
                            }
                        } else {
                            //TODO: What to call these images
                        }
                        m2->getObjectPtr()->meta_.push_back(meta);
                    }
                }
            }
            
            if (this->next()->putq(m2) < 0) {
                m1->release();
                m2->release();
                return GADGET_FAIL;
            }
            
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
