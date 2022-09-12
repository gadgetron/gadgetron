#include "FatWaterGadget.h"
#include "fatwater.h"
#include "mri_core_def.h"
#include <ismrmrd/xml.h>
#include "correct_frequency_shift.h"

namespace Gadgetron {

    using namespace std::complex_literals;

    FatWaterGadget::FatWaterGadget() {

    }

    int FatWaterGadget::process_config(ACE_Message_Block *mb) {
        ISMRMRD::IsmrmrdHeader hdr;
        ISMRMRD::deserialize(mb->rd_ptr(), hdr);

        if (hdr.sequenceParameters.is_present()) {
            if (hdr.sequenceParameters->TE.is_present()) {
                for (auto &te: *(hdr.sequenceParameters->TE)) {
                    this->echoTimes_.push_back(te * 1e-3);
                }
            } else {
                GERROR("No echo times found in sequence parameters\n");
                return GADGET_FAIL;
            }
        } else {
            GERROR("Sequence parameters are required to do water fat separations\n");
            return GADGET_FAIL;
        }

        for (auto &te: echoTimes_) {
            GDEBUG("Echo time: %f\n", te);
        }

        if (hdr.acquisitionSystemInformation.is_present() &&
            hdr.acquisitionSystemInformation->systemFieldStrength_T.is_present()) {
            this->fieldStrength_ = *(hdr.acquisitionSystemInformation->systemFieldStrength_T);
            GDEBUG("Field strength: %f\n", this->fieldStrength_);
        } else {
            GERROR("Field strength not defined. Required for fat-water separation\n");
            return GADGET_FAIL;
        }


//        config.frequency_range = {range_frequency_offset.value()[0],range_frequency_offset.value()[1]};
        float omega = 1/(this->echoTimes_[1]-this->echoTimes_[0]);
        config.frequency_range = {-omega,omega};
        config.lambda = regularization_lambda;
        config.lambda_extra = regularization_offset;
        config.number_of_frequency_samples = number_of_frequency_offsets;
        config.r2_range = {range_r2star.value()[0],range_r2star.value()[1]};
        config.number_of_r2_samples = number_of_r2stars;
        config.number_of_r2_fine_samples = number_of_r2stars_fine;
        config.do_gradient_descent = do_gradient_descent;
        config.downsamples = downsample_data;


        return GADGET_OK;
    }

    int FatWaterGadget::process(GadgetContainerMessage<IsmrmrdImageArray> *m1) {

        GDEBUG("In FW process\n");

        //Grab a reference to the buffer containing the imaging data
        IsmrmrdImageArray &imagearr = *m1->getObjectPtr();

        //7D, fixed order [X, Y, Z, CHA, N, S, LOC]
        uint16_t X = imagearr.data_.get_size(0);
        uint16_t Y = imagearr.data_.get_size(1);
        uint16_t Z = imagearr.data_.get_size(2);
        uint16_t CHA = imagearr.data_.get_size(3);
        uint16_t N = imagearr.data_.get_size(4);
        uint16_t S = imagearr.data_.get_size(5);
        uint16_t LOC = imagearr.data_.get_size(6);

        FatWater::Parameters parameters;

        parameters.echo_times_s = this->echoTimes_;

        parameters.field_strength_T = this->fieldStrength_;
        parameters.sample_time_us = sample_time_us;

        if (echoTimes_.size() < 3) {
            GERROR("At least 3 echo times are required for fw separation\n");
            return GADGET_FAIL;
        }

        FatWater::ChemicalSpecies water = {"water", {{1.0, 0.0}}};
        FatWater::ChemicalSpecies fat = {"fat",
                                       {
                                               {0.07960f - 0.0510if, -3.8415},
                                               {0.64660f           , -3.4860},
                                               {0.09570f + 0.0140if, -2.7583},
                                               {-0.0047f - 0.0221if, -1.8762},
                                               {0.01600f - 0.0150if, -0.5047},
                                               {0.08490f - 0.0244if,  0.5260}
                                       }};

//        FatWater::ChemicalSpecies fat = {"fat",
//                                         {
//                                                 {0.048, 5.3-4.7},
//                                                 {0.039           , 4.31-4.7},
//                                                 {0.004, 2.76-4.7},
//                                                 {0.128, 2.1-4.7},
//                                                 {0.693, 1.3-4.7},
//                                                 {0.087, 0.9-4.7}
//                                         }};
        parameters.species = {water, fat};


        //try {
        //This should return the images
        auto output = FatWater::fatwater_separation(imagearr.data_, parameters, config);
        auto& wfimages = output.images;
        auto& field_map = output.field_map;
        auto& r2star_map = output.r2star_map;

        if (sample_time_us > 0){
            correct_frequency_shift(wfimages,parameters);
        }

        //Now let's make an image array to return the f/w images + any additional recon products
        //Right now, we will be using a copy of the original data
        //TODO: Remove this data copy
        //wfimages = imagearr.data_;

        uint16_t n_images = wfimages.get_size(4);
        uint16_t s_images = wfimages.get_size(5); //S-dimention is the image dimension
        uint16_t loc_images = wfimages.get_size(6);

        if (n_images != N || loc_images != LOC) {
            GERROR("Wrong number of N or LOC images received from fat water separation\n");
            m1->release();
            return GADGET_FAIL;
        }


        GadgetContainerMessage<IsmrmrdImageArray> *m2 = FatWaterImageArray(parameters, std::move(wfimages),
                                                                           m1->getObjectPtr()->headers_,
                                                                           m1->getObjectPtr()->meta_);

        if (this->next()->putq(m2) < 0) {
            m1->release();
            m2->release();
            return GADGET_FAIL;
        }


        if (save_field_map) {
            GadgetContainerMessage<ISMRMRD::ImageHeader> *m3 = MakeImageMessage(parameters, std::move(field_map),
                                                                                m1->getObjectPtr()->headers_, m1->getObjectPtr()->headers_[0].image_series_index+200);
            if (this->next()->putq(m3) < 0){
                m3->release();
                return GADGET_FAIL;
            }
        }

        if (save_r2star_map) {
            GadgetContainerMessage<ISMRMRD::ImageHeader> *m3 = MakeImageMessage(parameters, std::move(r2star_map),
                                                                                m1->getObjectPtr()->headers_, m1->getObjectPtr()->headers_[0].image_series_index+400);
            if (this->next()->putq(m3) < 0){
                m3->release();
                return GADGET_FAIL;
            }
        }

        //Pass the image down the chain
        if (this->next()->putq(m1) < 0) {
            return GADGET_FAIL;
        }

        return GADGET_OK;

    }

    GadgetContainerMessage <IsmrmrdImageArray> *
FatWaterGadget::FatWaterImageArray(const FatWater::Parameters &parameters, hoNDArray <std::complex<float>> &&wfimages,
                                   const hoNDArray <ISMRMRD::ImageHeader> &headers,
                                   const std::vector<ISMRMRD::MetaContainer> &metadata) const {


        uint16_t n_images = wfimages.get_size(4);
        uint16_t s_images = wfimages.get_size(5); //S-dimention is the image dimension
        uint16_t loc_images = wfimages.get_size(6);

        auto m2 = new GadgetContainerMessage<IsmrmrdImageArray>();
        m2->getObjectPtr()->data_ = std::move(wfimages);
        m2->getObjectPtr()->headers_.create(n_images, s_images, loc_images);
        for (uint16_t loc = 0; loc < loc_images; loc++) {
            for (uint16_t s = 0; s < s_images; s++) {
                for (uint16_t n = 0; n < n_images; n++) {
                    m2->getObjectPtr()->headers_[loc * s_images * n_images + s * n_images +
                                                 n] =headers[loc * s_images * n_images +
                                                                                   0 * n_images + n];

                    m2->getObjectPtr()->headers_[loc * s_images * n_images + s * n_images + n].image_series_index =
                            m2->getObjectPtr()->headers_[loc * s_images * n_images + s * n_images +
                                                         n].image_series_index + 100;

                    ISMRMRD::MetaContainer meta = metadata[loc * n_images * s_images + n];
                    //TODO: These sepcies and image type specifiers should come from the toolbox
                    if (s < parameters.species.size()) {
                        if (parameters.species[s].name == "water") {
                            meta.set(GADGETRON_DATA_ROLE, GADGETRON_IMAGE_WATER);
                        } else if (parameters.species[s].name == "fat") {
                            meta.set(GADGETRON_DATA_ROLE, GADGETRON_IMAGE_FAT);
                        }
                    } else {
                        //TODO: What to call these images
                    }
                    m2->getObjectPtr()->meta_.push_back(meta);
                }
            }
        }
        return m2;
    }

    GadgetContainerMessage<ISMRMRD::ImageHeader> *
    FatWaterGadget::MakeImageMessage(FatWater::Parameters parameters, hoNDArray<float> &&array,
                                     const hoNDArray<ISMRMRD::ImageHeader> &headers,
                                     uint16_t image_index) const {


        auto image = new GadgetContainerMessage<hoNDArray<float>>(std::move(array));

        auto header_message = new GadgetContainerMessage<ISMRMRD::ImageHeader>(headers[0]);
        header_message->getObjectPtr()->data_type = ISMRMRD::ISMRMRD_FLOAT;

        header_message->getObjectPtr()->image_series_index = image_index;

        header_message->cont(image);

        return header_message;
    }

    GADGET_FACTORY_DECLARE(FatWaterGadget)
}
