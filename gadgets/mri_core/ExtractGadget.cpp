/*
 * ExtractMagnitudeGadget.cpp
 *
 *  Created on: Nov 8, 2011
 *      Author: Michael S. Hansen
 */

#include <bitset>
#include <unordered_map>
#include <cpu/math/hoNDArray_math.h>
#include "GadgetIsmrmrdReadWrite.h"
#include "ExtractGadget.h"


namespace Gadgetron {

    namespace {
        using IMTYPE = ISMRMRD::ISMRMRD_ImageTypes;


        static const std::unordered_map<IMTYPE,std::function<float(std::complex<float>)>> extract_functions = {
                {IMTYPE::ISMRMRD_IMTYPE_MAGNITUDE, [](std::complex<float> v) { return std::abs(v); }},
                {IMTYPE::ISMRMRD_IMTYPE_REAL,[](std::complex<float> v) { return std::real(v); }},
                {IMTYPE::ISMRMRD_IMTYPE_IMAG,[](std::complex<float> v) { return std::imag(v); }},
                {IMTYPE::ISMRMRD_IMTYPE_PHASE,[](std::complex<float> v) { return std::arg(v); }}
        };

        static const std::unordered_map<IMTYPE,size_t> series_offset{
            {IMTYPE::ISMRMRD_IMTYPE_MAGNITUDE,0} ,
            {IMTYPE::ISMRMRD_IMTYPE_REAL,1000},
            {IMTYPE::ISMRMRD_IMTYPE_IMAG,2000},
            {IMTYPE::ISMRMRD_IMTYPE_PHASE,3000}};
    }

    ExtractGadget::ExtractGadget() {


    }

    ExtractGadget::~ExtractGadget() {

    }

    int ExtractGadget::process_config(ACE_Message_Block *mb) {


        if (int(extract_mask) > 0){
            const auto bitmask = std::bitset<4>(int(extract_mask));
            for (int imtype = IMTYPE::ISMRMRD_IMTYPE_MAGNITUDE; imtype <= IMTYPE::ISMRMRD_IMTYPE_PHASE; imtype++){
                if (bitmask[imtype-1]) image_types.push_back(IMTYPE(imtype));
            }
        } else {
            if (extract_magnitude) image_types.push_back(IMTYPE::ISMRMRD_IMTYPE_MAGNITUDE);
            if (extract_real) image_types.push_back(IMTYPE::ISMRMRD_IMTYPE_REAL);
            if (extract_imag) image_types.push_back(IMTYPE::ISMRMRD_IMTYPE_IMAG);
            if (extract_phase) image_types.push_back(IMTYPE::ISMRMRD_IMTYPE_PHASE);
        }

        if (image_types.empty()) throw std::runtime_error("ExtractGadget: No valid extract functions specified");

        return GADGET_OK;
    }

    int ExtractGadget::process(GadgetContainerMessage<ISMRMRD::ImageHeader> *m1,
                               GadgetContainerMessage<hoNDArray<std::complex<float>>> *m2) {

        for (IMTYPE image_type : image_types) {

            GadgetContainerMessage<ISMRMRD::ImageHeader> *cm1 =
                    new GadgetContainerMessage<ISMRMRD::ImageHeader>();

            //Copy the header
            *cm1->getObjectPtr() = *m1->getObjectPtr();

            auto cm2 = new GadgetContainerMessage<hoNDArray<float> >(m2->getObjectPtr()->get_dimensions());

            std::complex<float> *src = m2->getObjectPtr()->get_data_ptr();
            float *dst = cm2->getObjectPtr()->get_data_ptr();

            for (unsigned long i = 0; i < cm2->getObjectPtr()->get_number_of_elements(); i++) {
                dst[i] = extract_functions.at(image_type)(src[i]);
            }

            if (force_positive){
                *cm2->getObjectPtr() -=  min(cm2->getObjectPtr());
            }


            cm1->cont(cm2);
            cm1->getObjectPtr()->data_type = ISMRMRD::ISMRMRD_FLOAT;//GADGET_IMAGE_REAL_FLOAT;
            cm1->getObjectPtr()->image_type = image_type;
            cm1->getObjectPtr()->image_series_index += series_offset.at(image_type);

            if (this->next()->putq(cm1) == -1) {
                m1->release();
                GDEBUG("Unable to put extracted images on next gadgets queue");
                return GADGET_FAIL;
            }
        }


        m1->release(); //We have copied all the data in this case
        return GADGET_OK;
    }


    GADGET_FACTORY_DECLARE(ExtractGadget)

}
