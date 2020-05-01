//
// Created by dchansen on 3/9/20.
//
#include "demons_registration.h"

#include "PureGadget.h"
#include "mri_core_data.h"
#include "t1fit.h"
#include "hoNDArray_math.h"

namespace Gadgetron {

    class T1Gadget : Core::PureGadget<Core::Image<float>, IsmrmrdImageArray> {

    public:
        T1Gadget(const Core::Context& context, const Core::GadgetProperties& properties)
            : Core::PureGadget<Core::Image<float>, IsmrmrdImageArray>(context, properties), TIs{*(context.header.sequenceParameters->TI)} {


        }

    private:
        Core::Image<float> process_function(IsmrmrdImageArray images) const final {



            auto vector_field = T1::t1_registration(images.data_,TIs);

            auto moco_images = T1::deform_groups(images.data_,vector_field);

            auto phase_corrected = T1::phase_correct(moco_images,TIs);

            auto [A,B,T1star] = T1::fit_T1_3param(phase_corrected,TIs);

            B /= A;
            B -= 1;

            auto T1 = T1star;
            T1 *= B;
            auto header = images.headers_[0];
            header.image_series_index = ISMRMRD::ISMRMRD_IMTYPE_REAL;
            return Core::Image<float>{header,T1,Core::none};


        }


        std::vector<float> exctract_MOLLI_TI(const hoNDArray<ISMRMRD::AcquisitionHeader>& acq_headers){

            auto dims = acq_headers.dimensions(); //Do we know what these dimensions are?
            std::cout << "MOLLI HEADER DIMENSIONS REMOVE THIS LINE";
            for (auto dim : dims)
                std::cout << ' ' << dim;
            std::cout << std::endl;

            return {};

        }


        std::vector<float> TIs;
    };
}