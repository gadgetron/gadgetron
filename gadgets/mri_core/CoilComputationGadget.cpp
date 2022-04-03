#include "CoilComputationGadget.h"
#include "mri_core_coil_map_estimation.h"

namespace Gadgetron {

    Core::Image<std::complex<float>> CoilComputationGadget::process_function(Core::Image<std::complex<float>> image) const {
        
		auto &header = std::get<ISMRMRD::ImageHeader>(image);
        auto &data = std::get<hoNDArray<std::complex<float>>>(image);

        hoNDArray<std::complex<float>> csm;

        if (data.get_number_of_dimensions() == 4)
        {
            coil_map_Inati<std::complex<float>>(data, csm, ks_, kz_, power_);
        }
        else
        {
            GERROR_STREAM("CoilComputationGadget, no 4D data was passed ... ");
        }

        return Core::Image<std::complex<float>>(
                std::get<ISMRMRD::ImageHeader>(image),
                csm,
                std::get<Core::optional<ISMRMRD::MetaContainer>>(image)
        );

    }

    GADGETRON_GADGET_EXPORT(CoilComputationGadget);
}


