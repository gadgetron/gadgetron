#ifndef MaxwellCorrectionGadget_H
#define MaxwellCorrectionGadget_H

#include "Gadget.h"
#include "GadgetMRIHeaders.h"
#include "hoNDArray.h"
#include "gadgetron_mricore_export.h"

#include <ismrmrd/ismrmrd.h>
#include <complex>

namespace Gadgetron{  

    class EXPORTGADGETSMRICORE MaxwellCorrectionGadget :
        public Gadget2< ISMRMRD::ImageHeader, hoNDArray< std::complex<float> > >
    {

    public:
        GADGET_DECLARE(MaxwellCorrectionGadget);
        MaxwellCorrectionGadget();
        virtual ~MaxwellCorrectionGadget();
        void patient_to_physical_coordinate(std::vector<float> &norm_vec, std::string patient_position);
        void find_flow_dir(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1);

    protected:
        virtual int process_config(ACE_Message_Block* mb);
        virtual int process(GadgetContainerMessage< ISMRMRD::ImageHeader >* m1,
            GadgetContainerMessage< hoNDArray< std::complex<float> > > * m2);
	
    private:
	std::vector<double> maxwell_coefficients_;
	std::vector<float> RO_dir_Physical_;
	std::vector<float> PE_dir_Physical_;
	std::vector<float> SLC_dir_Physical_;
	std::vector<float> SLC_position_Physical_;

    bool FlipPhaseDirection_ = false;
    int FlowDirection_ = 4; //flow encoding direction: 4 through plane, 2 RO direction, 1 PE direction


	std::string patient_position_;

	bool maxwell_coefficients_present_;
    };
}

#endif //MaxwellCorrectionGadget_H
