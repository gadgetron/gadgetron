#pragma once

#include "Gadget.h"
#include "hoNDArray.h"
#include "gadgetron_mricore_export.h"

#include <ismrmrd.h>
#include <complex>

namespace Gadgetron {

    class EXPORTGADGETSMRICORE NoiseAdjustGadget :
        public Gadget2<ISMRMRD::AcquisitionHeader,hoNDArray< std::complex<float> > >
    {
    public:

        typedef std::complex<float> ValueType;
        typedef std::complex<double> PerwhitenerValueType;

        NoiseAdjustGadget();

    protected:
        bool noise_decorrelation_calculated_;
        hoNDArray< PerwhitenerValueType > noise_covariance_matrix_;
        hoNDArray< ValueType > noise_covariance_matrixf_;
        unsigned long int number_of_noise_samples_;
        float noise_dwell_time_us_;
        float acquisition_dwell_time_us_;
        float noise_bw_scale_factor_;
        float receiver_noise_bandwidth_;
        bool is_configured_;
        hoNDArray< ValueType > prewhitened_buf_;

        virtual int process_config(ACE_Message_Block* mb);
        virtual int process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1,
            GadgetContainerMessage< hoNDArray< ValueType > >* m2);
    };
}
