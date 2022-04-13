#pragma once

#include "Gadget.h"
#include "hoNDArray.h"
#include "gadgetron_plplot_export.h"
#include "GadgetronTimer.h"

#include <ismrmrd/ismrmrd.h>
#include <ismrmrd/meta.h>
#include <ismrmrd/xml.h>
#include <complex>

namespace Gadgetron {

    class EXPORTPLPLOTGADGET NoiseCovariancePlottingGadget : public Gadget2<ISMRMRD::ImageHeader, hoNDArray< std::complex<float> > >
    {
    public:
        GADGET_DECLARE(NoiseCovariancePlottingGadget);

        typedef Gadget2<ISMRMRD::ImageHeader, hoNDArray< std::complex<float> > > BaseClass;

        NoiseCovariancePlottingGadget();
        virtual ~NoiseCovariancePlottingGadget();

        virtual int close(unsigned long flags);

    protected:

        GADGET_PROPERTY(noise_dependency_prefix, std::string, "Prefix of noise dependency file", "GadgetronNoiseCovarianceMatrix");
        GADGET_PROPERTY(xlabel, std::string, "Label for x axis", "Channel");
        GADGET_PROPERTY(ylabel, std::string, "Label for y axis", "Noise Standard Deviation");
        GADGET_PROPERTY(title, std::string, "Label for y axis", "Gadgetron, Noise STD Plot");
        GADGET_PROPERTY(xsize, int, "Noise variane plot, x size", 2048);
        GADGET_PROPERTY(ysize, int, "Noise variane plot, y size", 2048);
        GADGET_PROPERTY(series_num, int, "Noise plot image sereis number", 10000);

        bool noise_decorrelation_calculated_;
        hoNDArray< std::complex<float> > noise_covariance_matrixf_;

        float noise_dwell_time_us_;
        bool noiseCovarianceLoaded_;

        std::string noise_dependency_folder_;
        std::string noise_dependency_prefix_;
        std::string measurement_id_;
        std::string measurement_id_of_noise_dependency_;
        std::string full_name_stored_noise_dependency_;

        virtual int process_config(ACE_Message_Block* mb);

        virtual int process(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1,
            GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2);

        std::string generateNoiseDependencyFilename(const std::string& measurement_id);
        std::string generateMeasurementIdOfNoiseDependency(const std::string& noise_id);

        bool loadNoiseCovariance();
        void computeNoisePrewhitener();

        ISMRMRD::IsmrmrdHeader current_ismrmrd_header_;
        ISMRMRD::IsmrmrdHeader noise_ismrmrd_header_;

        ISMRMRD::ImageHeader curr_image_header_;
    };
}
