#pragma once

#include "Gadget.h"
#include "hoNDArray.h"
#include "gadgetron_mricore_export.h"
#include "GadgetronTimer.h"

#include <ismrmrd/ismrmrd.h>
#include <ismrmrd/xml.h>
#include <complex>

namespace Gadgetron {

  class EXPORTGADGETSMRICORE NoiseAdjustGadget :
    public Gadget2<ISMRMRD::AcquisitionHeader,hoNDArray< std::complex<float> > >
    {
    public:
      GADGET_DECLARE(NoiseAdjustGadget);

      typedef Gadget2<ISMRMRD::AcquisitionHeader,hoNDArray< std::complex<float> > > BaseClass;

      NoiseAdjustGadget();
      virtual ~NoiseAdjustGadget();

      virtual int close(unsigned long flags);

    protected:
      bool noise_decorrelation_calculated_;
      hoNDArray< std::complex<float> > noise_covariance_matrixf_;
      hoNDArray< std::complex<float> > noise_prewhitener_matrixf_;
      hoNDArray< std::complex<float> > noise_covariance_matrixf_once_;
      std::vector<unsigned int> scale_only_channels_;

      unsigned long long number_of_noise_samples_;
      unsigned long long number_of_noise_samples_per_acquisition_;
      float noise_dwell_time_us_;
      float noise_dwell_time_us_preset_;
      float acquisition_dwell_time_us_;
      float noise_bw_scale_factor_;
      float receiver_noise_bandwidth_;
      bool noiseCovarianceLoaded_;
      bool perform_noise_adjust_;

      std::string noise_dependency_folder_;
      std::string noise_dependency_prefix_;
      std::string measurement_id_;
      std::string measurement_id_of_noise_dependency_;
      std::string full_name_stored_noise_dependency_;

      virtual int process_config(ACE_Message_Block* mb);
      virtual int process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1,
			  GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2);

      std::string generateNoiseDependencyFilename(const std::string& measurement_id);
      std::string generateMeasurementIdOfNoiseDependency(const std::string& noise_id);

      bool loadNoiseCovariance();
      bool saveNoiseCovariance();
      void computeNoisePrewhitener();

      //We will store/load a copy of the noise scans XML header to enable us to check which coil layout, etc.
      ISMRMRD::IsmrmrdHeader current_ismrmrd_header_;
      ISMRMRD::IsmrmrdHeader noise_ismrmrd_header_;

    };
}
