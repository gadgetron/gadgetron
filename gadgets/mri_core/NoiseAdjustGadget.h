#pragma once

#include "Gadget.h"
#include "hoNDArray.h"
#include "gadgetron_mricore_export.h"
#include "GadgetronTimer.h"

#include <ismrmrd/ismrmrd.h>
#include <complex>

namespace Gadgetron {

  class EXPORTGADGETSMRICORE NoiseAdjustGadget :
    public Gadget2<ISMRMRD::AcquisitionHeader,hoNDArray< std::complex<float> > >
    {
    public:
      GADGET_DECLARE(NoiseAdjustGadget);

      typedef Gadget2<ISMRMRD::AcquisitionHeader,hoNDArray< std::complex<float> > > BaseClass;

      typedef std::complex<float> ValueType;
      typedef std::complex<double> PerwhitenerValueType;

      NoiseAdjustGadget();
      virtual ~NoiseAdjustGadget();

      virtual int close(unsigned long flags);

    protected:
      bool noise_decorrelation_calculated_;
      hoNDArray< ValueType > noise_covariance_matrixf_;
      hoNDArray< ValueType > noise_covariance_matrixf_once_;

      hoNDArray< ValueType > data_prewhitened_;

      hoNDArray< ValueType > readout_;

      unsigned long long number_of_noise_samples_;
      unsigned long long number_of_noise_samples_per_acquisition_;
      float noise_dwell_time_us_;
      float acquisition_dwell_time_us_;
      float noise_bw_scale_factor_;
      float receiver_noise_bandwidth_;
      bool is_configured_;
      bool computed_in_close_;

      std::string noise_dependency_folder_;

      std::string noise_dependency_prefix_;

      std::string patient_id_;
      std::string study_id_;
      std::string measurement_id_;
      std::string measurement_id_of_noise_dependency_;

      std::string full_name_stored_noise_dependency_;

      float noise_dwell_time_us_preset_;

      bool perform_noise_adjust_;

      virtual int process_config(ACE_Message_Block* mb);
      virtual int process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1,
			  GadgetContainerMessage< hoNDArray< ValueType > >* m2);

      std::string generateNoiseDependencyFilename(const std::string& measurement_id);
      std::string generateMeasurementIdOfNoiseDependency(const std::string& noise_id);

      bool use_stored_noise_prewhitener_;
      bool loadNoisePrewhitener(float& noise_dwell_time_us, hoNDArray< ValueType >& noise_covariance_matrixf);
      bool saveNoisePrewhitener(const std::string& full_name_stored_noise_dependency, float& noise_dwell_time_us, hoNDArray< ValueType >& noise_covariance_matrixf);

      void computeNoisePrewhitener(bool savePrewhitener=true);

      Gadgetron::GadgetronTimer gt_timer_;
      bool performTiming_;
    };
}
