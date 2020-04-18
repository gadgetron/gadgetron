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
      GADGET_PROPERTY(noise_dependency_prefix, std::string, "Prefix of noise depencency file", "GadgetronNoiseCovarianceMatrix");
      GADGET_PROPERTY(perform_noise_adjust, bool, "Whether to actually perform the noise adjust", true);
      GADGET_PROPERTY(pass_nonconformant_data, bool, "Whether to pass data that does not conform", false);
      GADGET_PROPERTY(noise_dwell_time_us_preset, float, "Preset dwell time for noise measurement", 0.0);
      GADGET_PROPERTY(scale_only_channels_by_name, std::string, "List of named channels that should only be scaled", "");
      GADGET_PROPERTY(clean_items_older_than_thres, double, "Clean items older than this number of hours", 24.0);

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
      bool pass_nonconformant_data_;
      bool saved_;

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
      std::vector<size_t> coil_order_of_data_in_noise_;

      void print_coil_label(const std::vector<ISMRMRD::CoilLabel>& coils);

      // compare coil labels of noise and data
      // if number of channels are different, return false and order.size()==0
      // if any channels in noise cannot be found in data, return false and order.size()==0
      // if all channels in noise exist in data, but order is incorrect, return false, but  and order.size()==CHA
      // if all channels in nosie match channels in data, return true
      // order gives the matching order for src and dst coils
      // e.g. [2 1 0 3] means coil 0 of src matches coil 2 of dst etc.
      bool compare_coil_label(const std::vector<ISMRMRD::CoilLabel>& src_coils, const std::vector<ISMRMRD::CoilLabel>& dst_coils, std::vector<size_t>& order);
    };
}
