#ifndef CUSBGADGET_H
#define CUSBGADGET_H
#pragma once

#include "gadgetronsbsense_export.h"
#include "cuNFFT.h"
#include "cuSbcCgSolver.h"
#include "Gadget.h"
#include "GadgetMRIHeaders.h"
#include "cuNonCartesianSenseOperator.h"
#include "cuPartialDerivativeOperator.h"
#include "cuCgPreconditioner.h"
#include "cuSenseRHSBuffer.h"
#include "cuImageOperator.h"
#include "ismrmrd.h"

#include <ace/Synch.h>
#include <ace/Mutex.h>
#include <complex>

namespace Gadgetron {

  class EXPORTGADGETSSBSENSE cuSbGadget : public Gadget2< ISMRMRD::AcquisitionHeader, hoNDArray< std::complex<float> > >
  {

  public:
    cuSbGadget();
    virtual ~cuSbGadget();

  protected:
    virtual int parameter_changed(std::string name, std::string new_value, std::string old_value);

    bool position_equal(float* position) {
      for (unsigned int i = 0; i < 3; i++) {
	if (position_[i] != position[i]) return false;
      }
      return true;
    }

    bool read_dir_equal(float* cosines) {
      for (unsigned int i = 0; i < 3; i++) {
	if (read_dir_[i] != cosines[i]) return false;
      }
      return true;
    }

    bool phase_dir_equal(float* cosines) {
      for (unsigned int i = 0; i < 3; i++) {
	if (phase_dir_[i] != cosines[i]) return false;
      }
      return true;
    }

    bool slice_dir_equal(float* cosines) {
      for (unsigned int i = 0; i < 3; i++) {
	if (slice_dir_[i] != cosines[i]) return false;
      }
      return true;
    }

    virtual int process( GadgetContainerMessage< ISMRMRD::AcquisitionHeader >* m1, GadgetContainerMessage< hoNDArray< std::complex<float> > > * m2 );
    virtual int process_config( ACE_Message_Block* mb );

    virtual boost::shared_ptr< cuNDArray<floatd2> > calculate_trajectory() = 0;
    virtual boost::shared_ptr< cuNDArray<float> > calculate_density_compensation() = 0;

    virtual int copy_samples_for_profile( float* host_base_ptr, std::complex<float>* data_base_ptr, int profile_no, int channel_no );
    virtual int configure_channels();
    virtual boost::shared_ptr< cuNDArray<float_complext> > upload_samples();

    ACE_Message_Queue<ACE_MT_SYNCH> buffer_;
    int slice_no_;
    int profiles_per_frame_;
    int shared_profiles_;
    int channels_;
    int samples_per_profile_;
    int device_number_;
    uintd2 matrix_size_;
    uintd2 matrix_size_os_;
    unsigned int number_of_cg_iterations_;
    unsigned int number_of_sb_iterations_;
    double cg_limit_;
    double oversampling_;
    double kernel_width_;
    double mu_;
    double lambda_;
    double alpha_;

    float position_[3];
    float read_dir_[3];
    float phase_dir_[3];
    float slice_dir_[3];

    int current_profile_offset_;
    int allocated_samples_;
    float *data_host_ptr_;

    bool is_configured_;

    // Define constraint Split Bregman solver
    cuSbcCgSolver<float_complext> sb_;

    // Define non-Cartesian Sense Encofing operator
    boost::shared_ptr< cuNonCartesianSenseOperator<float,2> > E_;

    // Define preconditioner
    boost::shared_ptr< cuCgPreconditioner<float_complext> > D_;

    // Define regularization operators
    boost::shared_ptr< cuPartialDerivativeOperator<float_complext,2> > Rx1_;
    boost::shared_ptr< cuPartialDerivativeOperator<float_complext,2> > Rx2_;
    boost::shared_ptr< cuPartialDerivativeOperator<float_complext,2> > Ry1_;
    boost::shared_ptr< cuPartialDerivativeOperator<float_complext,2> > Ry2_;
    //boost::shared_ptr< cuPartialDerivativeOperator<float_complext,3> > Rz1_;
    //boost::shared_ptr< cuPartialDerivativeOperator<float_complext,3> > Rz2_;

    // Define rhs operator (for regularization)
    boost::shared_ptr< cuSenseRHSBuffer<float,2> > rhs_buffer_;

    // Density compensation weights
    bool dcw_computed_;

    int image_series_;
    int image_counter_;

    ACE_Thread_Mutex mutex_;
  };
}

#endif //CUSBGADGET
