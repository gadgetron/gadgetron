#ifndef GPUCGGADGET_H
#define GPUCGGADGET_H
#pragma once

#include <complex>

#include "gadgetron_export.h"
#include "Gadget.h"
#include "GadgetMRIHeaders.h"
#include "cuCGSolver.h"
#include "cuNonCartesianSenseOperator.h"
#include "cuCGPrecondWeights.h"
#include "NFFT.h"
#include "cuSenseRHSBuffer.h"
#include "cuImageOperator.h"

class EXPORTGADGETSCGSENSE GPUCGGadget : public Gadget2< GadgetMessageAcquisition, hoNDArray< std::complex<float> > >
{

public:
  GPUCGGadget();
  virtual ~GPUCGGadget();

protected:

  virtual int process( GadgetContainerMessage< GadgetMessageAcquisition >* m1, GadgetContainerMessage< hoNDArray< std::complex<float> > > * m2 );
  virtual int process_config( ACE_Message_Block* mb );

  virtual boost::shared_ptr< cuNDArray<floatd2::Type> > calculate_trajectory() = 0;
  virtual boost::shared_ptr< cuNDArray<float> > calculate_density_compensation() = 0;

  virtual int copy_samples_for_profile( float* host_base_ptr, std::complex<float>* data_base_ptr, int profile_no, int channel_no );
  virtual int configure_channels();
  virtual boost::shared_ptr< cuNDArray<float_complext::Type> > upload_samples();

  ACE_Message_Queue<ACE_MT_SYNCH> buffer_;
  int slice_no_;
  int profiles_per_frame_;
  int shared_profiles_;
  int channels_;
  int samples_per_profile_;
  int device_number_;
  uintd2::Type matrix_size_;
  uintd2::Type matrix_size_os_;
  unsigned int number_of_iterations_;
  double cg_limit_;
  double oversampling_;
  double kernel_width_;
  double kappa_;

  int current_profile_offset_;
  int allocated_samples_;
  float *data_host_ptr_;

  bool is_configured_;

  // Define conjugate gradient solver
  cuCGSolver<float, float_complext::Type> cg_;

  // Define non-Cartesian Sense Encofing operator
  boost::shared_ptr< cuNonCartesianSenseOperator<float,2> > E_;

  // Define preconditioner
  boost::shared_ptr< cuCGPrecondWeights<float_complext::Type> > D_;

  // Define regularization image operator
  boost::shared_ptr< cuImageOperator<float,float_complext::Type> > R_;

  // Define rhs operator (for regularization)
  boost::shared_ptr< cuSenseRHSBuffer<float,2> > rhs_buffer_;

  // Density compensation weights
  bool dcw_computed_;
};

#endif //GPUCGGADGET
