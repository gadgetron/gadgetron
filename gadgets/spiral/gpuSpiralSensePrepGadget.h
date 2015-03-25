#ifndef gpuSpiralSensePrepGadget_H
#define gpuSpiralSensePrepGadget_H
#pragma once

#include "gadgetron_spiral_export.h"
#include "Gadget.h"
#include "GadgetMRIHeaders.h"
#include "cuCgSolver.h"
#include "cuNonCartesianSenseOperator.h"
#include "cuCgPreconditioner.h"
#include "cuNFFT.h"
#include "hoNDArray.h"
#include "vector_td.h"
#include "cuNFFT.h"

#include <ismrmrd/ismrmrd.h>
#include <complex>
#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>

namespace Gadgetron{

  class EXPORTGADGETS_SPIRAL gpuSpiralSensePrepGadget :
    public Gadget2< ISMRMRD::AcquisitionHeader, hoNDArray< std::complex<float> > >
  {

  public:
    GADGET_DECLARE(gpuSpiralSensePrepGadget);

    gpuSpiralSensePrepGadget();
    virtual ~gpuSpiralSensePrepGadget();

  protected:
    GADGET_PROPERTY(deviceno, int, "GPU device number", 0);
    GADGET_PROPERTY(propagate_csm_from_set, int, "Which set to use for CSM", -1);
    GADGET_PROPERTY(buffer_using_solver, bool, "Use solver for buffer", false);
    GADGET_PROPERTY(use_multiframe_grouping, bool, "Use multiframe grouping", false);
    GADGET_PROPERTY(buffer_convolution_kernel_width, float, "Convolution kernel width for buffer", 5.5);
    GADGET_PROPERTY(buffer_convolution_oversampling_factor, float, "Oversampling used in buffer convolution", 1.25);
    GADGET_PROPERTY(reconstruction_os_factor_x, float, "Oversampling for reconstruction in x-direction", 1.0);
    GADGET_PROPERTY(reconstruction_os_factor_y, float, "Oversampling for reconstruction in y-direction", 1.0);

    virtual int process_config(ACE_Message_Block* mb);
    
    virtual int process(GadgetContainerMessage< ISMRMRD::AcquisitionHeader >* m1,
			GadgetContainerMessage< hoNDArray< std::complex<float> > > * m2);
    
    virtual GadgetContainerMessage<ISMRMRD::AcquisitionHeader>*
      duplicate_profile( GadgetContainerMessage<ISMRMRD::AcquisitionHeader> *profile );
    
  private:
    int samples_to_skip_start_;
    int samples_to_skip_end_;
    int samples_per_interleave_;
    int interleaves_;
    int slices_;
    int sets_;
    boost::shared_array<long> image_counter_;
    int device_number_;

    long    Tsamp_ns_;
    long    Nints_;
    boost::shared_array<long> interleaves_counter_singleframe_;
    boost::shared_array<long> interleaves_counter_multiframe_;
    long    acceleration_factor_;
    double  gmax_;
    double  smax_;
    double  krmax_;
    double  fov_;

    bool prepared_;
    bool use_multiframe_grouping_;
    bool buffer_using_solver_;

    int propagate_csm_from_set_;

    float kernel_width_;
    float oversampling_factor_;

    boost::shared_ptr< hoNDArray<floatd2> > host_traj_;
    boost::shared_ptr< hoNDArray<float> > host_weights_;
    
    boost::shared_array< hoNDArray<float_complext> > host_data_buffer_;
    boost::shared_ptr< cuNDArray<float> > dcw_buffer_;

    std::vector<size_t> fov_vec_;
    std::vector<size_t> image_dimensions_recon_;
    uint64d2 image_dimensions_recon_os_;

    cuNFFT_plan<float,2> nfft_plan_;
    cuCgSolver<float_complext> cg_;
    boost::shared_ptr< cuNDArray<float_complext> > csm_;
    boost::shared_ptr< cuNonCartesianSenseOperator<float,2> > E_;
    boost::shared_ptr< cuCgPreconditioner<float_complext> > D_;

    boost::shared_array< ACE_Message_Queue<ACE_MT_SYNCH> > buffer_;
    boost::shared_array< ACE_Message_Queue<ACE_MT_SYNCH> > image_headers_queue_;
  };
}
#endif //gpuSpiralSensePrepGadget_H
