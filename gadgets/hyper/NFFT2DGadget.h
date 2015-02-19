#pragma once

#include "Gadget.h"
#include "hoNDArray.h"
#include "cuNDArray.h"
#include "gadgetron_hyper_export.h"

#include <ismrmrd/ismrmrd.h>
#include <complex>

namespace Gadgetron{

  class EXPORTGADGETSHYPER NFFT2DGadget :
    public Gadget3< ISMRMRD::AcquisitionHeader, hoNDArray< std::complex<float> >, hoNDArray<float> >
  {
    
  public:
    //GADGET_DECLARE(NFFT2DGadget);
    
    NFFT2DGadget() {}
    ~NFFT2DGadget() {}
    
  protected:
    
    virtual int process_config(ACE_Message_Block* mb);

    virtual int process(GadgetContainerMessage< ISMRMRD::AcquisitionHeader > *m1,        // header
                        GadgetContainerMessage< hoNDArray< std::complex<float> > > *m2,  // data
                        GadgetContainerMessage< hoNDArray<float> > *m3 );                // traj/dcw

  protected:
        
    template<class T> GadgetContainerMessage< hoNDArray<T> >* 
      duplicate_array( GadgetContainerMessage< hoNDArray<T> > *array );        
    
    boost::shared_ptr< hoNDArray<float_complext> > 
      extract_samples_from_queue ( ACE_Message_Queue<ACE_MT_SYNCH> *queue );
    
    boost::shared_ptr< hoNDArray<float> > 
      extract_trajectory_from_queue ( ACE_Message_Queue<ACE_MT_SYNCH> *queue );
    
    void extract_trajectory_and_dcw_from_queue
      ( ACE_Message_Queue<ACE_MT_SYNCH> *queue, cuNDArray<floatd2> *traj, cuNDArray<float> *dcw );

  protected:
    
    boost::shared_ptr< ACE_Message_Queue<ACE_MT_SYNCH> > frame_readout_queue_;
    boost::shared_ptr< ACE_Message_Queue<ACE_MT_SYNCH> > frame_traj_queue_;
    std::vector<size_t> dimensions_;
    std::vector<float> field_of_view_;
    size_t repetitions_;
    size_t samples_per_readout_;
    size_t num_coils_;
    size_t num_trajectory_dims_; // 2 for trajectories only, 3 for both trajectories + dcw
  };
}
