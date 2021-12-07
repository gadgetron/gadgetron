#pragma once

#include "Gadget.h"
#include "hoNDArray.h"
#include "cuNDArray.h"
#include "gadgetron_hyper_export.h"

#include "CMRTOperator.h"
#include <ismrmrd/ismrmrd.h>
#include <complex>
#include <queue>

namespace Gadgetron{

  class EXPORTGADGETSHYPER CMRTGadget :
    public Gadget3< ISMRMRD::AcquisitionHeader, hoNDArray< std::complex<float> >, hoNDArray<float> >
  {
    
  public:
    using ReadoutMessage = GadgetContainerMessage<hoNDArray<std::complex<float>>>;
    using TrajectoryMessage = GadgetContainerMessage<hoNDArray<float>>;

    CMRTGadget(): num_frames(0) {
    }
    ~CMRTGadget() {}
    
  protected:
    
    GADGET_PROPERTY(golden_ratio,bool,"Use golden ratio trajectories",false);
    GADGET_PROPERTY(use_TV,bool,"Use total variation",false);
    GADGET_PROPERTY(projections_per_recon,int,"Number of projections per reconstruction",0);
    GADGET_PROPERTY(iterations,int,"Number of iterations",30);
    virtual int process_config(ACE_Message_Block* mb);

    virtual int process(GadgetContainerMessage< ISMRMRD::AcquisitionHeader > *m1,        // header
                        GadgetContainerMessage< hoNDArray< std::complex<float> > > *m2,  // data
                        GadgetContainerMessage< hoNDArray<float> > *m3 );                // traj/dcw

  protected:
        
    template<class T> GadgetContainerMessage< hoNDArray<T> >* 
      duplicate_array( GadgetContainerMessage< hoNDArray<T> > *array );        
    
    boost::shared_ptr< hoNDArray<float_complext> > 
      extract_samples_from_queue ( std::queue<std::unique_ptr<ReadoutMessage>> &queue );
    
    boost::shared_ptr< hoNDArray<float> >
      extract_trajectory_from_queue ( std::queue<std::unique_ptr<TrajectoryMessage>> &queue );
    
    void extract_trajectory_and_dcw_from_queue
      ( std::queue<std::unique_ptr<TrajectoryMessage>> &queue, boost::shared_ptr< hoNDArray<floatd2> > & traj, boost::shared_ptr< hoNDArray<float> > & dcw  );

    /***
     * Combines all stored frames and resets the frame buffer
     */
    boost::shared_ptr<cuNDArray<float_complext> > get_combined_frames();

  protected:

    std::vector<size_t> image_space_dimensions_3D_;
    unsigned int projections_per_recon_;


    std::queue<std::unique_ptr<ReadoutMessage>> frame_readout_queue_;
    std::queue<std::unique_ptr<TrajectoryMessage>> frame_traj_queue_;

    std::vector<size_t> dimensions_;
    std::vector<float> field_of_view_;
    size_t repetitions_;
    size_t samples_per_readout_;
    size_t num_coils_;
    size_t num_trajectory_dims_; // 2 for trajectories only, 3 for both trajectories + dcw

    std::vector<boost::shared_ptr<hoNDArray<float_complext> > > frames;
    boost::shared_ptr<hoNDArray<float> > dcw;
    boost::shared_ptr<hoNDArray<floatd2> > traj;
    unsigned int num_frames;
    unsigned int iterations_;
    bool golden_ratio_;
    bool use_TV_;
  };
}
