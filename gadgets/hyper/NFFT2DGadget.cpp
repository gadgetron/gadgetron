#include "NFFT2DGadget.h"
#include "cuNFFT.h"
#include "vector_td_utilities.h"
#include "hoNDArray_utils.h"
#include "hoNDArray_fileio.h"
#include "cuNDArray_utils.h"
#include "ismrmrd/xml.h"
#include "radial_utilities.h"
#include <cmath>

namespace Gadgetron{

  int NFFT2DGadget::process_config(ACE_Message_Block* mb)
  {

  	ISMRMRD::IsmrmrdHeader h;
  	ISMRMRD::deserialize(mb->rd_ptr(),h);


    if (h.encoding.size() != 1) {
      GDEBUG("This Gadget only supports one encoding space\n");
      return GADGET_FAIL;
    }

    ISMRMRD::EncodingSpace e_space = h.encoding[0].encodedSpace;
    ISMRMRD::EncodingSpace r_space = h.encoding[0].reconSpace;
    ISMRMRD::EncodingLimits e_limits = h.encoding[0].encodingLimits;


    GDEBUG("Matrix size: %d, %d\n", e_space.matrixSize.x, e_space.matrixSize.y, e_space.matrixSize.z);
    dimensions_.push_back(r_space.matrixSize.x);
    dimensions_.push_back(r_space.matrixSize.y);

    field_of_view_.push_back(e_space.fieldOfView_mm.x);
    field_of_view_.push_back(e_space.fieldOfView_mm.y);
    GDEBUG("FOV: %f, %f\n", r_space.fieldOfView_mm.x, r_space.fieldOfView_mm.y);

    repetitions_ = e_limits.repetition.is_present() ? e_limits.repetition.get().maximum + 1 : 1;
    GDEBUG("#Repetitions: %d\n", repetitions_);

    // Allocate readout and trajectory/dcw queues
    //

    frame_readout_queue_ = boost::shared_ptr< ACE_Message_Queue<ACE_MT_SYNCH> >(new ACE_Message_Queue<ACE_MT_SYNCH>());
    frame_traj_queue_ = boost::shared_ptr< ACE_Message_Queue<ACE_MT_SYNCH> >(new ACE_Message_Queue<ACE_MT_SYNCH>());
    
    size_t bsize = sizeof(GadgetContainerMessage< hoNDArray< std::complex<float> > >)*dimensions_[0]*10;
    
    frame_readout_queue_->high_water_mark(bsize);
    frame_readout_queue_->low_water_mark(bsize);
    frame_traj_queue_->high_water_mark(bsize);
    frame_traj_queue_->low_water_mark(bsize);

    return GADGET_OK;
  }

  int NFFT2DGadget::process(GadgetContainerMessage< ISMRMRD::AcquisitionHeader > *m1,        // header
                            GadgetContainerMessage< hoNDArray< std::complex<float> > > *m2,  // data
                            GadgetContainerMessage< hoNDArray<float> > *m3 )                 // traj/dcw
  {    
    // Throw away any noise samples if they have been allowed to pass this far down the chain...
    //
    
  	bool is_noise = m1->getObjectPtr()->isFlagSet(ISMRMRD::ISMRMRD_ACQ_IS_NOISE_MEASUREMENT);
    if (is_noise) { 
      m1->release();
      return GADGET_OK;
    }
    
    // First pass initialization
    //
    
    if (frame_readout_queue_->message_count() == 0 ) {      
      samples_per_readout_ = m1->getObjectPtr()->number_of_samples;
      num_coils_ = m1->getObjectPtr()->active_channels;      
      dimensions_.push_back(m1->getObjectPtr()->active_channels);
      dimensions_.push_back(repetitions_);
      num_trajectory_dims_ = m3->getObjectPtr()->get_size(0); // 2 for trajectories only, 3 for both trajectories + dcw
    }

    int samples = m1->getObjectPtr()->number_of_samples;
    int readout = m1->getObjectPtr()->idx.kspace_encode_step_1;
    int repetition = m1->getObjectPtr()->idx.kspace_encode_step_2;

    // Enqueue incoming readouts and trajectories
    //

    frame_readout_queue_->enqueue_tail(duplicate_array(m2));
    frame_traj_queue_->enqueue_tail(duplicate_array(m3));
    
    // If the last readout for a slice has arrived then perform a reconstruction
    //

    bool is_last_scan_in_repetition = 
    		m1->getObjectPtr()->isFlagSet(ISMRMRD::ISMRMRD_ACQ_LAST_IN_REPETITION);

    if (is_last_scan_in_repetition) {


      // Define the image header
      //

      GadgetContainerMessage<ISMRMRD::ImageHeader> *cm1 = 
        new GadgetContainerMessage<ISMRMRD::ImageHeader>();      
      
      GadgetContainerMessage< hoNDArray< std::complex<float> > > *cm2 = 
        new GadgetContainerMessage<hoNDArray< std::complex<float> > >();
      
      cm1->getObjectPtr()->flags = 0;
      cm1->cont(cm2);
    
      cm1->getObjectPtr()->matrix_size[0]     = dimensions_[0];
      cm1->getObjectPtr()->matrix_size[1]     = dimensions_[1];
      cm1->getObjectPtr()->matrix_size[2]     = 1;
      cm1->getObjectPtr()->field_of_view[0]   = field_of_view_[0];
      cm1->getObjectPtr()->field_of_view[1]   = field_of_view_[1];
      cm1->getObjectPtr()->channels           = num_coils_;
      cm1->getObjectPtr()->repetition         = m1->getObjectPtr()->idx.repetition;

      memcpy(cm1->getObjectPtr()->position,
             m1->getObjectPtr()->position,
             sizeof(float)*3);

      memcpy(cm1->getObjectPtr()->read_dir,
             m1->getObjectPtr()->read_dir,
             sizeof(float)*3);

      memcpy(cm1->getObjectPtr()->phase_dir,
             m1->getObjectPtr()->phase_dir,
             sizeof(float)*3);

      memcpy(cm1->getObjectPtr()->slice_dir,
             m1->getObjectPtr()->slice_dir,
             sizeof(float)*3);

      memcpy(cm1->getObjectPtr()->patient_table_position,
             m1->getObjectPtr()->patient_table_position, sizeof(float)*3);

      cm1->getObjectPtr()->data_type = ISMRMRD::ISMRMRD_CXFLOAT;
      cm1->getObjectPtr()->image_index = 0;
      cm1->getObjectPtr()->image_series_index = 0;

      //
      // Perform reconstruction of repetition
      //
      
      // Get samples for frame
      //

      cuNDArray<float_complext> samples( extract_samples_from_queue( frame_readout_queue_.get()).get() );

      // Get trajectories/dcw for frame
      //
      
      boost::shared_ptr< cuNDArray<floatd2> > traj(new cuNDArray<floatd2>);
      boost::shared_ptr<cuNDArray<float> > dcw(new cuNDArray<float>);

      extract_trajectory_and_dcw_from_queue( frame_traj_queue_.get(), traj.get(), dcw.get() );
      //traj = compute_radial_trajectory_golden_ratio_2d<float>(samples_per_readout_,dimensions_[1],1,0,GR_ORIGINAL);

      unsigned int num_profiles = samples.get_number_of_elements()/samples_per_readout_;
      dcw = compute_radial_dcw_golden_ratio_2d<float>(samples_per_readout_,num_profiles,1.0,1.0f/samples_per_readout_/num_profiles,0,GR_ORIGINAL);
      // Create output array
      //


      std::vector<size_t> img_dims(2);
      img_dims[0] = dimensions_[0];
      img_dims[1] = dimensions_[1];
      cm2->getObjectPtr()->create(img_dims);
      cuNDArray<float_complext> image(&img_dims);
      
      // Initialize plan
      //
      
      const float kernel_width = 5.5f;
      cuNFFT_impl<float,2> plan( from_std_vector<size_t,2>(img_dims), from_std_vector<size_t,2>(img_dims)<<1, kernel_width );
      plan.preprocess( traj.get(), NFFT_prep_mode::NC2C );
/*
      if( dcw->get_number_of_elements() == 0 ){
        std::vector<size_t> dcw_dims; dcw_dims.push_back(samples_per_readout_);
        hoNDArray<float> host_dcw( dcw_dims );
        for( int i=0; i<(int)dcw_dims[0]; i++ )
          host_dcw.get_data_ptr()[i]=abs(i-(int)dcw_dims[0]/2);
        host_dcw.get_data_ptr()[dcw_dims[0]/2] = 0.25f; // ad hoc value (we do not want a DC component of 0)        
        dcw = expand(&host_dcw, traj->get_number_of_elements()/samples_per_readout_);
      }
*/
      // Gridder
      //
      
      plan.compute( samples,image,
                    (dcw->get_number_of_elements()>0) ? dcw.get() : 0x0,
                    NFFT_comp_mode::BACKWARDS_NC2C );


      // Download to host
      //

      image.to_host( (hoNDArray<float_complext>*)cm2->getObjectPtr() );
      // Pass on data down the gadget chain
      //

      if (this->next()->putq(cm1) < 0) {
        return GADGET_FAIL;
      }
    }

    m1->release();
    return GADGET_OK;
  }
  
  template<class T> GadgetContainerMessage< hoNDArray<T> >*
  NFFT2DGadget::duplicate_array( GadgetContainerMessage< hoNDArray<T> > *array )
  {
    GadgetContainerMessage< hoNDArray<T> > *copy = new GadgetContainerMessage< hoNDArray<T> >();   
    *(copy->getObjectPtr()) = *(array->getObjectPtr());
    return copy;
  }
  
  boost::shared_ptr< hoNDArray<float_complext> > 
  NFFT2DGadget::extract_samples_from_queue ( ACE_Message_Queue<ACE_MT_SYNCH> *queue )                                             
  {    
    if(!queue) {
      GDEBUG("Illegal queue pointer, cannot extract samples\n");
      throw std::runtime_error("NFFT2DGadget::extract_samples_from_queue: illegal queue pointer");	
    }

    unsigned int readouts_buffered = queue->message_count();
    
    std::vector<size_t> dims;
    dims.push_back(samples_per_readout_*readouts_buffered);
    dims.push_back(num_coils_);
    
    boost::shared_ptr< hoNDArray<float_complext> > host_samples(new hoNDArray<float_complext>(dims));
    
    for (unsigned int p=0; p<readouts_buffered; p++) {
      
      ACE_Message_Block* mbq;
      if (queue->dequeue_head(mbq) < 0) {
        GDEBUG("Message dequeue failed\n");
        throw std::runtime_error("NFFT2DGadget::extract_samples_from_queue: dequeing failed");	
      }
      
      GadgetContainerMessage< hoNDArray< std::complex<float> > > *daq = AsContainerMessage<hoNDArray< std::complex<float> > >(mbq);
	
      if (!daq) {
        GDEBUG("Unable to interpret data on message queue\n");
        throw std::runtime_error("NFFT2DGadget::extract_samples_from_queue: failed to interpret data");	
      }
	
      for (unsigned int c = 0; c < num_coils_; c++) {
	
        float_complext *data_ptr = host_samples->get_data_ptr();
        data_ptr += c*samples_per_readout_*readouts_buffered+p*samples_per_readout_;
	    
        std::complex<float> *r_ptr = daq->getObjectPtr()->get_data_ptr();
        r_ptr += c*daq->getObjectPtr()->get_size(0);
	  
        memcpy(data_ptr, r_ptr, samples_per_readout_*sizeof(float_complext));
      }

      mbq->release();
    }
    
    return host_samples;
  }

  boost::shared_ptr< hoNDArray<float> > 
  NFFT2DGadget::extract_trajectory_from_queue ( ACE_Message_Queue<ACE_MT_SYNCH> *queue )
  {    
    if(!queue) {
      GDEBUG("Illegal queue pointer, cannot extract trajectory\n");
      throw std::runtime_error("NFFT2DGadget::extract_trajectory_from_queue: illegal queue pointer");	
    }

    unsigned int readouts_buffered = queue->message_count();
    
    std::vector<size_t> dims;
    dims.push_back(num_trajectory_dims_); // 2 for trajectories only, 3 for both trajectories + dcw
    dims.push_back(samples_per_readout_);
    dims.push_back(readouts_buffered);
    
    boost::shared_ptr< hoNDArray<float> > host_traj(new hoNDArray<float>(dims));
    
    for (unsigned int p=0; p<readouts_buffered; p++) {      
      ACE_Message_Block* mbq;
      if (queue->dequeue_head(mbq) < 0) {
        GDEBUG("Message dequeue failed\n");
        throw std::runtime_error("NFFT2DGadget::extract_trajectory_from_queue: dequeing failed");	
      }
      
      GadgetContainerMessage< hoNDArray<float> > *daq = AsContainerMessage<hoNDArray<float> >(mbq);
	
      if (!daq) {
        GDEBUG("Unable to interpret data on message queue\n");
        throw std::runtime_error("NFFT2DGadget::extract_trajectory_from_queue: failed to interpret data");	
      }

      float *data_ptr = host_traj->get_data_ptr();
      data_ptr += num_trajectory_dims_*samples_per_readout_*p;
      
      float *r_ptr = daq->getObjectPtr()->get_data_ptr();
      
      memcpy(data_ptr, r_ptr, num_trajectory_dims_*samples_per_readout_*sizeof(float));
      
      mbq->release();
    }
    
    return host_traj;
  }

  void NFFT2DGadget::extract_trajectory_and_dcw_from_queue
  ( ACE_Message_Queue<ACE_MT_SYNCH> *queue, cuNDArray<floatd2> *traj, cuNDArray<float> *dcw )
  {
    // Extract trajectory and (if present) density compensation weights.
    // They are stored as a float array of dimensions: {2,3} x #samples_per_readout x #readouts.
    // We need
    // - a floatd2 trajectory array 
    // - a float dcw array 
    //
        
    if( num_trajectory_dims_ == 2 ){
      boost::shared_ptr< hoNDArray<float> > host_traj = extract_trajectory_from_queue( queue );
      std::vector<size_t> dims_1d; dims_1d.push_back(host_traj->get_size(1)*host_traj->get_size(2));
      hoNDArray<floatd2> host_traj2(dims_1d,(floatd2*)host_traj->get_data_ptr());
      *traj = cuNDArray<floatd2>(host_traj2);

    }
    else{

      boost::shared_ptr< hoNDArray<float> > host_traj_dcw = extract_trajectory_from_queue( queue );

      std::vector<size_t> order;
      order.push_back(1); order.push_back(2); order.push_back(0);
      
      auto host_traj_dcw_shifted = permute( *host_traj_dcw, order );
      
      std::vector<size_t> dims_1d;
      dims_1d.push_back(host_traj_dcw_shifted.get_size(0)*host_traj_dcw_shifted.get_size(1));
      
      hoNDArray<float> tmp(dims_1d, host_traj_dcw_shifted.get_data_ptr()+2*dims_1d[0]);
      *dcw = tmp;
      
      std::vector<size_t> dims_2d = dims_1d; dims_2d.push_back(2);
      order.clear(); order.push_back(1); order.push_back(0);
      
      tmp.create(dims_2d, host_traj_dcw_shifted.get_data_ptr());
      auto _traj = permute( tmp, order );
      hoNDArray<floatd2> tmp2(dims_1d,(floatd2*)_traj.get_data_ptr());
      
      *traj = cuNDArray<floatd2>(tmp2);
    }

    std::vector<size_t >dims_2d;
    dims_2d.push_back(traj->get_number_of_elements());
    dims_2d.push_back(1); // Number of frames

    traj->reshape(&dims_2d);
    if( num_trajectory_dims_ == 3 ) dcw->reshape(&dims_2d);
  }

  GADGET_FACTORY_DECLARE(NFFT2DGadget)
}
