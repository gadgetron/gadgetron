#include "gpuCgKtSenseGadget.h"
#include "cuNDArray_operators.h"
#include "cuNDArray_elemwise.h"
#include "cuNDArray_blas.h"
#include "cuNDArray_utils.h"
#include "cuNDArray_reductions.h"
#include "cuNDFFT.h"
#include "GadgetMRIHeaders.h"
#include "b1_map.h"
#include "GPUTimer.h"
#include "vector_td_utilities.h"
#include "ismrmrd/xml.h"

namespace Gadgetron{

  gpuCgKtSenseGadget::gpuCgKtSenseGadget()
    : is_configured_(false)
    , channels_(0)
    , frame_counter_(0)
  {
    matrix_size_ = uint64d2(0,0);
    matrix_size_os_ = uint64d2(0,0);
    matrix_size_seq_ = uint64d2(0,0);
  }

  gpuCgKtSenseGadget::~gpuCgKtSenseGadget() {}

  int gpuCgKtSenseGadget::process_config( ACE_Message_Block* mb )
  {
    //GDEBUG("gpuCgKtSenseGadget::process_config\n");

    device_number_ = deviceno.value();

    int number_of_devices = 0;
    if (cudaGetDeviceCount(&number_of_devices)!= cudaSuccess) {
      GDEBUG( "Error: unable to query number of CUDA devices.\n" );
      return GADGET_FAIL;
    }

    if (number_of_devices == 0) {
      GDEBUG( "Error: No available CUDA devices.\n" );
      return GADGET_FAIL;
    }

    if (device_number_ >= number_of_devices) {
      GDEBUG("Adjusting device number from %d to %d\n", device_number_,  (device_number_%number_of_devices));
      device_number_ = (device_number_%number_of_devices);
    }

    if (cudaSetDevice(device_number_)!= cudaSuccess) {
      GDEBUG( "Error: unable to set CUDA device.\n" );
      return GADGET_FAIL;
    }

    set_number_ = setno.value();
    slice_number_ = sliceno.value();
    number_of_iterations_ = number_of_iterations.value();
    cg_limit_ = cg_limit.value();
    oversampling_factor_ = oversampling_factor.value();
    kernel_width_ = kernel_width.value();
    kappa_ = kappa.value();
    shutter_radius_ = training_data_shutter_radius.value();
    rotations_to_discard_ = rotations_to_discard.value();
    output_convergence_ = output_convergence.value();

    if( (rotations_to_discard_%2) == 1 ){
      GDEBUG("#rotations to discard must be even.\n");
      return GADGET_FAIL;
    }

    // Get the Ismrmrd header
    //
    ISMRMRD::IsmrmrdHeader h;
    ISMRMRD::deserialize(mb->rd_ptr(),h);
    
    
    if (h.encoding.size() != 1) {
      GDEBUG("This Gadget only supports one encoding space\n");
      return GADGET_FAIL;
    }
    
    // Get the encoding space and trajectory description
    ISMRMRD::EncodingSpace e_space = h.encoding[0].encodedSpace;
    ISMRMRD::EncodingSpace r_space = h.encoding[0].reconSpace;
    ISMRMRD::EncodingLimits e_limits = h.encoding[0].encodingLimits;

    matrix_size_seq_ = uint64d2( r_space.matrixSize.x, r_space.matrixSize.y );

    if (!is_configured_) {

      if (h.acquisitionSystemInformation) {
	channels_ = h.acquisitionSystemInformation->receiverChannels ? *h.acquisitionSystemInformation->receiverChannels : 1;
      } else {
	channels_ = 1;
      }

      // Allocate encoding operator for non-Cartesian Sense
      E_ = boost::shared_ptr< cuNonCartesianKtSenseOperator<float,2> >( new cuNonCartesianKtSenseOperator<float,2>() );

      // Allocate preconditioner
      D_ = boost::shared_ptr< cuCgPreconditioner<float_complext> >( new cuCgPreconditioner<float_complext>() );

      // Allocate regularization image operator
      R_ = boost::shared_ptr< cuImageOperator<float_complext> >( new cuImageOperator<float_complext>() );
      R_->set_weight( kappa_ );

      // Setup solver
      cg_.set_encoding_operator( E_ );        // encoding matrix
      cg_.add_regularization_operator( R_ );  // regularization matrix
      cg_.set_preconditioner( D_ );           // preconditioning matrix
      cg_.set_max_iterations( number_of_iterations_ );
      cg_.set_tc_tolerance( cg_limit_ );
      cg_.set_output_mode( (output_convergence_) ? cuCgSolver<float_complext>::OUTPUT_VERBOSE : cuCgSolver<float_complext>::OUTPUT_SILENT );

      is_configured_ = true;
    }

    return GADGET_OK;
  }

  int gpuCgKtSenseGadget::process(GadgetContainerMessage<ISMRMRD::ImageHeader> *m1, GadgetContainerMessage<GenericReconJob> *m2)
  {
    // Is this data for this gadget's set/slice?
    //
    
    if( m1->getObjectPtr()->set != set_number_ || m1->getObjectPtr()->slice != slice_number_ ) {      
      // No, pass it downstream...
      return this->next()->putq(m1);
    }
    
    //GDEBUG("gpuCgKtSenseGadget::process\n");
    //GPUTimer timer("gpuCgKtSenseGadget::process");

    if (!is_configured_) {
      GDEBUG("Data received before configuration was completed\n");
      return GADGET_FAIL;
    }

    GenericReconJob* j = m2->getObjectPtr();

    // Some basic validation of the incoming Sense job
    if (!j->csm_host_.get() || !j->dat_host_.get() || !j->tra_host_.get() || !j->dcw_host_.get()) {
      GDEBUG("Received an incomplete Sense job\n");
      return GADGET_FAIL;
    }

    unsigned int samples = j->dat_host_->get_size(0);
    unsigned int channels = j->dat_host_->get_size(1);
    unsigned int rotations = samples / j->tra_host_->get_number_of_elements();
    unsigned int frames = j->tra_host_->get_size(1)*rotations;

    if( samples%j->tra_host_->get_number_of_elements() ) {
      GDEBUG("Mismatch between number of samples (%d) and number of k-space coordinates (%d).\nThe first should be a multiplum of the latter.\n", 
		    samples, j->tra_host_->get_number_of_elements());
      return GADGET_FAIL;
    }

    boost::shared_ptr< cuNDArray<floatd2> > traj(new cuNDArray<floatd2> (*j->tra_host_));
    boost::shared_ptr< cuNDArray<float> > dcw(new cuNDArray<float> (*j->dcw_host_));
    sqrt_inplace(dcw.get()); //Take square root to use for weighting
    boost::shared_ptr< cuNDArray<float_complext> > csm(new cuNDArray<float_complext> (*j->csm_host_));
    boost::shared_ptr< cuNDArray<float_complext> > device_samples(new cuNDArray<float_complext> (*j->dat_host_));

    cudaDeviceProp deviceProp;
    if( cudaGetDeviceProperties( &deviceProp, device_number_ ) != cudaSuccess) {
      GDEBUG( "Error: unable to query device properties.\n" );
      return GADGET_FAIL;
    }
    
    unsigned int warp_size = deviceProp.warpSize;
    
    matrix_size_ = uint64d2( j->reg_host_->get_size(0), j->reg_host_->get_size(1) );    

    matrix_size_os_ =
      uint64d2(((static_cast<unsigned int>(std::ceil(matrix_size_[0]*oversampling_factor_))+warp_size-1)/warp_size)*warp_size,
	     ((static_cast<unsigned int>(std::ceil(matrix_size_[1]*oversampling_factor_))+warp_size-1)/warp_size)*warp_size);
    
    GDEBUG("Matrix size    : [%d,%d] \n", matrix_size_[0], matrix_size_[1]);    
    GDEBUG("Matrix size OS : [%d,%d] \n", matrix_size_os_[0], matrix_size_os_[1]);

    std::vector<size_t> image_dims = to_std_vector(matrix_size_);
    image_dims.push_back(frames);
    
    E_->set_domain_dimensions(image_dims);
    E_->set_codomain_dimensions(device_samples->get_dimensions());
    E_->set_dcw(dcw);
    E_->set_csm(csm);

    E_->setup( matrix_size_, matrix_size_os_, static_cast<float>(kernel_width_) );
    E_->preprocess(traj.get());
        
    R_->compute(compute_regularization_image(j).get());

    // Define preconditioning weights
    boost::shared_ptr< cuNDArray<float> > __precon_weights = sum(abs_square(csm.get()).get(), 2);
    auto _precon_weights = expand<float>( *__precon_weights, frames );
    boost::shared_ptr<cuNDArray<float> > R_diag = R_->get();
    *R_diag *= float(kappa_);
    _precon_weights += *R_diag;
    R_diag.reset();
    reciprocal_sqrt_inplace(&_precon_weights);

    boost::shared_ptr< cuNDArray<float_complext> > precon_weights = real_to_complex<float_complext>( &_precon_weights );
    __precon_weights.reset();
    D_->set_weights( precon_weights );

    *device_samples *= *dcw;
    // Invoke solver
    // 

    boost::shared_ptr< cuNDArray<float_complext> > cgresult;
    
    {
      GPUTimer timer("gpuCgKtSenseGadget::solve()");
      cgresult = cg_.solve(device_samples.get());
    }

    if (!cgresult.get()) {
      GDEBUG("Iterative_sense_compute failed\n");
      return GADGET_FAIL;
    }

    // Goto from x-f to x-t space
    cuNDFFT<float>::instance()->fft( cgresult.get(), 2,true );

    // If the recon matrix size exceeds the sequence matrix size then crop
    if( matrix_size_seq_ != matrix_size_ )
      *cgresult = crop<float_complext,2>( (matrix_size_-matrix_size_seq_)>>1, matrix_size_seq_, *cgresult );
    
    // Now pass on the reconstructed images
    //

    unsigned int frames_per_rotation = frames/rotations;

    if( rotations == 1 ){ // this is the case for golden ratio
      rotations = frames;
      frames_per_rotation = 1;
    }

    for( unsigned int frame=0; frame<frames; frame++ ){

      unsigned int rotation_idx = frame/frames_per_rotation;

      // Check if we should discard this frame
      if( rotation_idx < (rotations_to_discard_>>1) || rotation_idx >= rotations-(rotations_to_discard_>>1) )
	continue;
            
      GadgetContainerMessage<ISMRMRD::ImageHeader> *m = 
	new GadgetContainerMessage<ISMRMRD::ImageHeader>();

      GadgetContainerMessage< hoNDArray< std::complex<float> > > *cm = 
	new GadgetContainerMessage< hoNDArray< std::complex<float> > >();      

      *m->getObjectPtr() = j->image_headers_[frame];
      m->cont(cm);
      
      std::vector<size_t> img_dims(2);
      img_dims[0] = matrix_size_seq_[0];
      img_dims[1] = matrix_size_seq_[1];

      cm->getObjectPtr()->create(img_dims);

      size_t data_length = prod(matrix_size_seq_);

      cudaMemcpy(cm->getObjectPtr()->get_data_ptr(),
		 cgresult->get_data_ptr()+frame*data_length,
		 data_length*sizeof(std::complex<float>),
		 cudaMemcpyDeviceToHost);

      cudaError_t err = cudaGetLastError();
      if( err != cudaSuccess ){
	GDEBUG("Unable to copy result from device to host: %s\n", cudaGetErrorString(err));
	m->release();
	return GADGET_FAIL;
      }

      m->getObjectPtr()->matrix_size[0] = matrix_size_seq_[0];
      m->getObjectPtr()->matrix_size[1] = matrix_size_seq_[1];
      m->getObjectPtr()->matrix_size[2] = 1;
      m->getObjectPtr()->channels       = 1;
      m->getObjectPtr()->image_index    = frame_counter_ + frame;
      
      if (this->next()->putq(m) < 0) {
	GDEBUG("Failed to put result image on to queue\n");
	m->release();
	return GADGET_FAIL;
      }
    }
    
    frame_counter_ += frames;

    m1->release();
    return GADGET_OK;
  }

  boost::shared_ptr< cuNDArray<float_complext> > gpuCgKtSenseGadget::
  compute_regularization_image( GenericReconJob *job )
  {
    // 
    // Estimate training data
    // 

    unsigned int num_samples = job->dat_host_->get_size(0);
    unsigned int num_coils = job->dat_host_->get_size(1);
    unsigned int num_rotations = num_samples / job->tra_host_->get_number_of_elements();
    unsigned int frames_per_reconstruction = job->tra_host_->get_size(1)*num_rotations;

    std::vector<size_t> dims = to_std_vector(matrix_size_os_);
    dims.push_back(frames_per_reconstruction); 
    dims.push_back(num_coils); 

    cuNDArray<float_complext> image_os(dims);
    cuNDArray<float_complext> data(*job->dat_host_);
    cuNDArray<float> dcw(*job->dcw_host_);
  
    // Convolve to Cartesian k-space
    //
    data *= dcw;
    E_->get_plan()->convolve( data, image_os,  NFFT_conv_mode::NC2C );

    // Apply shutter
    //

    if( shutter_radius_ < 0.0001 ){ // If not specified in the configuration then try to make an estimation

      // #profiles/frame : this is just an estimate (we dont have the exact value at this stage)
      unsigned int profiles_per_frame = num_samples / (frames_per_reconstruction*matrix_size_os_[0]);
      shutter_radius_ = ((float)matrix_size_os_[0]/(float)matrix_size_[0])*(float)profiles_per_frame/(float)M_PI;
      GDEBUG("Estimated training data shutter radius: %f\n", shutter_radius_);
    }

    fill_border<float_complext,2>( shutter_radius_, image_os );
    E_->get_plan()->fft( image_os, NFFT_fft_mode::BACKWARDS );
    E_->get_plan()->deapodize( image_os );

    // Remove oversampling
    //

    dims = to_std_vector(matrix_size_);
    dims.push_back(frames_per_reconstruction); 
    dims.push_back(num_coils);
    cuNDArray<float_complext> image(dims);
    crop<float_complext,2>( (matrix_size_os_-matrix_size_)>>1, matrix_size_, image_os, image );

    // Compute regularization image
    //

    dims.pop_back();
    boost::shared_ptr< cuNDArray<float_complext> > reg_image( new cuNDArray<float_complext>(dims) );

    E_->mult_csm_conj_sum( &image, reg_image.get() );
    cuNDFFT<float>::instance()->ifft( reg_image.get(), 2, true );

    return reg_image;
  }

  GADGET_FACTORY_DECLARE(gpuCgKtSenseGadget)
}
