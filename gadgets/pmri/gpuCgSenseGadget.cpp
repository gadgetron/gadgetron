#include "gpuCgSenseGadget.h"
#include "cuNDArray_operators.h"
#include "cuNDArray_elemwise.h"
#include "cuNDArray_blas.h"
#include "cuNDArray_utils.h"
#include "cuNDArray_reductions.h"
#include "GadgetMRIHeaders.h"
#include "b1_map.h"
#include "GPUTimer.h"
#include "vector_td_utilities.h"
#include "hoNDArray_fileio.h"
#include "ismrmrd/xml.h"

namespace Gadgetron{

  gpuCgSenseGadget::gpuCgSenseGadget()
    : is_configured_(false)
    , matrix_size_reported_(0),
      gpuSenseGadget()
  {
    matrix_size_ = uint64d2(0,0);
    matrix_size_os_ = uint64d2(0,0);
    matrix_size_seq_ = uint64d2(0,0);
  }

  gpuCgSenseGadget::~gpuCgSenseGadget() {}

  int gpuCgSenseGadget::process_config( ACE_Message_Block* mb )
  {
  	gpuSenseGadget::process_config(mb);
    //GDEBUG("gpuCgSenseGadget::process_config\n");
    number_of_iterations_ = number_of_iterations.value();
    kappa_ = kappa.value();
    cg_limit_ = cg_limit.value();



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
      E_ = boost::shared_ptr< cuNonCartesianSenseOperator<float,2> >( new cuNonCartesianSenseOperator<float,2>() );

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
      cg_.set_output_mode( (output_convergence_) ? cuCgSolver<float_complext>::OUTPUT_VERBOSE : cuCgSolver<float_complext>::OUTPUT_SILENT);
      is_configured_ = true;
    }

    return GADGET_OK;
  }

  int gpuCgSenseGadget::process(GadgetContainerMessage<ISMRMRD::ImageHeader> *m1, GadgetContainerMessage<GenericReconJob> *m2)
  {
    // Is this data for this gadget's set/slice?
    //
    
    if( m1->getObjectPtr()->set != set_number_ || m1->getObjectPtr()->slice != slice_number_ ) {      
      // No, pass it downstream...
      return this->next()->putq(m1);
    }
    
    //GDEBUG("gpuCgSenseGadget::process\n");

    boost::shared_ptr<GPUTimer> process_timer;
    if( output_timing_ )
      process_timer = boost::shared_ptr<GPUTimer>( new GPUTimer("gpuCgSenseGadget::process()") );
    
    if (!is_configured_) {
      GDEBUG("Data received before configuration was completed\n");
      return GADGET_FAIL;
    }

    GenericReconJob* j = m2->getObjectPtr();

    // Some basic validation of the incoming Sense job
    if (!j->csm_host_.get() || !j->dat_host_.get() || !j->tra_host_.get() || !j->dcw_host_.get() || !j->reg_host_.get()) {
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

    boost::shared_ptr< cuNDArray<floatd2> > traj(new cuNDArray<floatd2> (j->tra_host_.get()));
    boost::shared_ptr< cuNDArray<float> > dcw(new cuNDArray<float> (j->dcw_host_.get()));
    sqrt_inplace(dcw.get()); //Take square root to use for weighting
    boost::shared_ptr< cuNDArray<float_complext> > csm(new cuNDArray<float_complext> (j->csm_host_.get()));
    boost::shared_ptr< cuNDArray<float_complext> > device_samples(new cuNDArray<float_complext> (j->dat_host_.get()));
    
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

    if( !matrix_size_reported_ ) {
      GDEBUG("Matrix size    : [%d,%d] \n", matrix_size_[0], matrix_size_[1]);
      GDEBUG("Matrix size OS : [%d,%d] \n", matrix_size_os_[0], matrix_size_os_[1]);
      matrix_size_reported_ = true;
    }

    std::vector<size_t> image_dims = to_std_vector(matrix_size_);
    image_dims.push_back(frames);
    
    E_->set_domain_dimensions(&image_dims);
    E_->set_codomain_dimensions(device_samples->get_dimensions().get());
    E_->set_dcw(dcw);
    E_->set_csm(csm);

    E_->setup( matrix_size_, matrix_size_os_, static_cast<float>(kernel_width_) );
    E_->preprocess(traj.get());

    boost::shared_ptr< cuNDArray<float_complext> > reg_image(new cuNDArray<float_complext> (j->reg_host_.get()));
    R_->compute(reg_image.get());

    // Define preconditioning weights
    boost::shared_ptr< cuNDArray<float> > _precon_weights = sum(abs_square(csm.get()).get(), 2);
    boost::shared_ptr<cuNDArray<float> > R_diag = R_->get();
    *R_diag *= float(kappa_);
    *_precon_weights += *R_diag;
    R_diag.reset();
    reciprocal_sqrt_inplace(_precon_weights.get());	
    boost::shared_ptr< cuNDArray<float_complext> > precon_weights = real_to_complex<float_complext>( _precon_weights.get() );
    _precon_weights.reset();
    D_->set_weights( precon_weights );
    
    //Apply dcw weights
    *device_samples *= *dcw;


    /*{
      static int counter = 0;
      char filename[256];
      sprintf((char*)filename, "_traj_%d.real", counter);
      write_nd_array<floatd2>( traj->to_host().get(), filename );
      sprintf((char*)filename, "_dcw_%d.real", counter);
      write_nd_array<float>( dcw->to_host().get(), filename );
      sprintf((char*)filename, "_csm_%d.cplx", counter);
      write_nd_array<float_complext>( csm->to_host().get(), filename );
      sprintf((char*)filename, "_samples_%d.cplx", counter);
      write_nd_array<float_complext>( device_samples->to_host().get(), filename );
      sprintf((char*)filename, "_reg_%d.cplx", counter);
      write_nd_array<float_complext>( reg_image->to_host().get(), filename );
      counter++; 
      }*/

    // Invoke solver
    // 
    boost::shared_ptr< cuNDArray<float_complext> > cgresult;

    {
      boost::shared_ptr<GPUTimer> solve_timer;
      if( output_timing_ )
        solve_timer = boost::shared_ptr<GPUTimer>( new GPUTimer("gpuCgSenseGadget::solve()") );
      
      cgresult = cg_.solve(device_samples.get());
      
      if( output_timing_ )
        solve_timer.reset();
    }
    
    if (!cgresult.get()) {
      GDEBUG("Iterative_sense_compute failed\n");
      return GADGET_FAIL;
    }

    /*
      static int counter = 0;
      char filename[256];
      sprintf((char*)filename, "recon_%d.real", counter);
      write_nd_array<float>( abs(cgresult.get())->to_host().get(), filename );
      counter++; 
    */

    // If the recon matrix size exceeds the sequence matrix size then crop
    if( matrix_size_seq_ != matrix_size_ )
      *cgresult = crop<float_complext,2>( (matrix_size_-matrix_size_seq_)>>1, matrix_size_seq_, *cgresult );
    
    // Now pass on the reconstructed images
    //
    put_frames_on_que(frames,rotations,j,cgresult.get());
    
    frame_counter_ += frames;

    if( output_timing_ )
      process_timer.reset();

    m1->release();
    return GADGET_OK;
  }

  GADGET_FACTORY_DECLARE(gpuCgSenseGadget)
}
