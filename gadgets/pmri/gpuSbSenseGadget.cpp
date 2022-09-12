#include "gpuSbSenseGadget.h"
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

#include <boost/thread/mutex.hpp>
#include <boost/make_shared.hpp>

namespace Gadgetron{

#define max_number_of_gpus 10
  static boost::mutex _mutex[max_number_of_gpus];
  
  gpuSbSenseGadget::gpuSbSenseGadget()
    : is_configured_(false)
    , prepared_(false)
    , gpuSenseGadget()
  {
  }

  gpuSbSenseGadget::~gpuSbSenseGadget() {}

  int gpuSbSenseGadget::process_config( ACE_Message_Block* mb )
  {
    gpuSenseGadget::process_config(mb);
    
    number_of_sb_iterations_ = number_of_sb_iterations.value();
    number_of_cg_iterations_ = number_of_cg_iterations.value();
    cg_limit_ = cg_limit.value();
    mu_ = mu.value();
    lambda_ = lambda.value();
    lambdaT_ = lambdaT.value();
    alpha_ = alpha.value();
    gamma_ = gamma.value();
    exclusive_access_ = exclusive_access.value();
    is_cyclic_= is_cyclic.value();

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
      E_->set_weight(mu_);

      // Allocate preconditioner
      D_ = boost::shared_ptr< cuCgPreconditioner<float_complext> >( new cuCgPreconditioner<float_complext>() );

      Rx1_ = boost::shared_ptr< cuPartialDerivativeOperator<float_complext,3> >
        ( new cuPartialDerivativeOperator<float_complext,3>(0) );
      Rx1_->set_weight( (1.0-alpha_)*lambda_ );

      Ry1_ = boost::shared_ptr< cuPartialDerivativeOperator<float_complext,3> >
        ( new cuPartialDerivativeOperator<float_complext,3>(1) );
      Ry1_->set_weight( (1.0-alpha_)*lambda_ );

      Rz1_ = boost::shared_ptr< cuPartialDerivativeOperator<float_complext,3> >
        ( new cuPartialDerivativeOperator<float_complext,3>(2) );
      Rz1_->set_weight( (1.0-alpha_)*lambda_*lambdaT_ );


      Rt1_ = boost::shared_ptr< cuPartialDerivativeOperator2<float_complext,3> >
              ( new cuPartialDerivativeOperator2<float_complext,3>() );
      Rt1_->set_weight( (1.0-alpha_)*lambda_ *lambdaT_ );

      Rx2_ = boost::shared_ptr< cuPartialDerivativeOperator<float_complext,3> >
        ( new cuPartialDerivativeOperator<float_complext,3>(0) );
      Rx2_->set_weight( alpha_*lambda_ );

      Ry2_ = boost::shared_ptr< cuPartialDerivativeOperator<float_complext,3> >
        ( new cuPartialDerivativeOperator<float_complext,3>(1) );
      Ry2_->set_weight( alpha_*lambda_ );

      Rz2_ = boost::shared_ptr< cuPartialDerivativeOperator<float_complext,3> >
        ( new cuPartialDerivativeOperator<float_complext,3>(2) );
      Rz2_->set_weight( alpha_*lambda_*lambdaT_ );

      Rt2_ = boost::shared_ptr< cuPartialDerivativeOperator2<float_complext,3> >
              ( new cuPartialDerivativeOperator2<float_complext,3>() );
      Rt2_->set_weight( alpha_*lambda_*lambdaT_  );

      W_ = boost::make_shared<cuDWTOperator<float_complext,3>>();
      W_->set_weight(gamma_);
      W2_ = boost::make_shared<cuDWTOperator<float_complext,3>>();
      W2_->set_weight(gamma_);

      // Setup split-Bregman solver
      sb_.set_encoding_operator( E_ );
            
      sb_.set_max_outer_iterations(number_of_sb_iterations_);
      sb_.set_max_inner_iterations(1);
      sb_.set_output_mode( (output_convergence_) ? cuSbcCgSolver<float_complext>::OUTPUT_VERBOSE : cuSbcCgSolver<float_complext>::OUTPUT_SILENT );
      
      sb_.get_inner_solver()->set_max_iterations( number_of_cg_iterations_ );
      sb_.get_inner_solver()->set_tc_tolerance( cg_limit_ );
      //sb_.get_inner_solver()->set_output_mode( (output_convergence_) ? cuCgSolver<float_complext>::OUTPUT_VERBOSE : cuCgSolver<float_complext>::OUTPUT_SILENT );
      sb_.get_inner_solver()->set_preconditioner( D_ );

      is_configured_ = true;
    }

    GDEBUG("gpuSbSenseGadget::end of process_config\n");

    return GADGET_OK;
  }

  int gpuSbSenseGadget::process(GadgetContainerMessage<ISMRMRD::ImageHeader> *m1, GadgetContainerMessage<GenericReconJob> *m2)
  {
    // Is this data for this gadget's set/slice?
    //
    
    if( m1->getObjectPtr()->set != set_number_ || m1->getObjectPtr()->slice != slice_number_ ) {      
      // No, pass it downstream...
      return this->next()->putq(m1);
    }

    //GDEBUG("gpuSbSenseGadget::process\n");
    //GPUTimer timer("gpuSbSenseGadget::process");

    if (!is_configured_) {
      GDEBUG("\nData received before configuration complete\n");
      return GADGET_FAIL;
    }

    GenericReconJob* j = m2->getObjectPtr();

    // Let's first check that this job has the required data...
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

    boost::shared_ptr< cuNDArray<floatd2> > traj(new cuNDArray<floatd2> (j->tra_host_.get()));
    boost::shared_ptr< cuNDArray<float> > dcw(new cuNDArray<float> (j->dcw_host_.get()));
    sqrt_inplace(dcw.get());
    boost::shared_ptr< cuNDArray<float_complext> > csm(new cuNDArray<float_complext> (j->csm_host_.get()));
    boost::shared_ptr< cuNDArray<float_complext> > device_samples(new cuNDArray<float_complext> (j->dat_host_.get()));
    
    if( !prepared_){

      // Take the reconstruction matrix size from the regulariaztion image. 
      // It could be oversampled from the sequence specified size...
      
      matrix_size_ = uint64d2( j->reg_host_->get_size(0), j->reg_host_->get_size(1) );
      
      cudaDeviceProp deviceProp;
      if( cudaGetDeviceProperties( &deviceProp, device_number_ ) != cudaSuccess) {
        GDEBUG( "\nError: unable to query device properties.\n" );
        return GADGET_FAIL;
      }

      unsigned int warp_size = deviceProp.warpSize;

      matrix_size_os_ =
        uint64d2(((static_cast<unsigned int>(std::ceil(matrix_size_[0]*oversampling_factor_))+warp_size-1)/warp_size)*warp_size,
                 ((static_cast<unsigned int>(std::ceil(matrix_size_[1]*oversampling_factor_))+warp_size-1)/warp_size)*warp_size);
      
      GDEBUG("Matrix size    : [%d,%d] \n", matrix_size_[0], matrix_size_[1]);
      GDEBUG("Matrix size OS : [%d,%d] \n", matrix_size_os_[0], matrix_size_os_[1]);

      std::vector<size_t> image_dims = to_std_vector(matrix_size_);
      image_dims.push_back(frames);
      
      E_->set_domain_dimensions(&image_dims);
      E_->set_codomain_dimensions(device_samples->get_dimensions().get());
            
      reg_image_ = boost::shared_ptr< cuNDArray<float_complext> >(new cuNDArray<float_complext>(&image_dims));
      
      // These operators need their domain/codomain set before being added to the solver
      //

      Rx1_->set_domain_dimensions(&image_dims);
      Rx1_->set_codomain_dimensions(&image_dims);
      
      Ry1_->set_domain_dimensions(&image_dims);
      Ry1_->set_codomain_dimensions(&image_dims);
      
      Rz1_->set_domain_dimensions(&image_dims);
      Rz1_->set_codomain_dimensions(&image_dims);
      
      Rt1_->set_domain_dimensions(&image_dims);
      Rt1_->set_codomain_dimensions(&image_dims);

      Rx2_->set_domain_dimensions(&image_dims);
      Rx2_->set_codomain_dimensions(&image_dims);
      
      Ry2_->set_domain_dimensions(&image_dims);
      Ry2_->set_codomain_dimensions(&image_dims);

      Rt2_->set_domain_dimensions(&image_dims);
      Rt2_->set_codomain_dimensions(&image_dims);

      Rz2_->set_domain_dimensions(&image_dims);
      Rz2_->set_codomain_dimensions(&image_dims);

      W_->set_domain_dimensions(&image_dims);
      W_->set_codomain_dimensions(&image_dims);
      W2_->set_domain_dimensions(&image_dims);
      W2_->set_codomain_dimensions(&image_dims);
      W2_->set_shift(2);
      
      // Add "TV" regularization
      // 
      
      if( alpha_<1.0 ){
        sb_.add_regularization_group_operator( Rx1_ ); 
        sb_.add_regularization_group_operator( Ry1_ ); 
        if(frames>1)
        	if (is_cyclic_)
        		sb_.add_regularization_group_operator( Rz1_ );
        	else
        		sb_.add_regularization_group_operator( Rt1_ );
        sb_.add_group();
      }
      
      // Add "PICCS" regularization
      //

      if( alpha_ > 0.0 ){
        sb_.add_regularization_group_operator( Rx2_ ); 
        sb_.add_regularization_group_operator( Ry2_ ); 
        if(frames>1)
        	if (is_cyclic_)
        		sb_.add_regularization_group_operator( Rz2_ );
        	else
        		sb_.add_regularization_group_operator( Rt2_ );
        sb_.add_group(reg_image_);
      }
      
      if (gamma_ > 0.0){
    	  sb_.add_regularization_operator(W_);
    	  sb_.add_regularization_operator(W2_);
      }
      prepared_ = true;
    }
    
    E_->set_dcw(dcw);
    E_->set_csm(csm);    
    E_->setup( matrix_size_, matrix_size_os_, static_cast<float>(kernel_width_) );
    E_->preprocess(traj.get());

    // Expand the average image to the number of frames
    //

    {
      cuNDArray<float_complext> tmp(*j->reg_host_);
      *reg_image_ = expand( tmp, frames );
    }

    // Define preconditioning weights
    //

    boost::shared_ptr< cuNDArray<float> > _precon_weights = sum(abs_square(csm.get()).get(), 2);
    reciprocal_sqrt_inplace(_precon_weights.get());	
    boost::shared_ptr< cuNDArray<float_complext> > precon_weights = real_to_complex<float_complext>( _precon_weights.get() );
    _precon_weights.reset();
    D_->set_weights( precon_weights );
    precon_weights.reset();
    
    //Apply weights
    *device_samples *= *dcw;

    // Invoke solver
    //

    boost::shared_ptr< cuNDArray<float_complext> > sbresult;
    {
      GDEBUG("Running split Bregman solver\n");
      GPUTimer timer("Running split Bregman solver");

      // Optionally, allow exclusive (per device) access to the solver
      // This may not matter much in terms of speed, but it can in terms of memory consumption
      //

      if( exclusive_access_ )
        _mutex[device_number_].lock();

      sbresult = sb_.solve(device_samples.get());

      if( exclusive_access_ )
        _mutex[device_number_].unlock();
    }

    // Provide some info about the scaling between the regularization and reconstruction.
    // If it is not close to one, PICCS does not work optimally...
    // 

    if( alpha_ > 0.0 ){
      cuNDArray<float_complext> gpureg(j->reg_host_.get());
      boost::shared_ptr< cuNDArray<float_complext> > gpurec = sum(sbresult.get(),2);
      *gpurec /= float(sbresult->get_size(2));
      float scale = abs(dot(gpurec.get(), gpurec.get())/dot(gpurec.get(),&gpureg));
      GDEBUG("Scaling factor between regularization and reconstruction is %f.\n", scale);
    }
    
    if (!sbresult.get()) {
      GDEBUG("\nSplit Bregman solver failed\n");
      return GADGET_FAIL;
    }
    
    /*
      static int counter = 0;
      char filename[256];
      sprintf((char*)filename, "recon_sb_%d.cplx", counter);
      write_nd_array<float_complext>( sbresult->to_host().get(), filename );
      counter++; */

    // If the recon matrix size exceeds the sequence matrix size then crop
    if( matrix_size_seq_ != matrix_size_ )
      *sbresult = crop<float_complext,2>( (matrix_size_-matrix_size_seq_)>>1, matrix_size_seq_, *sbresult );
        
    // Now pass on the reconstructed images
    //

	put_frames_on_que(frames,rotations,j,sbresult.get());

    frame_counter_ += frames;
    m1->release();
    return GADGET_OK;
  }

  GADGET_FACTORY_DECLARE(gpuSbSenseGadget)
}

