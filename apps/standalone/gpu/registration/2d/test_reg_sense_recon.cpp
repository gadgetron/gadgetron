#define PAD_Z

/*

  This is an example of how to use optical flow image registration 
  and the image resampling operator for image reconstruction.
  
  This example uses golden ratio Sense MRI for demonstration. 
  It was tested with a free-breathing cardiac acquisition.

  !!! Note !!!
  ------------
  No cardiac phase binning is performed.
  And since the registration has trouble handling large, 
  non-rigid deformations such as the heart contraction
  it serves only for demonstration purposes. 

  An actual application should bin the cardiac phases and use the 
  registration to correct for respiratory motion only.
*/

#include "cuCKOpticalFlowSolver.h"
#include "cuLinearResampleOperator.h"
#include "cuNDArray.h"
#include "cuNDArray_elemwise.h"
#include "cuNDArray_operators.h"
#include "cuNDArray_utils.h"
#include "cuNDArray_reductions.h"
#include "hoNDArray_fileio.h"
#include "parameterparser.h"
#include "radial_utilities.h"
#include "cuNonCartesianSenseOperator.h"
#include "cuSenseBuffer.h"
#include "cuImageOperator.h"
#include "cuCgPreconditioner.h"
#include "cuCgSolver.h"
#include "b1_map.h"
#include "GPUTimer.h"

#include <iostream>

using namespace Gadgetron;
using namespace std;

// Define desired precision
//

typedef float _real; 
typedef complext<_real> _complext;
typedef reald<_real,2>::Type _reald2;

//
// Define matrix operator for "registration reconstruction" using non-Cartesian Sense
// For simplicity we assume that the respective operators have been setup from outside
//

template<class REAL, unsigned int D> class registrationReconOperator
  : public linearOperator< cuNDArray< complext<REAL> > >
{
public:
  
  registrationReconOperator() : linearOperator< cuNDArray< complext< REAL> > >() {}
  virtual ~registrationReconOperator() {}
  
  inline void set_encoding_operator( boost::shared_ptr< cuNonCartesianSenseOperator<REAL,D> > E ){
    E_ = E;
  }
  
  inline void set_resampling_operator( boost::shared_ptr< cuLinearResampleOperator<complext<REAL>,D> > R ){
    R_ = R;
  }
  
  virtual void mult_M( cuNDArray< complext<REAL> > *in, cuNDArray< complext<REAL> > *out, bool accumulate = false)
  {
    if( !in || !out || !R_->get_displacement_field() ){
      throw cuda_error("registrationReconOperator::mult_M failed (1)");
    }
    
    // Allocate intermediate image
    std::vector<size_t> tmp_dims = *R_->get_displacement_field()->get_dimensions(); tmp_dims.pop_back();
    cuNDArray< complext<REAL> > tmp_in_out;

    tmp_in_out.create(&tmp_dims);
    
    // Deform the input image into multiple frames by applying the registration vector field
    R_->mult_M( in, &tmp_in_out );

    // Apply non-Cartesian Sense encoding
    E_->mult_M( &tmp_in_out, out );
  }
  
  virtual void mult_MH( cuNDArray< complext<REAL> > *in, cuNDArray< complext<REAL> > *out, bool accumulate = false )
  {
    if( !in || !out || !R_->is_preprocessed() ){
      throw cuda_error("registrationReconOperator::mult_MH failed (1)");
    }
    
    // Allocate intermediate image
    std::vector<size_t> tmp_dims = *R_->get_displacement_field()->get_dimensions().get(); tmp_dims.pop_back();
    cuNDArray< complext<REAL> > tmp_in_out(&tmp_dims); 

    // Apply adjoint non-Cartesian Sense encoding
    E_->mult_MH( in, &tmp_in_out);
  
    // Apply adjoint registration
    R_->mult_MH( &tmp_in_out, out );
  }
  
  virtual void mult_MH_M( cuNDArray< complext<REAL> > *in, cuNDArray< complext<REAL> > *out, bool accumulate = false )
  {
    if( !in || !out || !R_->get_displacement_field() ){
      throw cuda_error("registrationReconOperator::mult_MH_M failed (1)");
    }

    // Allocate intermediate image
    std::vector<size_t> tmp_dims = *R_->get_displacement_field()->get_dimensions().get(); tmp_dims.pop_back();
    cuNDArray< complext<REAL> > tmp_in_out1(&tmp_dims), tmp_in_out2(&tmp_dims); 
    
    // Deform the input image into multiple frames by applying the registration vector field
    R_->mult_M( in, &tmp_in_out1 );

    // Apply non-Cartesian Sense encoding _iteration_
    E_->mult_MH_M( &tmp_in_out1, &tmp_in_out2 );
    
    // Apply adjoint registration
    R_->mult_MH( &tmp_in_out2, out );
  }
  
  virtual boost::shared_ptr< linearOperator< cuNDArray< complext<REAL> > > > clone() {
    return linearOperator< cuNDArray<complext< REAL > > >::clone(this);
  }
  
private:
  boost::shared_ptr< cuNonCartesianSenseOperator<REAL,D> > E_;
  boost::shared_ptr< cuLinearResampleOperator<complext<REAL>,D> > R_;
};


//
// Utility to upload samples for one reconstruction from host to device
//

boost::shared_ptr< cuNDArray<_complext> > 
upload_data( unsigned int reconstruction, unsigned int samples_per_reconstruction, unsigned int total_samples_per_coil, unsigned int num_coils, hoNDArray<_complext> *host_data, unsigned int offset = 0 )
{
  vector<size_t> dims; dims.push_back(samples_per_reconstruction); dims.push_back(num_coils);
  cuNDArray<_complext> *data = new cuNDArray<_complext>(); data->create( &dims );
  for( unsigned int i=0; i<num_coils; i++ )
    cudaMemcpy( data->get_data_ptr()+i*samples_per_reconstruction, 
		host_data->get_data_ptr()+i*total_samples_per_coil+reconstruction*samples_per_reconstruction+offset, 
		samples_per_reconstruction*sizeof(_complext), cudaMemcpyHostToDevice );

  return boost::shared_ptr< cuNDArray<_complext> >(data);
}

int main(int argc, char** argv)
{

  // Parse command line
  //

  ParameterParser parms;
  parms.add_parameter( 'd', COMMAND_LINE_STRING, 1, "MRI sample data file name", true, "fb_data.cplx" );
  parms.add_parameter( 'r', COMMAND_LINE_STRING, 1, "Reconstruction result file name", true, "result.real" );

  // Parameters for the initial Sense reconstruction
  //

  parms.add_parameter( 'm', COMMAND_LINE_INT,    1, "Matrix size", true, "256" );
  parms.add_parameter( 'o', COMMAND_LINE_INT,    1, "Oversampled matrix size", true, "384" );
  parms.add_parameter( 'p', COMMAND_LINE_INT,    1, "Profiles per frame", true, "16" );
  parms.add_parameter( 'i', COMMAND_LINE_INT,    1, "Number of iterations", true, "15" );
  parms.add_parameter( 'k', COMMAND_LINE_FLOAT,  1, "Kernel width", true, "5.5" );
  parms.add_parameter( 'K', COMMAND_LINE_FLOAT,  1, "Kappa", true, "0.1" );

  // Parameters for the registration
  //

  parms.add_parameter( 'a', COMMAND_LINE_FLOAT,  1, "Registration regularization weight (alpha)", true, "0.05" );
  parms.add_parameter( 'b', COMMAND_LINE_FLOAT,  1, "Registration regularization weight (beta)", true, "1.0" );
  
  parms.parse_parameter_list(argc, argv);
  if( parms.all_required_parameters_set() ){
    cout << " Running registration with the following parameters: " << endl;
    parms.print_parameter_list();
  }
  else{
    cout << " Some required parameters are missing: " << endl;
    parms.print_parameter_list();
    parms.print_usage();
    return 1;
  }
  
  GPUTimer *timer = new GPUTimer("\nPerforming Sense reconstruction");

  //
  // First perform the Sense reconstruction, 
  // resulting in aliased presumably...
  //
  
  // Load sample data from disk
  //
  
  boost::shared_ptr< hoNDArray<_complext> > host_data = read_nd_array<_complext>((char*)parms.get_parameter('d')->get_string_value());
   
  if( !(host_data->get_number_of_dimensions() == 3) ){
    cout << endl << "Input data is not three-dimensional (#samples/profile x #profiles x #coils). Quitting!" << endl;
    return 1;
  }
  
  // Configuration from the host data
  //

  unsigned int samples_per_profile = host_data->get_size(0);
  unsigned int num_profiles = host_data->get_size(1);
  unsigned int num_coils = host_data->get_size(2);
  
  // Configuration from the command line
  //

  uint64d2 matrix_size = uint64d2(parms.get_parameter('m')->get_int_value(), parms.get_parameter('m')->get_int_value());
  uint64d2 matrix_size_os = uint64d2(parms.get_parameter('o')->get_int_value(), parms.get_parameter('o')->get_int_value());
  _real kernel_width = parms.get_parameter('k')->get_float_value();
  _real kappa = parms.get_parameter('K')->get_float_value();
  unsigned int num_iterations = parms.get_parameter('i')->get_int_value();
  unsigned int profiles_per_frame = parms.get_parameter('p')->get_int_value();
  unsigned int frames_per_reconstruction = 1;

  // Silent correction of invalid command line parameters (clamp to valid range)
  //

  if( profiles_per_frame > num_profiles ) profiles_per_frame = num_profiles;
  if( frames_per_reconstruction < 0 ) frames_per_reconstruction = num_profiles / profiles_per_frame;
  if( frames_per_reconstruction*profiles_per_frame > num_profiles ) frames_per_reconstruction = num_profiles / profiles_per_frame;
  
  unsigned int profiles_per_reconstruction = frames_per_reconstruction*profiles_per_frame;
  unsigned int samples_per_frame = profiles_per_frame*samples_per_profile;
  unsigned int samples_per_reconstruction = profiles_per_reconstruction*samples_per_profile;

  // Set density compensation weights
  //

  boost::shared_ptr< cuNDArray<_real> > dcw = compute_radial_dcw_golden_ratio_2d
    ( samples_per_profile, profiles_per_frame, (_real)matrix_size_os.vec[0]/(_real)matrix_size.vec[0], 
      _real(1)/((_real)samples_per_profile/(_real)max(matrix_size.vec[0],matrix_size.vec[1])) );

  // Define encoding matrix for non-Cartesian SENSE
  //

  boost::shared_ptr< cuNonCartesianSenseOperator<_real,2> > E( new cuNonCartesianSenseOperator<_real,2>() );  
  E->setup( matrix_size, matrix_size_os, kernel_width );
  
  std::vector<size_t> tmp_vec = to_std_vector(matrix_size);
  tmp_vec.push_back(frames_per_reconstruction);
  E->set_domain_dimensions( &tmp_vec );

  // Notify encoding operator of dcw
  //
  
  E->set_dcw(dcw);
  
  // Define rhs buffer
  //

  boost::shared_ptr< cuSenseBuffer<_real,2> > rhs_buffer( new cuSenseBuffer<_real,2>() );
  rhs_buffer->setup( matrix_size, matrix_size_os, kernel_width, num_coils, 8, 16 );
  rhs_buffer->set_dcw(dcw);
  
  // Fill rhs buffer (go through all the data...)
  //
    
  for( unsigned int iteration = 0; iteration < num_profiles/profiles_per_frame; iteration++ ) {

    // Define trajectories
    boost::shared_ptr< cuNDArray<_reald2> > traj = compute_radial_trajectory_golden_ratio_2d<_real>
      ( samples_per_profile, profiles_per_frame, 1, iteration*profiles_per_frame );
    
    // Upload data
    boost::shared_ptr< cuNDArray<_complext> > csm_data = upload_data
      ( iteration, samples_per_frame, num_profiles*samples_per_profile, num_coils, host_data.get() );
    
    // Add frame to rhs buffer
    rhs_buffer->add_frame_data( csm_data.get(), traj.get() );
  }

  // Estimate CSM
  //

  boost::shared_ptr< cuNDArray<_complext> > acc_images = rhs_buffer->get_accumulated_coil_images();
  boost::shared_ptr< cuNDArray<_complext> > csm = estimate_b1_map<_real,2>( acc_images.get() );

  E->set_csm(csm);

  // Define regularization image operator 
  //

  std::vector<size_t> image_dims = to_std_vector(matrix_size);
  cuNDArray<_complext> *regul_image = new cuNDArray<_complext>(&image_dims);
  
  E->mult_csm_conj_sum( acc_images.get(), regul_image );
  acc_images.reset();

  boost::shared_ptr< cuImageOperator<_complext> > R( new cuImageOperator<_complext>() ); 
  R->set_weight( kappa );
  R->compute( regul_image ); 
  delete regul_image; regul_image = 0x0;

  // Define preconditioning weights
  //

  boost::shared_ptr< cuNDArray<_real> > _precon_weights = sum(abs_square(csm.get()).get(),2);
  boost::shared_ptr< cuNDArray<_real> > R_diag = R->get();
  *R_diag *= kappa;
  *_precon_weights += *R_diag;
  R_diag.reset();
  reciprocal_sqrt_inplace(_precon_weights.get());
  boost::shared_ptr< cuNDArray<_complext> > precon_weights = real_to_complex<_complext>( _precon_weights.get() );
  _precon_weights.reset();

  // Define preconditioning matrix
  //

  boost::shared_ptr< cuCgPreconditioner<_complext> > D( new cuCgPreconditioner<_complext>() );
  D->set_weights( precon_weights );
  precon_weights.reset();
  csm.reset();
  
  // Setup radial SENSE reconstructions (conjugate gradient solver)
  //
      
  cuCgSolver<_complext> *cg = new cuCgSolver<_complext>;
  cg->set_encoding_operator( E );  // encoding matrix
  cg->add_regularization_operator( R );  // regularization matrix
  cg->set_preconditioner ( D );  // preconditioning matrix
  cg->set_max_iterations( num_iterations );
  cg->set_tc_tolerance( 1e-6 );
  cg->set_output_mode( cuCgSolver<_complext>::OUTPUT_VERBOSE );

  // To save memory we allow only a certain number of frames
  unsigned int max_num_frames = 25;
  unsigned int reconstruction_offset = 100; // To find some respiratory movement in the test dataset
  unsigned int num_reconstructions = num_profiles / profiles_per_reconstruction;
  if( num_reconstructions<(max_num_frames+reconstruction_offset) ) reconstruction_offset = 0;
  if( num_reconstructions > max_num_frames ) num_reconstructions = max_num_frames;
  
  // Allocate space for aliased reconstruction
  //
  
  image_dims = to_std_vector(matrix_size); 
  image_dims.push_back(frames_per_reconstruction*num_reconstructions); 
  cuNDArray<_complext> *sense_result_cplx = new cuNDArray<_complext>; 
  GDEBUG_STREAM(std::endl << matrix_size[0] << " " << matrix_size[1] << " " << frames_per_reconstruction << " " << num_reconstructions);

  sense_result_cplx->create(&image_dims);
  
  // Loop and reconstruct 
  // 

  for( unsigned int reconstruction = 0; reconstruction<num_reconstructions; reconstruction++ ){
    
    // Determine trajectories
    //

    boost::shared_ptr< cuNDArray<_reald2> > traj = compute_radial_trajectory_golden_ratio_2d<_real>
      ( samples_per_profile, profiles_per_frame, frames_per_reconstruction, (reconstruction+reconstruction_offset)*profiles_per_reconstruction );
    
    // Upload data
    //

    boost::shared_ptr< cuNDArray<_complext> > data = upload_data
      ( reconstruction+reconstruction_offset, samples_per_reconstruction, num_profiles*samples_per_profile, num_coils, host_data.get() );
    
    // Set current trajectory and trigger NFFT preprocessing
    //
    
    E->preprocess(traj.get());
    
    // Form rhs (use sense_result_cplx array to save memory)
    //
    
    vector<size_t> rhs_dims = to_std_vector(matrix_size); 
    rhs_dims.push_back(frames_per_reconstruction);
    cuNDArray<_complext> rhs; 

    rhs.create( &rhs_dims, sense_result_cplx->get_data_ptr()+
		reconstruction*prod(matrix_size)*frames_per_reconstruction );

    E->mult_MH( data.get(), &rhs );
    
    // Conjugate gradient solver
    //

    boost::shared_ptr< cuNDArray<_complext> > cgresult = cg->solve(data.get());
    rhs = *(cgresult.get());
  }
  
  boost::shared_ptr< cuNDArray<_real> > sense_result = abs(sense_result_cplx);
  write_nd_array<_complext>(sense_result_cplx->to_host().get(), "_images_all.cplx");

  // We need all our device memory for the registration. Clean up after Sense.
  // E.reset(); D.reset(); R.reset();   -- we will reuse these below 

  rhs_buffer.reset();
  delete cg; delete sense_result_cplx;
  delete timer;

  // Determine fixed/moving image dimensions and create arrays
  //

#ifdef PAD_Z
  std::vector<size_t> _3d_dims = *(sense_result->get_dimensions());
  unsigned int last_dim = _3d_dims.back();
  _3d_dims.pop_back(); _3d_dims.push_back(1); _3d_dims.push_back(last_dim);
  sense_result->reshape( &_3d_dims );
#endif
  
  vector<size_t> multi_dims = *sense_result->get_dimensions();
  multi_dims.pop_back();
#ifdef PAD_Z
  multi_dims.push_back(sense_result->get_size(3)-1);
#else
  multi_dims.push_back(sense_result->get_size(2)-1);
#endif
  vector<size_t> single_dims = *sense_result->get_dimensions();
  single_dims.pop_back();
  
  cuNDArray<_real> 
    *multi_image = new cuNDArray<_real>, 
    *single_image = new cuNDArray<_real>;
  
  single_image->create( &single_dims, sense_result->get_data_ptr());
  multi_image->create( &multi_dims, sense_result->get_data_ptr()+prod(matrix_size));
  
  write_nd_array<_real>(multi_image->to_host().get(), "_images_multi.real");
  write_nd_array<_real>(single_image->to_host().get(), "_image_single.real");

  // Setup registration solver
  //
#ifdef PAD_Z
  cuCKOpticalFlowSolver<_real,3> *CK = new cuCKOpticalFlowSolver<_real,3>;
#else
  cuCKOpticalFlowSolver<_real,2> *CK = new cuCKOpticalFlowSolver<_real,2>;
#endif

  //CK->set_output_mode( cuCKOpticalFlowSolver<_real,2>::OUTPUT_VERBOSE );  
  CK->set_num_multires_levels( 1 );
  CK->set_max_num_iterations_per_level( 500 );
  CK->set_alpha((_real) parms.get_parameter('a')->get_float_value());
  CK->set_beta((_real) parms.get_parameter('b')->get_float_value());
  CK->set_limit(0.01f);
  
  // 
  // Perform "averaging by registration" type reconstruction
  //

  timer = new GPUTimer("\nReconstruction by optical flow averaging");

  // Run registration:
  // - multi_image -> single_image (many to one registration)
  // 

  // All to one
  boost::shared_ptr< cuNDArray<_real> > reg_result = CK->solve( single_image, multi_image );
  
  write_nd_array<_real>(reg_result->to_host().get(), "_reg1.real");

  // Deform the multi_image according to the deformation field and average
  //

  boost::shared_ptr< cuNDArray<_real> > regis_image = CK->deform( multi_image, reg_result );
#ifdef PAD_Z
  boost::shared_ptr< cuNDArray<_real> > regis_image_avg = sum<_real>( regis_image.get(), 3); 
#else
  boost::shared_ptr< cuNDArray<_real> > regis_image_avg = sum<_real>( regis_image.get(), 2); 
#endif
  write_nd_array<_real>(regis_image->to_host().get(), "_reg_avg.real");
  write_nd_array<_real>(regis_image_avg->to_host().get(), "_avg_recon.real");

  regis_image.reset(); regis_image_avg.reset(); reg_result.reset();

  delete timer;

  //
  // Perform "registration in cost function" type reconstruction
  //

  timer = new GPUTimer("\nRunning registration recon");

  // One to all
  reg_result = CK->solve( multi_image, single_image );
  
  write_nd_array<_real>(reg_result->to_host().get(), "_reg2.real");

  regis_image = CK->deform( single_image, reg_result );
  write_nd_array<_real>(regis_image->to_host().get(), "_multi_def.real");
  regis_image.reset(); 

  // Test iteration
  cuNDArray<_real> out; out.create(multi_image->get_dimensions().get());
  cuNDArray<_real> in; in.create(single_image->get_dimensions().get());
  
  // Release memory
  delete CK;
  exit(1);
  // Setup solver
  //

  // The non-Cartesian Sense operator is already setup, 
  // but the trajectories must be recomputed and preprocessed

  boost::shared_ptr< cuNDArray<_reald2> >traj = compute_radial_trajectory_golden_ratio_2d<_real>
    ( samples_per_profile, profiles_per_frame, frames_per_reconstruction*(num_reconstructions-1), 
      (1+reconstruction_offset)*profiles_per_reconstruction );
  
  E->preprocess(traj.get());

  // Define and preprocess resampling operator
  
  boost::shared_ptr< cuLinearResampleOperator<_complext,2> > resampler
    ( new cuLinearResampleOperator<_complext,2> );

  resampler->set_displacement_field(reg_result);
  resampler->mult_MH_preprocess();

  // Define registrationReconstruction encoding operator

  boost::shared_ptr< registrationReconOperator<_real,2> > 
    RR( new registrationReconOperator<_real,2>() );  

  std::vector<size_t> rhs_dims = to_std_vector(matrix_size); 
  RR->set_domain_dimensions( &rhs_dims );

  RR->set_encoding_operator( E );
  RR->set_resampling_operator( resampler );

  cg = new cuCgSolver<_complext>;
  cg->set_encoding_operator( RR );
  cg->add_regularization_operator( R );
  cg->set_preconditioner ( D ); 
  cg->set_max_iterations( num_iterations );
  cg->set_tc_tolerance( 1e-6 );
  cg->set_output_mode( cuCgSolver<_complext>::OUTPUT_VERBOSE );

  // Form rhs
  
  boost::shared_ptr< cuNDArray<_complext> > data = upload_data
    ( 0, samples_per_reconstruction*(num_reconstructions-1), 
      num_profiles*samples_per_profile, num_coils, host_data.get(), 
      (reconstruction_offset+1)*samples_per_reconstruction );
  
  cuNDArray<_complext> rhs(&rhs_dims); 
  RR->mult_MH( data.get(), &rhs );
  
  write_nd_array<_complext>(rhs.to_host().get(), "_rhs.cplx" );
  write_nd_array<_real>(abs(&rhs)->to_host().get(), "_rhs.real" );
 
  // Conjugate gradient solver
  //
  
  boost::shared_ptr< cuNDArray<_complext> > cgresult = cg->solve(data.get());

  boost::shared_ptr< hoNDArray<_real> > host_image = abs(cgresult.get())->to_host();
  write_nd_array<_real>(host_image.get(), "_reg_frame.real" );
  
  delete timer;
   
  return 0;
}
