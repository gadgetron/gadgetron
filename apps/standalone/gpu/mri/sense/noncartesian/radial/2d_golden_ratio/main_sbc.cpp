// Gadgetron includes
#include "hoNDArray_fileio.h"
#include "cuNDArray_elemwise.h"
#include "cuNDArray_utils.h"
#include "cuNDArray_reductions.h"
#include "radial_utilities.h"
#include "cuNonCartesianSenseOperator.h"
#include "cuSenseBuffer.h"
#include "cuCgPreconditioner.h"
#include "cuPartialDerivativeOperator.h"
#include "cuCgSolver.h"
#include "cuSbcCgSolver.h"
#include "b1_map.h"
#include "GPUTimer.h"
#include "parameterparser.h"

// Std includes
#include <iostream>

using namespace std;
using namespace Gadgetron;
// Define desired precision
typedef float _real; 
typedef complext<_real> _complext;
typedef reald<_real,2>::Type _reald2;

// Upload samples for one reconstruction from host to device
boost::shared_ptr< cuNDArray<_complext> > 
upload_data( unsigned int reconstruction, unsigned int samples_per_reconstruction, unsigned int total_samples_per_coil, unsigned int num_coils, hoNDArray<_complext> *host_data )
{
  vector<size_t> dims; dims.push_back(samples_per_reconstruction); dims.push_back(num_coils);
  cuNDArray<_complext> *data = new cuNDArray<_complext>(); data->create( &dims );
  for( unsigned int i=0; i<num_coils; i++ )
    cudaMemcpy( data->get_data_ptr()+i*samples_per_reconstruction, 
		host_data->get_data_ptr()+i*total_samples_per_coil+reconstruction*samples_per_reconstruction, 
		samples_per_reconstruction*sizeof(_complext), cudaMemcpyHostToDevice );

  return boost::shared_ptr< cuNDArray<_complext> >(data);
}

int main(int argc, char** argv)
{
  //
  // Parse command line
  //
 
  ParameterParser parms;
  parms.add_parameter( 'd', COMMAND_LINE_STRING, 1, "Sample data file name", true );
  parms.add_parameter( 'r', COMMAND_LINE_STRING, 1, "Result file name", true, "result.cplx" );
  parms.add_parameter( 'm', COMMAND_LINE_INT,    1, "Matrix size", true );
  parms.add_parameter( 'o', COMMAND_LINE_INT,    1, "Oversampled matrix size", true );
  parms.add_parameter( 'p', COMMAND_LINE_INT,    1, "Profiles per frame", true );
  parms.add_parameter( 'f', COMMAND_LINE_INT,    1, "Frames per reconstruction (negative meaning all)", true, "-1" );
  parms.add_parameter( 'i', COMMAND_LINE_INT,    1, "Number of cg iterations", true, "10" );
  parms.add_parameter( 'I', COMMAND_LINE_INT,    1, "Number of sb inner iterations", true, "1" );
  parms.add_parameter( 'O', COMMAND_LINE_INT,    1, "Number of sb outer iterations", true, "10" );
  parms.add_parameter( 'k', COMMAND_LINE_FLOAT,  1, "Kernel width", true, "5.5" );
  parms.add_parameter( 'M', COMMAND_LINE_FLOAT,  1, "Mu", true, "1.0" );
  parms.add_parameter( 'L', COMMAND_LINE_FLOAT,  1, "Lambda", true, "2.0" );
  parms.add_parameter( 'A', COMMAND_LINE_FLOAT,  1, "Alpha in [0;1] (for PICCS)", true, "0.5" );

  parms.parse_parameter_list(argc, argv);
  if( parms.all_required_parameters_set() ){
    GINFO_STREAM(" Running reconstruction with the following parameters: " << endl);
    parms.print_parameter_list();
  }
  else{
    GINFO_STREAM(" Some required parameters are missing: " << endl);
    parms.print_parameter_list();
    parms.print_usage();
    return 1;
  }
  
  GPUTimer *timer;
  
  // Load sample data from disk
  timer = new GPUTimer("\nLoading data");
  boost::shared_ptr< hoNDArray<_complext> > host_data = read_nd_array<_complext>((char*)parms.get_parameter('d')->get_string_value());
  delete timer;
   
  if( !(host_data->get_number_of_dimensions() == 3) ){
    GINFO_STREAM(endl << "Input data is not three-dimensional (#samples/profile x #profiles x #coils). Quitting!\n" << endl);
    return 1;
  }

  // Configuration from the host data
  unsigned int samples_per_profile = host_data->get_size(0);
  unsigned int num_profiles = host_data->get_size(1);
  unsigned int num_coils = host_data->get_size(2);
  
  // Configuration from the command line
  uint64d2 matrix_size = uint64d2(parms.get_parameter('m')->get_int_value(), parms.get_parameter('m')->get_int_value());
  uint64d2 matrix_size_os = uint64d2(parms.get_parameter('o')->get_int_value(), parms.get_parameter('o')->get_int_value());
  _real kernel_width = parms.get_parameter('k')->get_float_value();
  unsigned int num_cg_iterations = parms.get_parameter('i')->get_int_value();
  unsigned int num_sb_inner_iterations = parms.get_parameter('I')->get_int_value();
  unsigned int num_sb_outer_iterations = parms.get_parameter('O')->get_int_value();
  unsigned int profiles_per_frame = parms.get_parameter('p')->get_int_value();
  unsigned int frames_per_reconstruction = parms.get_parameter('f')->get_int_value();

  _real mu = (_real) parms.get_parameter('M')->get_float_value();
  _real lambda = (_real) parms.get_parameter('L')->get_float_value();
  _real alpha = (_real) parms.get_parameter('A')->get_float_value();

  if( alpha>1 ) alpha = 1;
  if( alpha<0 ) alpha = 0;

  // Silent correction of invalid command line parameters (clamp to valid range)
  if( profiles_per_frame > num_profiles ) profiles_per_frame = num_profiles;
  if( frames_per_reconstruction < 0 ) frames_per_reconstruction = num_profiles / profiles_per_frame;
  if( frames_per_reconstruction*profiles_per_frame > num_profiles ) frames_per_reconstruction = num_profiles / profiles_per_frame;
  
  unsigned int profiles_per_reconstruction = frames_per_reconstruction*profiles_per_frame;
  unsigned int samples_per_frame = profiles_per_frame*samples_per_profile;
  unsigned int samples_per_reconstruction = profiles_per_reconstruction*samples_per_profile;

  GINFO_STREAM(endl << "#samples/profile: " << samples_per_profile);
  GINFO_STREAM(endl << "#profiles/frame: " << profiles_per_frame);
  GINFO_STREAM(endl << "#profiles: " << num_profiles);
  GINFO_STREAM(endl << "#coils: " << num_coils);
  GINFO_STREAM(endl << "#frames/reconstruction " << frames_per_reconstruction);
  GINFO_STREAM(endl << "#profiles/reconstruction " << profiles_per_reconstruction);
  GINFO_STREAM(endl << "#samples/reconstruction " << samples_per_reconstruction << endl << endl);

  // Density compensation weights are constant throughout all reconstrutions
  boost::shared_ptr< cuNDArray<_real> > dcw = compute_radial_dcw_golden_ratio_2d
    ( samples_per_profile, profiles_per_frame, (_real)matrix_size_os[0]/(_real)matrix_size[0], 
      _real(1)/((_real)samples_per_profile/(_real)max(matrix_size[0],matrix_size[1])) );
  
  // Define encoding matrix for non-Cartesian SENSE
  boost::shared_ptr< cuNonCartesianSenseOperator<_real,2> > E( new cuNonCartesianSenseOperator<_real,2>() );  
  E->set_weight( mu );
  E->setup( matrix_size, matrix_size_os, kernel_width );


  // Define rhs buffer
  //

  boost::shared_ptr< cuSenseBuffer<_real,2> > rhs_buffer( new cuSenseBuffer<_real,2>() );

  rhs_buffer->setup( matrix_size, matrix_size_os, kernel_width, num_coils, 8, 16 );
  rhs_buffer->set_dcw(dcw);

  //
  // Compute CSM using accumulation in the rhs buffer
  // 
 
  timer = new GPUTimer("CSM and regularization estimation");
    
  // Go through all the data...
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
    
  // Estimate csm
  boost::shared_ptr< cuNDArray<_complext> > acc_images = rhs_buffer->get_accumulated_coil_images();
  *acc_images *= rhs_buffer->get_normalization_factor();
  boost::shared_ptr< cuNDArray<_complext> > csm = boost::make_shared<cuNDArray<_complext>>(estimate_b1_map<_real,2>( acc_images.get() ));
  E->set_csm(csm);

  std::vector<size_t> reg_dims = to_std_vector(matrix_size);
  cuNDArray<_complext> _reg_image = cuNDArray<_complext>(&reg_dims);
  E->mult_csm_conj_sum( acc_images.get(), &_reg_image );

  // Duplicate the regularization image to 'frames_per_reconstruction' frames
  auto reg_image = boost::make_shared<cuNDArray<_complext>>(expand( _reg_image, frames_per_reconstruction ));

  acc_images.reset();

  // Define preconditioning weights
  //

  boost::shared_ptr< cuNDArray<_real> > _precon_weights = sum(abs_square(csm.get()).get(),2);
  reciprocal_sqrt_inplace(_precon_weights.get());
  boost::shared_ptr< cuNDArray<_complext> > precon_weights = real_to_complex<_complext>( _precon_weights.get() );
  _precon_weights.reset();

  // Define preconditioning matrix
  boost::shared_ptr< cuCgPreconditioner<_complext> > D( new cuCgPreconditioner<_complext>() );
  D->set_weights( precon_weights );
  precon_weights.reset();
  csm.reset();

  boost::shared_ptr< std::vector<size_t> > recon_dims( new std::vector<size_t> );
  *recon_dims = to_std_vector(matrix_size); recon_dims->push_back(frames_per_reconstruction); 

  // Define regularization operators 
  // We need "a pair" for PICCS
  //

  boost::shared_ptr< cuPartialDerivativeOperator<_complext,3> >
    Rx( new cuPartialDerivativeOperator<_complext,3>(0) );
  Rx->set_weight( (1.0f-alpha)*lambda );
  Rx->set_domain_dimensions(recon_dims.get());
  Rx->set_codomain_dimensions(recon_dims.get());

  boost::shared_ptr< cuPartialDerivativeOperator<_complext,3> >
    Ry( new cuPartialDerivativeOperator<_complext,3>(1) );
  Ry->set_weight( (1.0f-alpha)*lambda );
  Ry->set_domain_dimensions(recon_dims.get());
  Ry->set_codomain_dimensions(recon_dims.get());
 
  boost::shared_ptr< cuPartialDerivativeOperator<_complext,3> >
    Rz( new cuPartialDerivativeOperator<_complext,3>(2) );
  Rz->set_weight( (1.0f-alpha)*lambda );
  Rz->set_domain_dimensions(recon_dims.get());
  Rz->set_codomain_dimensions(recon_dims.get());

  boost::shared_ptr< cuPartialDerivativeOperator<_complext,3> >
    Rx2( new cuPartialDerivativeOperator<_complext,3>(0) );
  Rx2->set_weight( alpha*lambda );
  Rx2->set_domain_dimensions(recon_dims.get());
  Rx2->set_codomain_dimensions(recon_dims.get());

  boost::shared_ptr< cuPartialDerivativeOperator<_complext,3> >
    Ry2( new cuPartialDerivativeOperator<_complext,3>(1) );
  Ry2->set_weight( alpha*lambda );
  Ry2->set_domain_dimensions(recon_dims.get());
  Ry2->set_codomain_dimensions(recon_dims.get());
 
  boost::shared_ptr< cuPartialDerivativeOperator<_complext,3> >
    Rz2( new cuPartialDerivativeOperator<_complext,3>(2) );
  Rz2->set_weight( alpha*lambda );
  Rz2->set_domain_dimensions(recon_dims.get());
  Rz2->set_codomain_dimensions(recon_dims.get());

  delete timer;
    
  // 
  // Setup radial SENSE reconstructions
  //

  vector<size_t> data_dims; 
  data_dims.push_back(samples_per_reconstruction); data_dims.push_back(num_coils);

  E->set_domain_dimensions(recon_dims.get());
  E->set_codomain_dimensions(&data_dims);

  sqrt_inplace(dcw.get());
  E->set_dcw(dcw);
  // Setup split-Bregman solver
  cuSbcCgSolver<_complext> sb;
  sb.set_encoding_operator( E );
  
  // Add "TV" regularization
  if( alpha<1.0 ){
    sb.add_regularization_group_operator( Rx ); 
    sb.add_regularization_group_operator( Ry ); 
    sb.add_regularization_group_operator( Rz ); 
    sb.add_group();
  }
  
  // Add "PICCS" regularization
  if( alpha > 0.0 ){
    sb.add_regularization_group_operator( Rx2 ); 
    sb.add_regularization_group_operator( Ry2 ); 
    sb.add_regularization_group_operator( Rz2 ); 
    sb.add_group(reg_image);
  }
  
  sb.set_max_outer_iterations(num_sb_outer_iterations);
  sb.set_max_inner_iterations(num_sb_inner_iterations);
  sb.set_output_mode( cuSbcCgSolver<_complext>::OUTPUT_VERBOSE );

  sb.get_inner_solver()->set_preconditioner ( D );
  sb.get_inner_solver()->set_max_iterations( num_cg_iterations );
  sb.get_inner_solver()->set_tc_tolerance( 1e-4 );
  sb.get_inner_solver()->set_output_mode( cuCgSolver<_complext>::OUTPUT_WARNINGS );
  
  unsigned int num_reconstructions = num_profiles / profiles_per_reconstruction;

  // Allocate space for result
  std::vector<size_t> res_dims = to_std_vector(matrix_size); 
  res_dims.push_back(frames_per_reconstruction*num_reconstructions); 
  cuNDArray<_complext> result = cuNDArray<_complext>(&res_dims);

  timer = new GPUTimer("Full SENSE reconstruction with TV regularization.");

  for( unsigned int reconstruction = 0; reconstruction<num_reconstructions; reconstruction++ ){

    // Determine trajectories
    boost::shared_ptr< cuNDArray<_reald2> > traj = compute_radial_trajectory_golden_ratio_2d<_real>
      ( samples_per_profile, profiles_per_frame, frames_per_reconstruction, reconstruction*profiles_per_reconstruction );
    
    // Upload data
    boost::shared_ptr< cuNDArray<_complext> > data = upload_data
      ( reconstruction, samples_per_reconstruction, num_profiles*samples_per_profile, num_coils, host_data.get() );
    
    // Set current trajectory and trigger NFFT preprocessing
    E->preprocess(traj.get());
        
    *data *= *dcw;
    //
    // Split-Bregman solver
    //

    boost::shared_ptr< cuNDArray<_complext> > sbresult;
    {
      GPUTimer timer("GPU constrained Split Bregman solve");
      sbresult = sb.solve(data.get());
    }

    vector<size_t> tmp_dims = to_std_vector(matrix_size); tmp_dims.push_back(frames_per_reconstruction);
    cuNDArray<_complext> tmp(&tmp_dims, result.get_data_ptr()+reconstruction*prod(matrix_size)*frames_per_reconstruction );

    // Copy sbresult to result (pointed to by tmp)
    tmp = *sbresult;
  }
  
  delete timer;

  // All done, write out the result

  timer = new GPUTimer("Writing out result");

  boost::shared_ptr< hoNDArray<_complext> > host_result = result.to_host();
  write_nd_array<_complext>(host_result.get(), (char*)parms.get_parameter('r')->get_string_value());
    
  boost::shared_ptr< hoNDArray<_real> > host_norm = abs(&result)->to_host();
  write_nd_array<_real>( host_norm.get(), "result.real" );
  
  delete timer;

  return 0;
}
