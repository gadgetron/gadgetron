// Gadgetron includes
#include "cuNDArray_elemwise.h"
#include "cuNDArray_utils.h"
#include "cuNDArray_reductions.h"
#include "hoNDArray_fileio.h"
#include "vector_td_utilities.h"
#include "cuImageOperator.h"
#include "radial_utilities.h"
#include "cuNonCartesianSenseOperator.h"
#include "cuSenseBuffer.h"
#include "cuCgPreconditioner.h"
#include "cuCgSolver.h"
#include "b1_map.h"
#include "parameterparser.h"
#include "GPUTimer.h"

// Std includes
#include <iostream>

using namespace std;
using namespace Gadgetron;

// Define desired precision
typedef float _real; 
typedef complext<_real> _complext;
typedef reald<_real,2>::Type _reald2;

const bool use_atomics = false;

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
  parms.add_parameter( 'i', COMMAND_LINE_INT,    1, "Number of iterations", true, "10" );
  parms.add_parameter( 'k', COMMAND_LINE_FLOAT,  1, "Kernel width", true, "5.5" );
  parms.add_parameter( 'K', COMMAND_LINE_FLOAT,  1, "Kappa", true, "0.3" );

  parms.parse_parameter_list(argc, argv);
  if( parms.all_required_parameters_set() ){
    cout << " Running reconstruction with the following parameters: " << endl;
    parms.print_parameter_list();
  }
  else{
    cout << " Some required parameters are missing: " << endl;
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
    cout << endl << "Input data is not three-dimensional (#samples/profile x #profiles x #coils). Quitting!\n" << endl;
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
  _real kappa = parms.get_parameter('K')->get_float_value();
  unsigned int num_iterations = parms.get_parameter('i')->get_int_value();
  unsigned int profiles_per_frame = parms.get_parameter('p')->get_int_value();
  unsigned int frames_per_reconstruction = parms.get_parameter('f')->get_int_value();

  // Silent correction of invalid command line parameters (clamp to valid range)
  if( profiles_per_frame > num_profiles ) profiles_per_frame = num_profiles;
  if( frames_per_reconstruction < 0 ) frames_per_reconstruction = num_profiles / profiles_per_frame;
  if( frames_per_reconstruction*profiles_per_frame > num_profiles ) frames_per_reconstruction = num_profiles / profiles_per_frame;
  
  unsigned int profiles_per_reconstruction = frames_per_reconstruction*profiles_per_frame;
  unsigned int samples_per_frame = profiles_per_frame*samples_per_profile;
  unsigned int samples_per_reconstruction = profiles_per_reconstruction*samples_per_profile;

  cout << endl << "#samples/profile: " << samples_per_profile;
  cout << endl << "#profiles/frame: " << profiles_per_frame;
  cout << endl << "#profiles: " << num_profiles;
  cout << endl << "#coils: " << num_coils;
  cout << endl << "#frames/reconstruction: " << frames_per_reconstruction;
  cout << endl << "#profiles/reconstruction: " << profiles_per_reconstruction;
  cout << endl << "#samples/reconstruction: " << samples_per_reconstruction << endl << endl;

  // Set density compensation weights
  boost::shared_ptr< cuNDArray<_real> > dcw = compute_radial_dcw_golden_ratio_2d
    ( samples_per_profile, profiles_per_frame, (_real)matrix_size_os[0]/(_real)matrix_size[0], 
      _real(1)/((_real)samples_per_profile/(_real)max(matrix_size[0],matrix_size[1])) );

  // Define encoding matrix for non-Cartesian SENSE
  boost::shared_ptr< cuNonCartesianSenseOperator<_real,2> > E
    ( new cuNonCartesianSenseOperator<_real,2>() );

  E->setup( matrix_size, matrix_size_os, kernel_width );



  // Define rhs buffer
  //

  boost::shared_ptr< cuSenseBuffer<_real,2> > rhs_buffer
    ( new cuSenseBuffer<_real,2>() );

  rhs_buffer->setup( matrix_size, matrix_size_os, kernel_width, num_coils, 8, 16 );
  rhs_buffer->set_dcw(dcw);

  // Fill rhs buffer
  //

  timer = new GPUTimer("Filling rhs buffer");
    
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

  delete timer;


  // Estimate CSM
  //

  timer = new GPUTimer("Estimating csm");

  boost::shared_ptr< cuNDArray<_complext> > acc_images = rhs_buffer->get_accumulated_coil_images();
  boost::shared_ptr< cuNDArray<_complext> > csm = boost::make_shared<cuNDArray<_complext>>(estimate_b1_map<_real,2>( *acc_images.get() ));
  E->set_csm(csm);

  delete timer;
  

  // Define regularization image operator 
  //

  timer = new GPUTimer("Computing regularization");

  std::vector<size_t> image_dims = to_std_vector(matrix_size);
  cuNDArray<_complext> reg_image = cuNDArray<_complext>(&image_dims);

  E->mult_csm_conj_sum( acc_images.get(), &reg_image );
  acc_images.reset();

  boost::shared_ptr< cuImageOperator<_complext> > R( new cuImageOperator<_complext>() );
  R->set_weight( kappa );
  R->compute( &reg_image );

  delete timer;

  // Define preconditioning weights
  //

  timer = new GPUTimer("Computing preconditioning weights");

  boost::shared_ptr< cuNDArray<_real> > _precon_weights = sum(abs_square(csm.get()).get(),2);
  boost::shared_ptr< cuNDArray<_real> > R_diag = R->get();
  *R_diag *= kappa;
  *_precon_weights += *R_diag;
  R_diag.reset();
  reciprocal_sqrt_inplace(_precon_weights.get());
  boost::shared_ptr< cuNDArray<_complext> > precon_weights = real_to_complex<_complext>( _precon_weights.get() );
  _precon_weights.reset();

  // Define preconditioning matrix
  boost::shared_ptr< cuCgPreconditioner<_complext> > D( new cuCgPreconditioner<_complext>() );
  D->set_weights( precon_weights );
  precon_weights.reset();
  csm.reset();

  delete timer;
  
  // 
  // Setup radial SENSE reconstructions
  //
  // Notify encoding operator of dcw
  sqrt_inplace(dcw.get());
	E->set_dcw(dcw);
  // Setup conjugate gradient solver
  cuCgSolver<_complext> cg;
  cg.set_preconditioner ( D );  // preconditioning matrix
  cg.set_max_iterations( num_iterations );
  cg.set_tc_tolerance( 1e-6 );
  cg.set_output_mode( cuCgSolver< _complext>::OUTPUT_VERBOSE );
  cg.set_encoding_operator( E );        // encoding matrix
  cg.add_regularization_operator( R );  // regularization matrix
  
  // Reconstruct all SENSE frames iteratively
  unsigned int num_reconstructions = num_profiles / profiles_per_reconstruction;
  
  // Allocate space for result
  image_dims.push_back(frames_per_reconstruction*num_reconstructions); 
  cuNDArray<_complext> result = cuNDArray<_complext>(&image_dims);
  
  timer = new GPUTimer("Full SENSE reconstruction.");
  
  // Define image dimensions
  image_dims = to_std_vector(matrix_size); 
  image_dims.push_back(frames_per_reconstruction);
  
  for( unsigned int reconstruction = 0; reconstruction<num_reconstructions; reconstruction++ ){

    // Determine trajectories
    boost::shared_ptr< cuNDArray<_reald2> > traj = compute_radial_trajectory_golden_ratio_2d<_real>
      ( samples_per_profile, profiles_per_frame, frames_per_reconstruction, reconstruction*profiles_per_reconstruction );
    
    // Upload data
    boost::shared_ptr< cuNDArray<_complext> > data = upload_data
      ( reconstruction, samples_per_reconstruction, num_profiles*samples_per_profile, num_coils, host_data.get() );
    
    // Pass image dimensions to encoding operator
    E->set_domain_dimensions(&image_dims);
    E->set_codomain_dimensions(data->get_dimensions().get());
  
    // Set current trajectory and trigger NFFT preprocessing
    E->preprocess(traj.get());
    
    *data *= *dcw;
    //
    // Invoke conjugate gradient solver
    //

    boost::shared_ptr< cuNDArray<_complext> > cgresult;
    {
      GPUTimer timer("GPU Conjugate Gradient solve");
      cgresult = cg.solve(data.get());
    }

    if( !cgresult.get() )
      return 1;

    // Copy cgresult to overall result
    cuNDArray<_complext> out(&image_dims, result.get_data_ptr()+reconstruction*prod(matrix_size)*frames_per_reconstruction );    
    out = *(cgresult.get());
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
