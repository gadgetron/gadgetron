// Gadgetron includes
#include "cuNDArray.h"
#include "hoNDArray_fileio.h"
#include "ndarray_vector_td_utilities.h"
#include "radial_utilities.h"
#include "cuNonCartesianKtSenseOperator.h"
#include "cuSenseRHSBuffer.h"
#include "cuImageOperator.h"
#include "cuCGPrecondWeights.h"
#include "cuCGSolver.h"
#include "cuNDFFT.h"
#include "b1_map.h"
#include "GPUTimer.h"
#include "parameterparser.h"

// Std includes
#include <iostream>
#include <math.h>

using namespace std;

// Define desired precision
typedef float _real; 
typedef complext<_real> _complext;
typedef reald<_real,2>::Type _reald2;

// Upload samples for one reconstruction from host to device
boost::shared_ptr< cuNDArray<_complext> > 
upload_data( unsigned int reconstruction, unsigned int samples_per_reconstruction, unsigned int total_samples_per_coil, unsigned int num_coils,
	     hoNDArray<_complext> *host_data )
{
  vector<unsigned int> dims; dims.push_back(samples_per_reconstruction); dims.push_back(num_coils);
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
  parms.add_parameter( 'K', COMMAND_LINE_FLOAT,  1, "Kappa", true, "0.1" );

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
  uintd2 matrix_size = uintd2(parms.get_parameter('m')->get_int_value(), parms.get_parameter('m')->get_int_value());
  uintd2 matrix_size_os = uintd2(parms.get_parameter('o')->get_int_value(), parms.get_parameter('o')->get_int_value());
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
  cout << endl << "#frames/reconstruction " << frames_per_reconstruction;
  cout << endl << "#profiles/reconstruction " << profiles_per_reconstruction;
  cout << endl << "#samples/reconstruction " << samples_per_reconstruction << endl << endl;

  // Density compensation weights are constant throughout all reconstrutions
  boost::shared_ptr< cuNDArray<_real> > dcw = compute_radial_dcw_golden_ratio_2d
    ( samples_per_profile, profiles_per_frame, (_real)matrix_size_os.vec[0]/(_real)matrix_size.vec[0], 
      _real(1)/((_real)samples_per_profile/(_real)max(matrix_size.vec[0],matrix_size.vec[1])) );
  
  // Define encoding matrix for non-Cartesian kt-SENSE
  boost::shared_ptr< cuNonCartesianKtSenseOperator<_real,2> > E( new cuNonCartesianKtSenseOperator<_real,2>() );

  if( E->setup( matrix_size, matrix_size_os, kernel_width ) < 0 ){
    cout << "Failed to setup non-Cartesian Sense operator" << endl;
    return 1;
  }

  // Notify encoding operator of dcw
  if( E->set_dcw(dcw) < 0 ) {
    cout << "Failed to set density compensation weights on encoding matrix" << endl;
    return 1;
  }

  // Use a rhs buffer to estimate the csm
  //

  boost::shared_ptr< cuSenseRHSBuffer<_real,2> > rhs_buffer( new cuSenseRHSBuffer<_real,2>() );

  rhs_buffer->set_num_coils(num_coils);

  if( rhs_buffer->set_sense_operator(E) < 0 ){
    cout << "Failed to set sense operator on rhs buffer" << endl;
  }
   
  // Fill rhs buffer
  //
 
  timer = new GPUTimer("CSM estimation");
    
  // Go through all the data...
  for( unsigned int iteration = 0; iteration < num_profiles/profiles_per_frame; iteration++ ) {

    // Define trajectories
    boost::shared_ptr< cuNDArray<_reald2> > traj = compute_radial_trajectory_golden_ratio_2d<_real>
      ( samples_per_profile, profiles_per_frame, 1, iteration*profiles_per_reconstruction );
    
    // Upload data
    boost::shared_ptr< cuNDArray<_complext> > csm_data = upload_data
      ( iteration, samples_per_frame, num_profiles*samples_per_profile, num_coils, host_data.get() );
        
    // Add frame to rhs buffer
    rhs_buffer->add_frame_data( csm_data.get(), traj.get() );
  }
  
  boost::shared_ptr< cuNDArray<_complext> > acc_images = rhs_buffer->get_acc_coil_images();

  boost::shared_ptr< cuNDArray<_complext> > csm = estimate_b1_map<_real,2>( acc_images.get() );

  if( E->set_csm(csm) < 0 ) {
    cout << "Failed to set csm on encoding matrix" << endl;
    return 1;
  }

  acc_images.reset();
  rhs_buffer.reset();
 
  delete timer;

  // 
  // Setup radial kt-SENSE reconstructions
  //
    
  // Define regularization image operator
  boost::shared_ptr< cuImageOperator<_real,_complext> > R( new cuImageOperator<_real,_complext>() ); 
  R->set_weight( kappa );

  // Define preconditioning operator
  boost::shared_ptr< cuCGPrecondWeights<_complext> > D( new cuCGPrecondWeights<_complext>() );
  boost::shared_ptr< cuNDArray<_real> > ___precon_weights = cuNDA_ss<_real,_complext>( csm.get(), 2 ); 
  boost::shared_ptr< cuNDArray<_real> > __precon_weights = cuNDA_expand<_real>( ___precon_weights.get(), frames_per_reconstruction );
  ___precon_weights.reset();

  // Setup conjugate gradient solver
  cuCGSolver<_real, _complext> cg;
  cg.add_matrix_operator( E );  // encoding matrix
  cg.add_matrix_operator( R );  // regularization matrix
  cg.set_preconditioner ( D );  // preconditioning matrix
  cg.set_max_iterations( num_iterations );
  cg.set_limit( 1e-6 );
  cg.set_output_mode( cuCGSolver<_real, _complext>::OUTPUT_VERBOSE );
      
  // Reconstruct all SENSE frames iteratively
  unsigned int num_reconstructions = num_profiles / profiles_per_reconstruction;
  
  // Allocate space for result
  vector<unsigned int> image_dims = uintd_to_vector<2>(matrix_size); 
  image_dims.push_back(frames_per_reconstruction*num_reconstructions); 

  cuNDArray<_complext> result; 
  if( result.create(&image_dims) == 0x0 ){
    cout << "Failed to allocate result " << endl;
    return 1;
  }
  
  // Define shutter for training data
  _reald2 shutter_radius = to_vector_td<_real,2>(((_real)matrix_size_os.vec[0]/(_real)matrix_size.vec[0])*(_real)profiles_per_frame/(_real)M_PI);
  
  vector<unsigned int> image_os_dims = uintd_to_vector<2>(matrix_size_os); 
  image_os_dims.push_back(frames_per_reconstruction); image_os_dims.push_back(num_coils);    
  cuNDArray<_complext> *image_os = new cuNDArray<_complext>(); 

  if( image_os->create(&image_os_dims) == 0x0 ){
    cout << "Failed to allocate image_os " << endl;
    return 1;
  }

  timer = new GPUTimer("Full SENSE reconstruction.");
  
  for( unsigned int reconstruction = 0; reconstruction<num_reconstructions; reconstruction++ ){

    // 
    // Estimate training data
    // 

    // Define trajectories
    boost::shared_ptr< cuNDArray<_reald2> > traj = compute_radial_trajectory_golden_ratio_2d<_real>
      ( samples_per_profile, profiles_per_frame, frames_per_reconstruction, reconstruction*profiles_per_reconstruction );
    
    // Preprocess
    E->preprocess( traj.get() );
    
    // Upload data
    boost::shared_ptr< cuNDArray<_complext> > data = upload_data
      ( reconstruction, samples_per_reconstruction, num_profiles*samples_per_profile, num_coils, host_data.get() );
    
    // Convolve to Cartesian k-space
    E->get_plan()->convolve( data.get(), image_os, dcw.get(), NFFT_plan<_real,2>::NFFT_BACKWARDS );

    // Apply shutter
    cuNDA_zero_fill_border<_real,_complext,2>( shutter_radius, image_os );
    E->get_plan()->fft( image_os, NFFT_plan<_real,2>::NFFT_BACKWARDS );
    E->get_plan()->deapodize( image_os );

    // Remove oversampling
    image_dims = uintd_to_vector<2>(matrix_size);
    image_dims.push_back(frames_per_reconstruction); image_dims.push_back(num_coils);
    cuNDArray<_complext> *image = new cuNDArray<_complext>(); 

    if( image->create(&image_dims) == 0x0 ){
      cout << "Failed to allocate image " << endl;
      return 1;
    }
    
    cuNDA_crop<_complext,2>( (matrix_size_os-matrix_size)>>1, image_os, image );
    
    // Compute regularization image
    cuNDArray<_complext> *reg_image = new cuNDArray<_complext>(); 

    image_dims.pop_back();
    if( reg_image->create(&image_dims) == 0x0 ){
      cout << "Failed to allocate regularization image " << endl;
      return 1;
    }

    E->mult_csm_conj_sum( image, reg_image );
    cuNDFFT<_complext>().ifft( reg_image, 2, true );
    R->compute( reg_image );

    delete reg_image; reg_image = 0x0;
    delete image; image = 0x0;
    
    // Define preconditioning weights
    cuNDArray<_real> _precon_weights(*__precon_weights.get());
    cuNDA_axpy<_real>( kappa, R->get(), &_precon_weights );  
    cuNDA_reciprocal_sqrt<_real>( &_precon_weights );
    boost::shared_ptr< cuNDArray<_complext> > precon_weights = 
      cuNDA_real_to_complext<_real>( &_precon_weights );
    
    // Define preconditioning matrix
    D->set_weights( precon_weights );
    precon_weights.reset();
      
    // Form rhs (use result array to save memory)
    cuNDArray<_complext> rhs; 

    if( rhs.create(&image_dims, result.get_data_ptr()+reconstruction*prod(matrix_size)*frames_per_reconstruction ) == 0x0 ){
      cout << "Failed to create rhs array" << endl;
      return 1;
    }

    E->mult_MH( data.get(), &rhs );
    
    //
    // Conjugate gradient solver
    //

    boost::shared_ptr< cuNDArray<_complext> > cgresult;
    {
      GPUTimer timer("GPU Conjugate Gradient solve");
      cgresult = cg.solve_from_rhs(&rhs);
    }

    // Goto from x-f to x-t space
    cuNDFFT<_complext>().fft( cgresult.get(), 2 );
    
    // Copy cgresult to result (pointed to by rhs)
    rhs = *(cgresult.get());  
  }
  
  delete timer;
  delete image_os; image_os = 0x0;
  csm.reset();

  // All done, write out the result

  timer = new GPUTimer("Writing out result");
  
  boost::shared_ptr< hoNDArray<_complext> > host_result = result.to_host();
  write_nd_array<_complext>(host_result.get(), (char*)parms.get_parameter('r')->get_string_value());
    
  boost::shared_ptr< hoNDArray<_real> > host_norm = cuNDA_cAbs<_real,_complext>(&result)->to_host();
  write_nd_array<_real>( host_norm.get(), "result.real" );
  
  delete timer;

  return 0;
}
