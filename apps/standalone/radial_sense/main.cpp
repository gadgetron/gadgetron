// Gadgetron includes
#include "cuNDArray.h"
#include "hoNDArray_fileio.h"
#include "vector_td_operators.h"
#include "vector_td_utilities.h"
#include "ndarray_vector_td_utilities.h"
#include "radial_utilities.h"
#include "cgOperatorNonCartesianSense.h"
#include "vector_td_identity_operator.h"
#include "cuCG.h"
#include "b1_map.h"
#include "GPUTimer.h"
#include "parameterparser.h"

// Std includes
#include <iostream>

using namespace std;

// Upload samples for one reconstruction from host to device
auto_ptr< cuNDArray<float_complext::Type> > 
upload_data( unsigned int reconstruction, unsigned int samples_per_reconstruction, unsigned int total_samples_per_coil, unsigned int num_coils,
	     hoNDArray<float_complext::Type> *host_data )
{
  vector<unsigned int> dims; dims.push_back(samples_per_reconstruction); dims.push_back(num_coils);
  cuNDArray<float_complext::Type> *data = new cuNDArray<float_complext::Type>(); data->create( dims );
  for( unsigned int i=0; i<num_coils; i++ )
    cudaMemcpy( data->get_data_ptr()+i*samples_per_reconstruction, 
		host_data->get_data_ptr()+i*total_samples_per_coil+reconstruction*samples_per_reconstruction, 
		samples_per_reconstruction*sizeof(float_complext::Type), cudaMemcpyHostToDevice );

  return auto_ptr< cuNDArray<float_complext::Type> >(data);
}

int main(int argc, char** argv)
{

  //
  // Parse command line
  //

  ParameterParser parms;
  parms.add_parameter( 'd', COMMAND_LINE_STRING, 1, "Sample data file name", true );
  parms.add_parameter( 'm', COMMAND_LINE_INT,    1, "Matrix size", true );
  parms.add_parameter( 'o', COMMAND_LINE_INT,    1, "Oversampled Matrix size", true );
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
  hoNDArray<float_complext::Type> host_data = read_nd_array<float_complext::Type>((char*)parms.get_parameter('d')->get_string_value());
  delete timer;
   
  if( !host_data.get_number_of_dimensions() == 3 ){
    cout << endl << "Input data is not three-dimensional (#samples/profile x #profiles x #coils). Quitting!\n" << endl;
    return 1;
  }

  // Configuration from the host data
  unsigned int samples_per_profile = host_data.get_size(0);
  unsigned int num_profiles = host_data.get_size(1);
  unsigned int num_coils = host_data.get_size(2);
  
  // Configuration from the command line
  uintd2 matrix_size = uintd2(parms.get_parameter('m')->get_int_value(), parms.get_parameter('m')->get_int_value());
  uintd2 matrix_size_os = uintd2(parms.get_parameter('o')->get_int_value(), parms.get_parameter('o')->get_int_value());
  float kernel_width = parms.get_parameter('k')->get_float_value();
  float kappa = parms.get_parameter('K')->get_float_value();
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
  auto_ptr< cuNDArray<float> > _dcw = compute_radial_dcw_golden_ratio_2d
    ( samples_per_profile, profiles_per_frame, (float)matrix_size_os.vec[0]/(float)matrix_size.vec[0], 
      get_one<float>()/((float)samples_per_profile/(float)max(matrix_size.vec[0],matrix_size.vec[1])) );
  vector<unsigned int> dcw_dims; dcw_dims.push_back(samples_per_frame); dcw_dims.push_back(frames_per_reconstruction);
  cuNDArray<float> dcw; dcw.create(dcw_dims);
  for(unsigned int i=0; i<frames_per_reconstruction; i++ ){
    cudaMemcpy( dcw.get_data_ptr()+i*samples_per_frame, _dcw->get_data_ptr(), samples_per_frame*sizeof(float), cudaMemcpyDeviceToDevice );
  }
  
  //
  // Compute CSM, preconditioner, and regularization image
  // 

  auto_ptr< cuNDArray<float_complext::Type> > csm;
  auto_ptr< cuNDArray<float_complext::Type> > preconditioner;
  //  auto_ptr< cuNDArray<float_complext::Type> > regularization;

  {
    timer = new GPUTimer("CSM estimation and preconditioning computation");
    
    vector<unsigned int> image_os_dims = uintd_to_vector<2>(matrix_size_os); 
    image_os_dims.push_back(frames_per_reconstruction); image_os_dims.push_back(num_coils);    
    cuNDArray<float_complext::Type> *image_os = new cuNDArray<float_complext::Type>(); image_os->create(image_os_dims);
    
    NFFT_plan<float, 2> plan( matrix_size, matrix_size_os, kernel_width );

    for( unsigned int iteration = 0; iteration < num_profiles/profiles_per_reconstruction; iteration++ ) {
  
      // Define trajectories
      auto_ptr< cuNDArray<floatd2::Type> > traj = compute_radial_trajectory_golden_ratio_2d<float>
	( samples_per_profile, profiles_per_frame, frames_per_reconstruction, iteration*profiles_per_reconstruction );

      // Preprocess
      plan.preprocess( traj.get(), NFFT_plan<float,2>::NFFT_PREP_BACKWARDS );
      traj.reset();
     
      // Upload data
      auto_ptr< cuNDArray<float_complext::Type> > csm_data = upload_data
	( iteration, samples_per_reconstruction, num_profiles*samples_per_profile, num_coils, &host_data );
      
      // Accumulate k-space for CSM estimation.
      plan.convolve( csm_data.get(), image_os, &dcw, NFFT_plan<float,2>::NFFT_BACKWARDS, (iteration==0) ? false : true );
      csm_data.reset();
    }
    
    // We now have 'frames_per_reconstruction' k-space images of each coil. Add these up.
    auto_ptr< cuNDArray<float_complext::Type> > acc_image_os = cuNDA_sum<float_complext::Type>( image_os, 2 );    
    delete image_os;

    // Complete gridding of k-space CSM image
    plan.fft( acc_image_os.get(), NFFT_plan<float,2>::NFFT_BACKWARDS );
    plan.deapodize( acc_image_os.get() );

    // Remove oversampling
    vector<unsigned int> image_dims = uintd_to_vector<2>(matrix_size); image_dims.push_back(num_coils);
    cuNDArray<float_complext::Type> *image = new cuNDArray<float_complext::Type>(); image->create(image_dims);
    cuNDA_crop<float_complext::Type,2>( (matrix_size_os-matrix_size)>>1, acc_image_os.get(), image );
    acc_image_os.reset();
    
    // Estimate CSM
    csm = estimate_b1_map<float,2>( image );
    delete image;
    
    // Define preconditioner
    preconditioner = cuNDA_reciprocal_rss<float_complext::Type>( csm.get(), 2 );

    delete timer;
  }
  
  //hoNDArray<float_complext::Type> out_csm = csm->to_host();
  //write_nd_array<float_complext::Type>(out_csm,"csm.cplx");

  // 
  // Setup radial SENSE reconstructions
  //

  // Define encoding matrix for non-Cartesian SENSE
  cgOperatorNonCartesianSense<float,2> E;
  
  E.setup( matrix_size, matrix_size_os, kernel_width );

  if( E.set_csm(csm.get()) < 0 ) {
    cout << "Failed to set csm on encoding matrix" << endl;
  }
  
  // Define preconditioning matrix
  cuCGPrecondWeight<float_complext::Type> D;
  D.set_weights( preconditioner.get() );
  
  // Define conjugate gradient solver
  cuCG<float, float_complext::Type> cg;

  // Define regularization matrix
  cuCGIdentityOperator<float_complext::Type> R( cg.get_cublas_handle() );

  // Setup solver
  cg.add_matrix_operator( &E, 1.0f );  // encoding matrix
  cg.add_matrix_operator( &R, kappa ); // regularization matrix
  cg.set_preconditioner ( &D );        // preconditioning matrix
  cg.set_iterations( num_iterations );
  cg.set_limit( 1e-5 );
  cg.set_output_mode( cuCG<float, float_complext::Type>::OUTPUT_VERBOSE );
  
  unsigned int reconstruction = 0;

  // Determine trajectories
  auto_ptr< cuNDArray<floatd2::Type> > traj = compute_radial_trajectory_golden_ratio_2d<float>
    ( samples_per_profile, profiles_per_frame, frames_per_reconstruction, reconstruction*profiles_per_reconstruction );
  
  // Upload data
  auto_ptr< cuNDArray<float_complext::Type> > data = upload_data
    ( reconstruction, samples_per_reconstruction, num_profiles*samples_per_profile, num_coils, &host_data );

  // Set current trajectory and trigger NFFT preprocessing
  if( E.set_trajectory(traj.get()) < 0 ) {
    cout << "Failed to set trajectory on encoding matrix" << endl;
  }

  // Associate density compensation weights
  if( E.set_weights(&dcw) < 0 ) {
    cout << "Failed to set density compensation weights on encoding matrix" << endl;
  }

  // Form rhs
  vector<unsigned int> rhs_dims = uintd_to_vector<2>(matrix_size);
  if( frames_per_reconstruction > 1 ) rhs_dims.push_back(frames_per_reconstruction);
  cuNDArray<float_complext::Type> rhs; rhs.create(rhs_dims);
  
  E.mult_MH( data.get(), &rhs );

  hoNDArray<float_complext::Type> out_rhs = rhs.to_host();
  write_nd_array<float_complext::Type>(out_rhs,"rhs.cplx");
  exit(1);

  auto_ptr< cuNDArray<float_complext::Type> > cgresult;
  {
    GPUTimer timer("GPU Conjugate Gradient solve");
    cgresult = cg.solve(&rhs);
  }

  hoNDArray<float_complext::Type> rho_out = cgresult->to_host();
  write_nd_array<float_complext::Type>(rho_out,"rho_out.cplx");

  hoNDArray<float> host_norm = cuNDA_norm<float,2>(cgresult.get())->to_host();
  write_nd_array<float>( host_norm, "rho_out.real" );

  cout << "Reconstruction done" << endl;
  
  return 0;
}
