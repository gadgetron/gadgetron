/*
  
  Sample application of the NFFT toolbox: using the NFFT matrix operator in a conjugate gradient solver
  
*/

#include "cuNFFT.h"
#include "radial_utilities.h"
#include "hoNDArray_fileio.h"
#include "cuNDArray.h"
#include "vector_td_utilities.h"
#include "parameterparser.h"
#include "../../../../../../toolboxes/nfft/NFFTOperator.h"
#include "cuCgSolver.h"
#include "GPUTimer.h"

#include <iostream>

using namespace std;
using namespace Gadgetron;

// Define desired precision
typedef float _real; 
typedef complext<_real> _complext;
typedef reald<_real,2>::Type _reald2;

int main( int argc, char** argv) 
{

  //
  // Parse command line
  //

  ParameterParser parms;
  parms.add_parameter( 'd', COMMAND_LINE_STRING, 1, "Input samples file name (.cplx)", true );
  parms.add_parameter( 'r', COMMAND_LINE_STRING, 1, "Output image file name (.cplx)", true, "result.cplx" );
  parms.add_parameter( 'm', COMMAND_LINE_INT,    1, "Matrix size", true );
  parms.add_parameter( 'o', COMMAND_LINE_INT,    1, "Oversampled matrix size", true );
  parms.add_parameter( 'k', COMMAND_LINE_FLOAT,  1, "Kernel width", true, "5.5" );
  parms.add_parameter( 'i', COMMAND_LINE_INT,    1, "Number of iterations", true, "10" );

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
  timer = new GPUTimer("Loading samples from disk");
  boost::shared_ptr< hoNDArray<_complext> > host_samples = read_nd_array<_complext>((char*)parms.get_parameter('d')->get_string_value());
  delete timer;
   
  if( !(host_samples->get_number_of_dimensions() == 2) ){
    cout << endl << "Samples ndarray is not two-dimensional (samples/profile x #profiles). Quitting.\n" << endl;
    return 1;
  }
  
  // Configuration from the command line
  uint64d2 matrix_size = uint64d2(parms.get_parameter('m')->get_int_value(), parms.get_parameter('m')->get_int_value());
  uint64d2 matrix_size_os = uint64d2(parms.get_parameter('o')->get_int_value(), parms.get_parameter('o')->get_int_value());
  _real kernel_width = parms.get_parameter('k')->get_float_value();
  unsigned int num_iterations = parms.get_parameter('i')->get_int_value();

  unsigned int num_profiles = host_samples->get_size(1);
  unsigned int samples_per_profile = host_samples->get_size(0);  
  _real alpha = (_real)matrix_size_os.vec[0]/(_real)matrix_size.vec[0];

  // Upload host data to device
  timer = new GPUTimer("Uploading samples to device");
  cuNDArray<_complext> samples(host_samples.get());
  delete timer;
  
  // Compute trajectories
  timer = new GPUTimer("Computing golden ratio radial trajectories");
  boost::shared_ptr< cuNDArray<_reald2> > traj = compute_radial_trajectory_golden_ratio_2d<_real>( samples_per_profile, num_profiles,  1 );
  delete timer;

  // Compute density compensation weights
  timer = new GPUTimer("Computing density compensation weights");
  boost::shared_ptr< cuNDArray<_real> > dcw = compute_radial_dcw_golden_ratio_2d
    ( samples_per_profile, num_profiles, alpha, _real(1)/((_real)samples_per_profile/(_real)matrix_size.vec[0]) );
  sqrt_inplace(dcw.get());
  delete timer;

  samples *= *dcw;

  // Define and setup NFFT encoding operator
  auto  E = boost::make_shared<NFFTOperator<cuNDArray,_real,2>>();
  
  E->setup( matrix_size, matrix_size_os, kernel_width );

  // Notify encoding operator of dcw
  E->set_dcw(dcw);
  
  // Set image dimensions
  vector<size_t> image_dims = to_std_vector(matrix_size);
  E->set_domain_dimensions(&image_dims);
  
  // Preprocess
  timer = new GPUTimer("NFFT preprocessing");
  E->preprocess( traj.get() );
  delete timer;

  // Setup conjugate gradient solver
  cuCgSolver< _complext> cg;
  cg.set_max_iterations( num_iterations );
  cg.set_tc_tolerance( 1e-6 );
  cg.set_output_mode( cuCgSolver<_complext>::OUTPUT_VERBOSE );
  cg.set_encoding_operator( E); 

  // Solve
  boost::shared_ptr< cuNDArray<_complext> > cgresult;
  {
    GPUTimer timer("GPU Conjugate Gradient solve");
    cgresult = cg.solve(&samples);
  }
  
  //
  // Output result
  //
  
  timer = new GPUTimer("Output result to disk");

  boost::shared_ptr< hoNDArray<_complext> > host_image = cgresult->to_host();
  write_nd_array<_complext>( host_image.get(), (char*)parms.get_parameter('r')->get_string_value());

  boost::shared_ptr< hoNDArray<_real> > host_norm = abs(cgresult.get())->to_host();
  write_nd_array<_real>( host_norm.get(), "result.real" );

  delete timer;

  return 0;
}
