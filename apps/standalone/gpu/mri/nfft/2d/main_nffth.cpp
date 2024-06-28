/*

  Sample application of the NFFT toolbox: standalone "gridding" example.

  -----------

  The nfft is written generically and templetized to
  - transform arbitrary trajectories
  - transform an arbitrary number of dimensions (currently instantiated for 1d/2d/3d/4d)
  - support both single and double precision

  General principles of the implementation can be found in:

  Accelerating the Non-equispaced Fast Fourier Transform on Commodity Graphics Hardware.
  T.S. Sørensen, T. Schaeffter, K.Ø. Noe, M.S. Hansen. 
  IEEE Transactions on Medical Imaging 2008; 27(4):538-547.

  Real-time Reconstruction of Sensitivity Encoded Radial Magnetic Resonance Imaging Using a Graphics Processing Unit.
  T.S. Sørensen, D. Atkinson, T. Schaeffter, M.S. Hansen.
  IEEE Transactions on Medical Imaging 2009; 28(12):1974-1985. 

  This example programme of the nnft utilizes golden ratio based radial trajectories 
  and outputs a gridded image from input ndarrays of the corresponding samples, trajectory, and density compensation weights.

*/

#include "cuNFFT.h"
#include "radial_utilities.h"
#include "vector_td_utilities.h"
#include "hoNDArray_fileio.h"
#include "cuNDArray_elemwise.h"
#include "GPUTimer.h"
#include "parameterparser.h"
#include "complext.h"

#include <iostream>

using namespace std;
using namespace Gadgetron;

// Define desired precision
typedef float _real; 
typedef complext<_real> _complext;
typedef reald<_real,2>::Type _reald2;
typedef cuNFFT_impl<_real,2> plan_type;

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
  timer = new GPUTimer("Loading samples from disk");
  boost::shared_ptr< hoNDArray<_complext> > host_samples = read_nd_array<_complext>((char*)parms.get_parameter('d')->get_string_value());
  delete timer;
   
  if( !(host_samples->get_number_of_dimensions() == 2) ){
    GINFO_STREAM(endl << "Samples ndarray is not two-dimensional (samples/profile x #profiles). Quitting.\n" << endl);
    return 1;
  }
  
  // Configuration from the command line
  uint64d2 matrix_size = uint64d2(parms.get_parameter('m')->get_int_value(), parms.get_parameter('m')->get_int_value());
  uint64d2 matrix_size_os = uint64d2(parms.get_parameter('o')->get_int_value(), parms.get_parameter('o')->get_int_value());
  _real kernel_width = parms.get_parameter('k')->get_float_value();

  unsigned int num_profiles = host_samples->get_size(1);
  unsigned int samples_per_profile = host_samples->get_size(0);  
  _real alpha = (_real)matrix_size_os.vec[0]/(_real)matrix_size.vec[0];

  // Upload host data to device
  timer = new GPUTimer("Uploading samples to device");
  cuNDArray<_complext> samples(*host_samples);
  delete timer;
  
  // Setup resulting image array
  vector<size_t> image_dims = to_std_vector(matrix_size);
  cuNDArray<_complext> image(image_dims);
  
  // Initialize plan
  timer = new GPUTimer("Initializing plan");
  plan_type plan( matrix_size, matrix_size_os, kernel_width );
  delete timer;

  // Compute trajectories
  timer = new GPUTimer("Computing golden ratio radial trajectories");
  boost::shared_ptr< cuNDArray<_reald2> > traj = compute_radial_trajectory_golden_ratio_2d<_real>( samples_per_profile, num_profiles,  1 );
  delete timer;
  
  // Preprocess
  timer = new GPUTimer("NFFT preprocessing");
  plan.preprocess( *traj, NFFT_prep_mode::NC2C );
  delete timer;

  // Compute density compensation weights
  timer = new GPUTimer("Computing density compensation weights");
  boost::shared_ptr< cuNDArray<_real> > dcw = compute_radial_dcw_golden_ratio_2d
    ( samples_per_profile, num_profiles, alpha, _real(1)/((_real)samples_per_profile/(_real)matrix_size.vec[0]) );
  delete timer;

  // Gridder
  timer = new GPUTimer("Computing adjoint nfft (gridding)");
  plan.compute( samples, image, dcw.get(), NFFT_comp_mode::BACKWARDS_NC2C );
  delete timer;

  //
  // Output result
  //
  
  timer = new GPUTimer("Output result to disk");

  boost::shared_ptr< hoNDArray<_complext> > host_image = image.to_host();
  write_nd_array<_complext>( host_image.get(), (char*)parms.get_parameter('r')->get_string_value());


  boost::shared_ptr< hoNDArray<_real> > host_norm = abs(&image)->to_host();
  write_nd_array<_real>( host_norm.get(), "result.real" );

  delete timer;

  return 0;
}
