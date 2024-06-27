/*

  Sample application of the NFFT toolbox: standalone "inverse gridding" example.

  -----------

  The nfft is written generically and templetized to

  - transform arbitrary trajectories
  - transform an "arbitrary" number of dimensions (currently instantiated for 1d/2d/3d/4d)
  - support both single and double precision

  General principles of the implementation can be found in:

  Accelerating the Non-equispaced Fast Fourier Transform on Commodity Graphics Hardware.
  T.S. Sørensen, T. Schaeffter, K.Ø. Noe, M.S. Hansen. 
  IEEE Transactions on Medical Imaging 2008; 27(4):538-547.

  Real-time Reconstruction of Sensitivity Encoded Radial Magnetic Resonance Imaging Using a Graphics Processing Unit.
  T.S. Sørensen, D. Atkinson, T. Schaeffter, M.S. Hansen.
  IEEE Transactions on Medical Imaging 2009; 28(12):1974-1985. 

  This example programme of the nnft utilizes golden ratio based radial trajectories 
  and outputs from an single precision input image ndarrays of the corresponding samples, trajectory, and density compensation weights.

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
  parms.add_parameter( 'd', COMMAND_LINE_STRING, 1, "Input image file name (.real)", true );
  parms.add_parameter( 'r', COMMAND_LINE_STRING, 1, "Result file name (.cplx)", true, "samples.cplx" );
  parms.add_parameter( 'o', COMMAND_LINE_INT,    1, "Oversampled matrix size", true );
  parms.add_parameter( 'p', COMMAND_LINE_INT,    1, "Number of profiles", true );
  parms.add_parameter( 's', COMMAND_LINE_INT,    1, "Samples per profiles", true );
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
  
  // Load image from disk
  timer = new GPUTimer("Loading image from disk");
  boost::shared_ptr< hoNDArray<_real> > host_image = read_nd_array<_real>((char*)parms.get_parameter('d')->get_string_value());
  delete timer;
   
  if( !(host_image->get_number_of_dimensions() == 2) ){
    GINFO_STREAM(endl << "Input image is not two-dimensional. Quitting.\n" << endl);
    return 1;
  }
  
  // Configuration from the command line
  uint64d2 matrix_size_os = uint64d2(parms.get_parameter('o')->get_int_value(), parms.get_parameter('o')->get_int_value());
  unsigned int num_profiles = parms.get_parameter('p')->get_int_value();
  unsigned int samples_per_profile = parms.get_parameter('s')->get_int_value();  
  _real kernel_width = parms.get_parameter('k')->get_float_value();

  uint64d2 matrix_size = from_std_vector<size_t,2>(*(host_image->get_dimensions().get()));
  _real alpha = (_real)matrix_size_os.vec[0]/(_real)matrix_size.vec[0];

  if( matrix_size.vec[0] != matrix_size.vec[1] ){
    GINFO_STREAM(
      endl << "For this samples application we only allow square input images. " <<
      endl << "The only reason being that only one oversampled matrix size is specified and the oversampling ratio must be consistent." << 
      endl);
  }
    
  // Upload host image to device, normalize, and convert to complex type
  timer = new GPUTimer("Uploading, normalizing and converting to complex");
  cuNDArray<_real> _image(*host_image);
  normalize( &_image, 1.0f );
  boost::shared_ptr< cuNDArray<_complext> > image = real_to_complex<_complext>( &_image );
  delete timer;
  
  // Setup resulting samples array
  vector<size_t> samples_dims; samples_dims.push_back( samples_per_profile ); samples_dims.push_back( num_profiles );
  cuNDArray<_complext> samples(samples_dims);
  
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
  plan.preprocess( *traj, NFFT_prep_mode::C2NC );
  delete timer;

  // Gridder
  timer = new GPUTimer("Computing nfft");
  plan.compute( *image, samples, 0, NFFT_comp_mode::BACKWARDS_NC2C );
  delete timer;

  //
  // Output result
  //
  
  timer = new GPUTimer("Output result to disk");
  boost::shared_ptr< hoNDArray<_complext> > host_samples = samples.to_host();
  write_nd_array<_complext>( host_samples.get(), (char*)parms.get_parameter('r')->get_string_value());
  delete timer;

  return 0;
}
