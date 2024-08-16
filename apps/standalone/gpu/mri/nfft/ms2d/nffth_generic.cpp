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
  and outputs gridded images from 2D multislice input ndarrays of the corresponding samples, trajectory, and density compensation weights.

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
  parms.add_parameter( 't', COMMAND_LINE_STRING, 1, "Input trajectories file name (.real)", true );
  parms.add_parameter( 'w', COMMAND_LINE_STRING, 1, "Input density compensation weights file name (.real)", true );
  parms.add_parameter( 'r', COMMAND_LINE_STRING, 1, "Output image file name (.cplx)", true, "result.cplx" );
  parms.add_parameter( 'm', COMMAND_LINE_INT,    1, "Matrix size", true );
  parms.add_parameter( 'o', COMMAND_LINE_INT,    1, "Oversampled matrix size", true );
  parms.add_parameter( 'f', COMMAND_LINE_INT,    1, "#frames/reconstruction (a negative value means all)", true, "-1" );
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
  
  // Load data from disk
  timer = new GPUTimer("Loading data from disk");
  boost::shared_ptr< hoNDArray<_complext> > host_samples = read_nd_array<_complext>((char*)parms.get_parameter('d')->get_string_value());
  boost::shared_ptr< hoNDArray<_reald2> >   host_traj    = read_nd_array<_reald2>  ((char*)parms.get_parameter('t')->get_string_value());
  boost::shared_ptr< hoNDArray<_real> >     host_dcw     = read_nd_array<_real>    ((char*)parms.get_parameter('w')->get_string_value());
  delete timer;

  if( !(host_samples->get_number_of_dimensions() == 2 && host_traj->get_number_of_dimensions() == 2) ){
    GINFO_STREAM(endl << "Samples/trajectory arrays must be two-dimensional: (dim 0: samples/profile x #profiles/frame; dim 1: #frames). Quitting.\n" << endl);
    return 1;
  }

  // Configuration from the command line
  uint64d2 matrix_size = uint64d2(parms.get_parameter('m')->get_int_value(), parms.get_parameter('m')->get_int_value());
  uint64d2 matrix_size_os = uint64d2(parms.get_parameter('o')->get_int_value(), parms.get_parameter('o')->get_int_value());
  int frames_per_reconstruction = parms.get_parameter('f')->get_int_value();  
  _real kernel_width = parms.get_parameter('k')->get_float_value();
  _real alpha = (_real)matrix_size_os.vec[0]/(_real)matrix_size.vec[0];
  
  unsigned int num_frames = host_traj->get_size(1);  

  if( frames_per_reconstruction < 0 ) frames_per_reconstruction = num_frames;
  if( (unsigned int)frames_per_reconstruction > num_frames ) frames_per_reconstruction = num_frames;
  
  // Setup resulting image array
  vector<size_t> image_dims = to_std_vector(matrix_size); 
  image_dims.push_back((num_frames/frames_per_reconstruction)*frames_per_reconstruction);
  cuNDArray<_complext> image(image_dims);
  clear(&image);
  
  // Initialize plan
  timer = new GPUTimer("Initializing plan");
  plan_type plan( matrix_size, matrix_size_os, kernel_width );
  delete timer;

  // Upload arrays to device
  cuNDArray<_complext> _samples(*host_samples);
  cuNDArray<_reald2> _trajectory(*host_traj);
  cuNDArray<_real> dcw(*host_dcw);

  std::vector<size_t> dims_recon;
  dims_recon.push_back(host_samples->get_size(0));
  dims_recon.push_back(frames_per_reconstruction);

  for( unsigned int iteration = 0; iteration < num_frames/frames_per_reconstruction; iteration++ ) {
    
    // Set samples/trajectory for sub-frames
    cuNDArray<_complext> samples( dims_recon, _samples.get_data_ptr()+iteration*dims_recon[0]*dims_recon[1] );
    cuNDArray<_reald2> trajectory( dims_recon, _trajectory.get_data_ptr()+iteration*dims_recon[0]*dims_recon[1] );

    // Preprocess
    timer = new GPUTimer("NFFT preprocessing");
    plan.preprocess( trajectory, NFFT_prep_mode::NC2C );
    delete timer;
    
    std::vector<size_t> image_dims = to_std_vector(matrix_size); 
    image_dims.push_back(frames_per_reconstruction);
    cuNDArray<_complext> tmp_image(image_dims, image.get_data_ptr()+iteration*prod(matrix_size)*frames_per_reconstruction);

    // Gridder
    timer = new GPUTimer("Computing adjoint nfft (gridding)");
    plan.compute( samples, tmp_image, &dcw, NFFT_comp_mode::BACKWARDS_NC2C );
    delete timer;
  }
  
  //
  // Output result
  //
  
  timer = new GPUTimer("Output result to disk");
  boost::shared_ptr< hoNDArray<_complext> > host_image = image.to_host();
  write_nd_array<_complext>( host_image.get(), (char*)parms.get_parameter('r')->get_string_value() );
  write_nd_array<_real>( abs(&image)->to_host().get(), "result.real" );
  delete timer;

  return 0;
}
