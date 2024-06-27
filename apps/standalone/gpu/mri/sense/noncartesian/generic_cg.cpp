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
#include "cuNDArray_operators.h"
#include "cuNDArray_elemwise.h"
#include "cuNDArray_blas.h"
#include "cuNDArray_utils.h"
#include "cuNDArray_reductions.h"
#include "cuNonCartesianSenseOperator.h"
#include "cuCgPreconditioner.h"
#include "cuImageOperator.h"
#include "cuCgSolver.h"
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
typedef cuNFFT_plan<_real,2> plan_type;

int main( int argc, char** argv) 
{

  //
  // Parse command line
  //

  ParameterParser parms;
  parms.add_parameter( 'd', COMMAND_LINE_STRING, 1, "Input samples file name (.cplx)", true );
  parms.add_parameter( 't', COMMAND_LINE_STRING, 1, "Input trajectories file name (.real)", true );
  parms.add_parameter( 'w', COMMAND_LINE_STRING, 1, "Input density compensation weights file name (.real)", true );
  parms.add_parameter( 'c', COMMAND_LINE_STRING, 1, "Input coil sensitivity maps file name (.cplx)", true );
  parms.add_parameter( 'g', COMMAND_LINE_STRING, 1, "Input regularization image file name (.cplx)", true );
  parms.add_parameter( 'r', COMMAND_LINE_STRING, 1, "Output image file name (.cplx)", true, "result.cplx" );
  parms.add_parameter( 'i', COMMAND_LINE_INT,    1, "Number of iterations", true, "10" );
  parms.add_parameter( 'l', COMMAND_LINE_FLOAT,  1, "Regularization weight", true, "0.3" );
  parms.add_parameter( 'k', COMMAND_LINE_FLOAT,  1, "Kernel width", true, "5.5" );
  parms.add_parameter( 'a', COMMAND_LINE_FLOAT,  1, "Oversampling factor", true, "2.0" );

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
  boost::shared_ptr< hoNDArray<_complext> > host_samples = read_nd_array<_complext> ((char*)parms.get_parameter('d')->get_string_value());
  boost::shared_ptr< hoNDArray<_reald2> >   host_traj    = read_nd_array<_reald2>   ((char*)parms.get_parameter('t')->get_string_value());
  boost::shared_ptr< hoNDArray<_real> >     host_dcw     = read_nd_array<_real>     ((char*)parms.get_parameter('w')->get_string_value());
  boost::shared_ptr< hoNDArray<_complext> > host_csm     = read_nd_array<_complext> ((char*)parms.get_parameter('c')->get_string_value());
  boost::shared_ptr< hoNDArray<_complext> > host_reg     = read_nd_array<_complext> ((char*)parms.get_parameter('g')->get_string_value());
  delete timer;
   
  /* {
    std::vector<size_t> dims;
    dims.push_back(host_traj->get_size(0));
    dims.push_back(host_samples->get_number_of_elements()/dims[0]);
    host_samples->reshape(&dims);
    } */

  if( !(host_samples->get_number_of_dimensions() == 2 && host_traj->get_number_of_dimensions() == 2) ){
    GINFO_STREAM(endl << "Samples/trajectory arrays must be two-dimensional: (dim 0: samples/profile x #profiles/frame; dim 1: #frames). Quitting.\n" << endl);
    return 1;
  }

  if( !(host_csm->get_number_of_dimensions() == 3 )){
    GINFO_STREAM(endl << "Coil sensitivity maps must be three-dimensional. Quitting.\n" << endl);
    return 1;
  }

  if( !(host_reg->get_number_of_dimensions() == 2 )){
    GINFO_STREAM(endl << "Regularization image must be two-dimensional. Quitting.\n" << endl);
    return 1;
  }

  // Configuration from the command line
  uint64d2 matrix_size = uint64d2(host_csm->get_size(0), host_csm->get_size(0));
  size_t _matrix_size_os = size_t((float)matrix_size[0]*parms.get_parameter('a')->get_float_value());
  uint64d2 matrix_size_os = uint64d2(_matrix_size_os, _matrix_size_os);
  int num_iterations = parms.get_parameter('i')->get_int_value();
  _real kernel_width = parms.get_parameter('k')->get_float_value();
  _real alpha = parms.get_parameter('a')->get_float_value();
  _real kappa = parms.get_parameter('l')->get_float_value();
  
  unsigned int num_frames = host_traj->get_size(1);  
  unsigned int num_coils = host_csm->get_size(2);

  std::vector<size_t> recon_dims = to_std_vector(matrix_size);
  recon_dims.push_back(num_frames);

  // Upload arrays to device
  cuNDArray<_complext> samples(*host_samples);
  cuNDArray<_reald2> trajectory(*host_traj);
  boost::shared_ptr< cuNDArray<_complext> > csm( new cuNDArray<_complext>(*host_csm));
  boost::shared_ptr< cuNDArray<_complext> > reg_image( new cuNDArray<_complext>(*host_reg));
  boost::shared_ptr< cuNDArray<_real> > dcw( new cuNDArray<_real>(*host_dcw));

  // Define encoding matrix for non-Cartesian SENSE
  boost::shared_ptr< cuNonCartesianSenseOperator<_real,2> > E( new cuNonCartesianSenseOperator<_real,2>() );  
  E->setup( matrix_size, matrix_size_os, kernel_width );
  E->set_dcw(dcw) ;
  E->set_csm(csm);
  E->set_domain_dimensions(&recon_dims);
  E->set_codomain_dimensions(samples.get_dimensions().get());
  E->preprocess(&trajectory);
  
  // Define regularization operator
  boost::shared_ptr< cuImageOperator<_complext> > R( new cuImageOperator<_complext>() );
  R->set_weight( kappa );
  R->compute( reg_image.get() );

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

  // Setup conjugate gradient solver
  cuCgSolver<_complext> cg;
  cg.set_preconditioner ( D );           // preconditioning matrix
  cg.set_max_iterations( num_iterations );
  cg.set_tc_tolerance( 1e-6 );
  cg.set_output_mode( cuCgSolver< _complext>::OUTPUT_VERBOSE );
  cg.set_encoding_operator( E );        // encoding matrix
  cg.add_regularization_operator( R );  // regularization matrix

  //
  // Invoke conjugate gradient solver
  //
  
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
  write_nd_array<_complext>( host_image.get(), (char*)parms.get_parameter('r')->get_string_value() );
  write_nd_array<_real>( abs(cgresult.get())->to_host().get(), "result.real" );
  delete timer;

  return 0;
}
