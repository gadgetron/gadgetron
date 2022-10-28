/*
  An example of how to register two 3d volumes using Cornelius-Kanade optical flow
*/

// Gadgetron includes
#include "cuCKOpticalFlowSolver.h"
#include "cuLinearResampleOperator.h"
#include "cuNDArray.h"
#include "hoNDArray_fileio.h"
#include "parameterparser.h"

// Std includes
#include <iostream>

using namespace Gadgetron;
using namespace std;

// Define desired precision
typedef float _real; 

int main(int argc, char** argv)
{
  //
  // Parse command line
  //

  ParameterParser parms;
  parms.add_parameter( 'f', COMMAND_LINE_STRING, 1, "Fixed image file name (.real)", true );
  parms.add_parameter( 'm', COMMAND_LINE_STRING, 1, "Moving image file name (.real)", true );
  parms.add_parameter( 'r', COMMAND_LINE_STRING, 1, "Result file name", true, "displacement_field.real" );
  parms.add_parameter( 'a', COMMAND_LINE_FLOAT,  1, "Regularization weight (alpha)", true, "0.05" );
  parms.add_parameter( 'b', COMMAND_LINE_FLOAT,  1, "Regularization weight (beta)", true, "1.0" );
  parms.add_parameter( 'l', COMMAND_LINE_INT,    1, "Number of multiresolution levels", true, "3" );
  
  parms.parse_parameter_list(argc, argv);
  if( parms.all_required_parameters_set() ){
    GINFO_STREAM(" Running registration with the following parameters: " << endl);
    parms.print_parameter_list();
  }
  else{
    GINFO_STREAM(" Some required parameters are missing: " << endl);
    parms.print_parameter_list();
    parms.print_usage();
    return 1;
  }
  
  // Load sample data from disk
  //
  
  boost::shared_ptr< hoNDArray<_real> > host_fixed = 
    read_nd_array<_real>((char*)parms.get_parameter('f')->get_string_value());

  boost::shared_ptr< hoNDArray<_real> > host_moving = 
    read_nd_array<_real>((char*)parms.get_parameter('m')->get_string_value());
  
  if( !host_fixed.get() || !host_moving.get() ){
    GINFO_STREAM(endl << "One of the input images is not found. Quitting!\n" << endl);
    return 1;
  }
  
  unsigned int num_fixed_dims = host_fixed->get_number_of_dimensions();
  unsigned int num_moving_dims = host_moving->get_number_of_dimensions();

  if( !(num_fixed_dims == 3 || num_fixed_dims == 4)  ){
    GINFO_STREAM(endl << "The fixed image is not three- or four-dimensional. Quitting!\n" << endl);
    return 1;
  }
  
  if( !(num_moving_dims == 3 || num_moving_dims == 4)  ){
    GINFO_STREAM(endl << "The moving image is not three- or four-dimensional. Quitting!\n" << endl);
    return 1;
  }

  // Upload host data to device
  //

  cuNDArray<_real> fixed_image(host_fixed.get());
  cuNDArray<_real> moving_image(host_moving.get());
  
  _real alpha = (_real) parms.get_parameter('a')->get_float_value();
  _real beta = (_real) parms.get_parameter('b')->get_float_value();

  unsigned int multires_levels = parms.get_parameter('l')->get_int_value();

  // Use trilinear interpolation for resampling
  //

  boost::shared_ptr< cuLinearResampleOperator<_real,3> > R( new cuLinearResampleOperator<_real,3>() );

  // Setup solver
  //
  
  cuCKOpticalFlowSolver<_real,3> CK;
  CK.set_interpolator( R );
  CK.set_output_mode( cuCKOpticalFlowSolver<_real,3>::OUTPUT_VERBOSE );  
  CK.set_max_num_iterations_per_level( 500 );
  CK.set_num_multires_levels( multires_levels );
  CK.set_alpha(alpha);
  CK.set_beta(beta);
  CK.set_limit(0.01f);
  
  // Run registration
  //

  boost::shared_ptr< cuNDArray<_real> > result = CK.solve( &fixed_image, &moving_image/*, true*/ );

  if( !result.get() ){
    GINFO_STREAM(endl << "Registration solver failed. Quitting!\n" << endl);
    return 1;
  }
  
  boost::shared_ptr< cuNDArray<_real> > deformed_moving = CK.deform( &moving_image, result );
  
  // All done, write out the result
  //

  boost::shared_ptr< hoNDArray<_real> > host_result = result->to_host();
  write_nd_array<_real>(host_result.get(), (char*)parms.get_parameter('r')->get_string_value());

  host_result = deformed_moving->to_host();
  write_nd_array<_real>(host_result.get(), "def_moving.real" );
  
  return 0;
}
