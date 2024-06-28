/*
  An example of how to register two 2d images using Horn-Schunk optical flow
*/

// Gadgetron includes
#include "cuHSOpticalFlowSolver.h"
#include "cuLinearResampleOperator.h"
#include "cuNDArray.h"
#include "hoNDArray_fileio.h"
#include "hoNDArray_utils.h"
#include "parameterparser.h"
#include "cuCGHSOFSolver.h"
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
  parms.add_parameter( 'a', COMMAND_LINE_FLOAT,  1, "Regularization weight (alpha)", true, "0.1" );
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

  if( !(num_fixed_dims == 2 || num_fixed_dims == 3)  ){
    GINFO_STREAM(endl << "The fixed image is not two- or three-dimensional. Quitting!\n" << endl);
    return 1;
  }
  
  if( !(num_moving_dims == 2 || num_moving_dims == 3)  ){
    GINFO_STREAM(endl << "The moving image is not two- or three-dimensional. Quitting!\n" << endl);
    return 1;
  }
  
  if( num_fixed_dims < num_moving_dims  ){
    *host_fixed = expand( *host_fixed, host_moving->get_size(2) );
    num_fixed_dims = host_fixed->get_number_of_dimensions();
  }

  if( num_moving_dims < num_moving_dims  ){
    *host_moving = expand( *host_moving, host_fixed->get_size(2) );
    num_moving_dims = host_moving->get_number_of_dimensions();
  }

  // Upload host data to device
  //

  cuNDArray<_real> fixed_image(*host_fixed);
  cuNDArray<_real> moving_image(*host_moving);
  
  _real alpha = (_real) parms.get_parameter('a')->get_float_value();

  unsigned int multires_levels = parms.get_parameter('l')->get_int_value();

  // Use bilinear interpolation for resampling
  //

  boost::shared_ptr< cuLinearResampleOperator<_real,2> > R( new cuLinearResampleOperator<_real,2>() );

  // Setup solver
  //
  
  cuCGHSOFSolver<_real,2> HS;
  HS.set_interpolator( R );
  HS.set_output_mode( cuCGHSOFSolver<_real,2>::OUTPUT_VERBOSE );
  HS.get_solver()->set_max_iterations( 100 );
  HS.get_solver()->set_output_mode(cuCgSolver<_real>::OUTPUT_VERBOSE);
  HS.set_num_multires_levels( multires_levels );
  HS.set_alpha(alpha);
  

  // Run registration
  //

  boost::shared_ptr< cuNDArray<_real> > result = HS.solve( &fixed_image, &moving_image );

  if( !result.get() ){
    GINFO_STREAM(endl << "Registration solver failed. Quitting!\n" << endl);
    return 1;
  }
  
  boost::shared_ptr< cuNDArray<_real> > deformed_moving = HS.deform( &moving_image, result );
  
  // All done, write out the result
  //

  boost::shared_ptr< hoNDArray<_real> > host_result = result->to_host();
  write_nd_array<_real>(host_result.get(), (char*)parms.get_parameter('r')->get_string_value());

  host_result = deformed_moving->to_host();
  write_nd_array<_real>(host_result.get(), "def_moving.real" );
  
  return 0;
}
