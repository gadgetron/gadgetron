/*
  An example of how to register two 2d images using Horn-Schunk optical flow
*/

// Gadgetron includes
#include "hoHSOpticalFlowSolver.h"
#include "hoLinearResampleOperator.h"
#include "hoNDArray.h"
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
  parms.add_parameter( 'a', COMMAND_LINE_FLOAT,  1, "Regularization weight (alpha)", true, "0.1" );
  
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
  
  boost::shared_ptr< hoNDArray<_real> > fixed_image = 
    read_nd_array<_real>((char*)parms.get_parameter('f')->get_string_value());
  
  boost::shared_ptr< hoNDArray<_real> > moving_image = 
    read_nd_array<_real>((char*)parms.get_parameter('m')->get_string_value());
  
  if( !moving_image.get() || !fixed_image.get() ){
    GINFO_STREAM(endl << "One of the input images is not found. Quitting!\n" << endl);
    return 1;
  }
  
  size_t num_fixed_dims = fixed_image->get_number_of_dimensions();
  size_t num_moving_dims = moving_image->get_number_of_dimensions();

  if( !(num_fixed_dims == 2 || num_fixed_dims == 3)  ){
    GINFO_STREAM(endl << "The fixed image is not two- or three-dimensional. Quitting!\n" << endl);
    return 1;
  }
  
  if( !(num_moving_dims == 2 || num_moving_dims == 3)  ){
    GINFO_STREAM(endl << "The moving image is not two- or three-dimensional. Quitting!\n" << endl);
    return 1;
  }
    
  _real alpha = (_real) parms.get_parameter('a')->get_float_value();

  // Use bilinear interpolation for resampling
  //

  boost::shared_ptr< hoLinearResampleOperator<_real,2> > R( new hoLinearResampleOperator<_real,2>() );

  // Setup solver
  //
  
  hoHSOpticalFlowSolver<_real,2> HS;
  HS.set_interpolator( R );
  HS.set_output_mode( hoHSOpticalFlowSolver<_real,2>::OUTPUT_VERBOSE );  
  HS.set_num_multires_levels( 4 );
  HS.set_max_num_iterations_per_level( 500 );
  HS.set_alpha(alpha);
  HS.set_limit(0.01f);
  
  // Run registration
  //

  boost::shared_ptr< hoNDArray<_real> > result = HS.solve( fixed_image.get(), moving_image.get() );

  if( !result.get() ){
    GINFO_STREAM(endl << "Registration solver failed. Quitting!\n" << endl);
    return 1;
  }
  
  boost::shared_ptr< hoNDArray<_real> > deformed_moving = HS.deform( moving_image.get(), result );
  
  // All done, write out the result
  //

  write_nd_array<_real>(result.get(), (char*)parms.get_parameter('r')->get_string_value());
  write_nd_array<_real>(deformed_moving.get(), "def_moving.real" );
  
  return 0;
}
