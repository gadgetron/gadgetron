/*
  Total variation denoising based on the paper 
  "The Split Bregman Method for L1-Regularized Problems" by Tom Goldstein and Stanley Osher. 
  Siam J. Imaging Sciences. Vol. 2, No. 2, pp. 323-343.
*/

// Gadgetron includes
#include "hoNDArray.h"
#include "hoNDArray_fileio.h"
#include "hoSbCgSolver.h"
#include "hoIdentityOperator.h"
#include "hoPartialDerivativeOperator.h"
#include "parameterparser.h"

// Std includes
#include <iostream>

using namespace std;
using namespace Gadgetron;

// Define desired precision
typedef float _real; 

int main(int argc, char** argv)
{
  //
  // Parse command line
  //

  ParameterParser parms;
  parms.add_parameter( 'd', COMMAND_LINE_STRING, 1, "Noisy image file name (.real)", true );
  parms.add_parameter( 'r', COMMAND_LINE_STRING, 1, "Result file name", true, "denoised_image_TV.real" );
  parms.add_parameter( 'i', COMMAND_LINE_INT,    1, "Number of cg iterations", true, "20" );
  parms.add_parameter( 'I', COMMAND_LINE_INT,    1, "Number of sb inner iterations", true, "1" );
  parms.add_parameter( 'O', COMMAND_LINE_INT,    1, "Number of sb outer iterations", true, "10" );
  parms.add_parameter( 'm', COMMAND_LINE_FLOAT,  1, "Regularization weight (mu)", true, "25.0" );

  parms.parse_parameter_list(argc, argv);
  if( parms.all_required_parameters_set() ){
    GINFO_STREAM(" Running denoising with the following parameters: " << endl);
    parms.print_parameter_list();
  }
  else{
    GINFO_STREAM(" Some required parameters are missing: " << endl);
    parms.print_parameter_list();
    parms.print_usage();
    return 1;
  }
    
  // Load sample data from disk
  boost::shared_ptr< hoNDArray<_real> > data = 
    read_nd_array<_real>((char*)parms.get_parameter('d')->get_string_value());

  if( !data.get() ){
    GINFO_STREAM(endl << "Input image not found. Quitting!\n" << endl);
    return 1;
  }
  
  if( data->get_number_of_dimensions() != 2 ){
    GINFO_STREAM(endl << "Input image is not two-dimensional. Quitting!\n" << endl);
    return 1;
  }
  
  _real mu = (_real) parms.get_parameter('m')->get_float_value();
  _real lambda = (_real)2.0*mu; // This is a good alround setting according to Goldstein et al.

  if( mu <= (_real) 0.0 ) {
    GINFO_STREAM(endl << "Regularization parameter mu should be strictly positive. Quitting!\n" << endl);
    return 1;
  }
  
  size_t num_cg_iterations = parms.get_parameter('i')->get_int_value();
  size_t num_inner_iterations = parms.get_parameter('I')->get_int_value();
  size_t num_outer_iterations = parms.get_parameter('O')->get_int_value();
  
  // Setup regularization operators
  boost::shared_ptr< hoPartialDerivativeOperator<_real,2> > Rx( new hoPartialDerivativeOperator<_real,2>(0) );
  Rx->set_weight( lambda );
  Rx->set_domain_dimensions(data->get_dimensions());
  Rx->set_codomain_dimensions(data->get_dimensions());
  
  boost::shared_ptr< hoPartialDerivativeOperator<_real,2> > Ry( new hoPartialDerivativeOperator<_real,2>(1) );
  Ry->set_weight( lambda );
  Ry->set_domain_dimensions(data->get_dimensions());
  Ry->set_codomain_dimensions(data->get_dimensions());
  
  // Define encoding operator (identity)
  boost::shared_ptr< identityOperator<hoNDArray<_real> > > E( new identityOperator<hoNDArray<_real> >() );
  E->set_weight( mu );
  E->set_domain_dimensions(data->get_dimensions());
  E->set_codomain_dimensions(data->get_dimensions());
  
  // Setup split-Bregman solver
  hoSbCgSolver<_real> sb;
  sb.set_encoding_operator( E );
  //sb.add_regularization_operator( Rx ); // Anisotropic denoising
  //sb.add_regularization_operator( Ry ); // Anisotropic denoising
  sb.add_regularization_group_operator( Rx ); // Isotropic denoising
  sb.add_regularization_group_operator( Ry); // Isotropic denoising
  sb.add_group();
  sb.set_max_outer_iterations(num_outer_iterations);
  sb.set_max_inner_iterations(num_inner_iterations);
  sb.set_output_mode( hoCgSolver<_real>::OUTPUT_VERBOSE );
  
  // Setup inner conjugate gradient solver
  sb.get_inner_solver()->set_max_iterations( num_cg_iterations );
  sb.get_inner_solver()->set_tc_tolerance( 1e-4 );
  sb.get_inner_solver()->set_output_mode( hoCgSolver<_real>::OUTPUT_WARNINGS );
  
  // Run split-Bregman solver
  boost::shared_ptr< hoNDArray<_real> > sbresult = sb.solve(data.get());
  
  // All done, write out the result
  write_nd_array<_real>(sbresult.get(), (char*)parms.get_parameter('r')->get_string_value());
  
  return 0;
}
