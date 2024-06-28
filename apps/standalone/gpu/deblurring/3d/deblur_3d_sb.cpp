/*
  Deblurring using conjugate gradient solver.
*/

// Gadgetron includes
#include "cuNDArray.h"
#include "hoNDArray_fileio.h"
#include "cuSbcCgSolver.h"
#include "cuCgSolver.h"
#include "cuPartialDerivativeOperator.h"
#include "cuConvolutionOperator.h"
#include "parameterparser.h"

// Std includes
#include <iostream>

using namespace std;
using namespace Gadgetron;

// Define desired precision
typedef float _real; 
typedef complext<_real> _complext;

int main(int argc, char** argv)
{
  //
  // Parse command line
  //

  ParameterParser parms;
  parms.add_parameter( 'd', COMMAND_LINE_STRING, 1, "Blurred image file name (.cplx)", true, "blurred_image.cplx" );
  parms.add_parameter( 'k', COMMAND_LINE_STRING, 1, "Kernel image file name (.cplx)", true, "kernel_image.cplx" );
  parms.add_parameter( 'r', COMMAND_LINE_STRING, 1, "Result file name", true, "sb_deblurred_image.cplx" );
  parms.add_parameter( 'i', COMMAND_LINE_INT,    1, "Number of cg iterations", true, "20" );
  parms.add_parameter( 'I', COMMAND_LINE_INT,    1, "Number of sb inner iterations", true, "1" );
  parms.add_parameter( 'O', COMMAND_LINE_INT,    1, "Number of sb outer iterations", true, "50" );
  parms.add_parameter( 'M', COMMAND_LINE_FLOAT,  1, "Mu", true, "1.0" );
  parms.add_parameter( 'L', COMMAND_LINE_FLOAT,  1, "Lambda", true, "1.0" );

  parms.parse_parameter_list(argc, argv);
  if( parms.all_required_parameters_set() ){
    GINFO_STREAM(" Running deblurring with the following parameters: " << endl);
    parms.print_parameter_list();
  }
  else{
    GINFO_STREAM(" Some required parameters are missing: " << endl);
    parms.print_parameter_list();
    parms.print_usage();
    return 1;
  }
    
  // Load sample data from disk
  boost::shared_ptr< hoNDArray<_complext> > host_data = 
    read_nd_array<_complext>((char*)parms.get_parameter('d')->get_string_value());

  boost::shared_ptr< hoNDArray<_complext> > host_kernel = 
    read_nd_array<_complext>((char*)parms.get_parameter('k')->get_string_value());
   
  if( !(host_data->get_number_of_dimensions() == 3) || !(host_kernel->get_number_of_dimensions() == 3) ){
    GINFO_STREAM(endl << "Input data (image/kernel) is not two-dimensional. Quitting!\n" << endl);
    return 1;
  }

  // Upload host data to device
  cuNDArray<_complext> data(*host_data);
  cuNDArray<_complext> kernel(*host_kernel);
  
  unsigned int num_cg_iterations = parms.get_parameter('i')->get_int_value();
  unsigned int num_inner_iterations = parms.get_parameter('I')->get_int_value();
  unsigned int num_outer_iterations = parms.get_parameter('O')->get_int_value();
  
  // Setup regularization operators
  boost::shared_ptr< cuPartialDerivativeOperator<_complext,3> > Rx( new cuPartialDerivativeOperator<_complext,3>(0) );
  boost::shared_ptr< cuPartialDerivativeOperator<_complext,3> > Ry( new cuPartialDerivativeOperator<_complext,3>(1) );
  boost::shared_ptr< cuPartialDerivativeOperator<_complext,3> > Rz( new cuPartialDerivativeOperator<_complext,3>(2) );
  
  _real mu = (_real) parms.get_parameter('M')->get_float_value();
  _real lambda = (_real) parms.get_parameter('L')->get_float_value();

  if( mu <= (_real) 0.0 ) {
    GINFO_STREAM(endl << "Regularization parameter mu should be strictly positive. Quitting!\n" << endl);
    return 1;
  }

  Rx->set_weight( lambda );
  Rx->set_domain_dimensions(data.get_dimensions());
  Rx->set_codomain_dimensions(data.get_dimensions());

  Ry->set_weight( lambda );
  Ry->set_domain_dimensions(data.get_dimensions());
  Ry->set_codomain_dimensions(data.get_dimensions());

  Rz->set_weight( lambda );
  Rz->set_domain_dimensions(data.get_dimensions());
  Rz->set_codomain_dimensions(data.get_dimensions());

  //
  // Setup conjugate gradients solver
  //

  // Define encoding matrix
  boost::shared_ptr< cuConvolutionOperator<_real,3> > E( new cuConvolutionOperator<_real,3>() );  
  E->set_kernel( &kernel );
  E->set_weight( mu );
  E->set_domain_dimensions(data.get_dimensions());
  E->set_codomain_dimensions(data.get_dimensions());
  
  // Setup split-Bregman solver
  cuSbcCgSolver<_complext> sb;
  sb.set_encoding_operator( E );
  sb.add_regularization_group_operator( Rx ); 
  sb.add_regularization_group_operator( Ry ); 
  sb.add_group();
  sb.add_regularization_operator( Rz ); 
  sb.set_max_outer_iterations(num_outer_iterations);
  sb.set_max_inner_iterations(num_inner_iterations);
  sb.set_output_mode( cuSbcCgSolver< _complext>::OUTPUT_VERBOSE );

  sb.get_inner_solver()->set_max_iterations( num_cg_iterations );
  sb.get_inner_solver()->set_tc_tolerance( 1e-8 );
  sb.get_inner_solver()->set_output_mode( cuCgSolver<_complext>::OUTPUT_WARNINGS );

  // Run split-Bregman solver
  boost::shared_ptr< cuNDArray<_complext> > sbresult = sb.solve(&data);

  // All done, write out the result
  
  boost::shared_ptr< hoNDArray<_complext> > host_result = sbresult->to_host();
  write_nd_array<_complext>(host_result.get(), (char*)parms.get_parameter('r')->get_string_value());
    
  boost::shared_ptr< hoNDArray<_real> > host_norm = abs(sbresult.get())->to_host();
  write_nd_array<_real>( host_norm.get(), "sb_deblurred_image.real" );  

  return 0;
}
