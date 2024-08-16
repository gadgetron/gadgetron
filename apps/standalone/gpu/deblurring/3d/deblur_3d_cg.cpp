/*
  Deblurring using conjugate gradient solver.
*/

// Gadgetron includes
#include "cuNDArray.h"
#include "hoNDArray_fileio.h"
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
  parms.add_parameter( 'r', COMMAND_LINE_STRING, 1, "Result file name", true, "cg_deblurred_image.cplx" );
  parms.add_parameter( 'i', COMMAND_LINE_INT,    1, "Number of iterations", true, "25" );
  parms.add_parameter( 'K', COMMAND_LINE_FLOAT,  1, "Regularization weight", true, "0.1" );

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
  
  _real kappa = (_real) parms.get_parameter('K')->get_float_value();
  unsigned int num_iterations = parms.get_parameter('i')->get_int_value();
  
  // Setup regularization operators
  boost::shared_ptr< cuPartialDerivativeOperator<_complext,3> > Rx( new cuPartialDerivativeOperator<_complext,3>(0) );
  boost::shared_ptr< cuPartialDerivativeOperator<_complext,3> > Ry( new cuPartialDerivativeOperator<_complext,3>(1) );
  boost::shared_ptr< cuPartialDerivativeOperator<_complext,3> > Rz( new cuPartialDerivativeOperator<_complext,3>(2) );

  Rx->set_weight( kappa );
  Ry->set_weight( kappa );
  Rz->set_weight( kappa );
     
  //
  // Setup conjugate gradients solver
  //

  // Define encoding matrix
  boost::shared_ptr< cuConvolutionOperator<_real,3> > E( new cuConvolutionOperator<_real,3>() );
  E->set_kernel( &kernel );
  E->set_domain_dimensions(data.get_dimensions());
    
  // Setup conjugate gradient solver
  cuCgSolver<_complext> cg;
  cg.set_encoding_operator( E );                         // encoding matrix
  if( kappa>0.0 ) cg.add_regularization_operator( Rx );  // regularization matrix
  if( kappa>0.0 ) cg.add_regularization_operator( Ry );  // regularization matrix
  if( kappa>0.0 ) cg.add_regularization_operator( Rz );  // regularization matrix
  cg.set_max_iterations( num_iterations );
  cg.set_tc_tolerance( 1e-12 );
  cg.set_output_mode( cuCgSolver<_complext>::OUTPUT_VERBOSE );
                
  // Form right hand side
  cuNDArray<_complext> rhs; rhs.create(data.get_dimensions());
  E->mult_MH( &data, &rhs );
  
  //
  // Conjugate gradient solver
  //
  
  boost::shared_ptr< cuNDArray<_complext> > cgresult = cg.solve_from_rhs(&rhs);

  // All done, write out the result
  
  boost::shared_ptr< hoNDArray<_complext> > host_result = cgresult->to_host();
  write_nd_array<_complext>(host_result.get(), (char*)parms.get_parameter('r')->get_string_value());
    
  boost::shared_ptr< hoNDArray<_real> > host_norm = abs(cgresult.get())->to_host();
  write_nd_array<_real>( host_norm.get(), "cg_deblurred_image.real" );  

  return 0;
}
