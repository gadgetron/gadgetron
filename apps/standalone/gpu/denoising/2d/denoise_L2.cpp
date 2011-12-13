/*
  Deblurring using conjugate gradient solver.
*/

// Gadgetron includes
#include "cuNDArray.h"
#include "hoNDArray_fileio.h"
#include "ndarray_vector_td_utilities.h"
#include "cuCGSolver.h"
#include "cuPartialDerivativeOperator.h"
#include "cuConvolutionOperator.h"
#include "parameterparser.h"

// Std includes
#include <iostream>

using namespace std;

// Define desired precision (note that decent deblurring of noisy images requires double precision)
typedef double _real; 
typedef complext<_real>::Type _complext;

int main(int argc, char** argv)
{
  //
  // Parse command line
  //

  ParameterParser parms;
  parms.add_parameter( 'd', COMMAND_LINE_STRING, 1, "Blurred image file name (.cplx)", true, "blurred_image.cplx" );
  parms.add_parameter( 'k', COMMAND_LINE_STRING, 1, "Kernel image file name (.cplx)", true, "kernel_image.cplx" );
  parms.add_parameter( 'r', COMMAND_LINE_STRING, 1, "Result file name", true, "cg_deblurred_image.cplx" );
  parms.add_parameter( 'i', COMMAND_LINE_INT,    1, "Number of iterations", true, "10" );
  parms.add_parameter( 'K', COMMAND_LINE_FLOAT,  1, "Regularization weight", true, "0.1" );

  parms.parse_parameter_list(argc, argv);
  if( parms.all_required_parameters_set() ){
    cout << " Running deblurring with the following parameters: " << endl;
    parms.print_parameter_list();
  }
  else{
    cout << " Some required parameters are missing: " << endl;
    parms.print_parameter_list();
    parms.print_usage();
    return 1;
  }
    
  // Load sample data from disk
  boost::shared_ptr< hoNDArray<_complext> > host_data = 
    read_nd_array<_complext>((char*)parms.get_parameter('d')->get_string_value());

  boost::shared_ptr< hoNDArray<_complext> > host_kernel = 
    read_nd_array<_complext>((char*)parms.get_parameter('k')->get_string_value());
   
  if( !(host_data->get_number_of_dimensions() == 2) || !(host_kernel->get_number_of_dimensions() == 2) ){
    cout << endl << "Input data (image/kernel) is not two-dimensional. Quitting!\n" << endl;
    return 1;
  }

  // Upload host data to device
  cuNDArray<_complext> data(host_data.get());
  cuNDArray<_complext> kernel(host_kernel.get());
  
  _real kappa = (_real) parms.get_parameter('K')->get_float_value();
  unsigned int num_iterations = parms.get_parameter('i')->get_int_value();
  
  // Setup regularization operators
  boost::shared_ptr< cuPartialDerivativeOperator<_real,_complext,2> > Rx( new cuPartialDerivativeOperator<_real,_complext,2>(0) ); 
  boost::shared_ptr< cuPartialDerivativeOperator<_real,_complext,2> > Ry( new cuPartialDerivativeOperator<_real,_complext,2>(1) ); 
  Rx->set_weight( kappa );
  Ry->set_weight( kappa );
     
  //
  // Setup conjugate gradients solver
  //

  // Define encoding matrix
  boost::shared_ptr< cuConvolutionOperator<_real> > 
    E( new cuConvolutionOperator<_real>() );
  E->set_kernel( &kernel );
    
  // Setup conjugate gradient solver
  cuCGSolver<_real, _complext> cg;
  cg.add_matrix_operator( E );                  // encoding matrix
  if( kappa>0.0 ) cg.add_matrix_operator( Rx );  // regularization matrix
  if( kappa>0.0 ) cg.add_matrix_operator( Ry );  // regularization matrix
  cg.set_iterations( num_iterations );
  cg.set_limit( 1e-12 );
  cg.set_output_mode( cuCGSolver<_real, _complext>::OUTPUT_VERBOSE );
                
  // Form right hand side
  cuNDArray<_complext> rhs; rhs.create(data.get_dimensions().get());
  E->mult_MH( &data, &rhs );
  
  //
  // Conjugate gradient solver
  //
  
  boost::shared_ptr< cuNDArray<_complext> > cgresult = cg.solve(&rhs);

  // All done, write out the result
  
  boost::shared_ptr< hoNDArray<_complext> > host_result = cgresult->to_host();
  write_nd_array<_complext>(host_result.get(), (char*)parms.get_parameter('r')->get_string_value());
    
  boost::shared_ptr< hoNDArray<_real> > host_norm = cuNDA_norm<_real,2>(cgresult.get())->to_host();
  write_nd_array<_real>( host_norm.get(), "cg_deblurred_image.real" );  

  boost::shared_ptr< hoNDArray<_real> > host_rhs = cuNDA_norm<_real,2>(&rhs)->to_host();
  write_nd_array<_real>( host_rhs.get(), "rhs.real" );  

  return 0;
}

