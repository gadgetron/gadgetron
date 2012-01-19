/*
	Total variation denoising based on the paper 
	"The Split Bregman Method for L1-Regularized Problems" by Tom Goldstein and Stanley Osher. 
	Siam J. Imaging Sciences. Vol. 2, No. 2, pp. 323–343.
*/

// Gadgetron includes
#include "cuNDArray.h"
#include "hoNDArray_fileio.h"
#include "ndarray_vector_td_utilities.h"
#include "cuSBSolver.h"
#include "cuCGSolver.h"
#include "cuIdentityOperator.h"
#include "cuPartialDerivativeOperator.h"
#include "parameterparser.h"

// Std includes
#include <iostream>

using namespace std;

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
    cout << " Running denoising with the following parameters: " << endl;
    parms.print_parameter_list();
  }
  else{
    cout << " Some required parameters are missing: " << endl;
    parms.print_parameter_list();
    parms.print_usage();
    return 1;
  }
    
  // Load sample data from disk
  boost::shared_ptr< hoNDArray<_real> > host_data = 
    read_nd_array<_real>((char*)parms.get_parameter('d')->get_string_value());

  if( !host_data.get() ){
    cout << endl << "Input image not found. Quitting!\n" << endl;
    return 1;
  }
  
  if( host_data->get_number_of_dimensions() != 2 ){
    cout << endl << "Input image is not two-dimensional. Quitting!\n" << endl;
    return 1;
  }
  
  // Upload host data to device
  cuNDArray<_real> data(host_data.get());
  
  _real mu = (_real) parms.get_parameter('m')->get_float_value();
  _real lambda = (_real)2.0*mu; // This is a good alround setting according to Goldstein et al.

  if( mu <= (_real) 0.0 ) {
    cout << endl << "Regularization parameter mu should be strictly positive. Quitting!\n" << endl;
    return 1;
  }

  unsigned int num_cg_iterations = parms.get_parameter('i')->get_int_value();
  unsigned int num_inner_iterations = parms.get_parameter('I')->get_int_value();
  unsigned int num_outer_iterations = parms.get_parameter('O')->get_int_value();
  
  // Setup regularization operators
  boost::shared_ptr< cuPartialDerivativeOperator<_real,_real,2> > Rx( new cuPartialDerivativeOperator<_real,_real,2>(0) ); 
  boost::shared_ptr< cuPartialDerivativeOperator<_real,_real,2> > Ry( new cuPartialDerivativeOperator<_real,_real,2>(1) ); 
  Rx->set_weight( lambda );
  Ry->set_weight( lambda );

  //
  // Setup conjugate gradient solver
  //
  
  // Define encoding matrix (identity)
  boost::shared_ptr< cuIdentityOperator<_real,_real> > E( new cuIdentityOperator<_real,_real>() );
  E->set_weight( mu );
  
  // Setup conjugate gradient solver
  boost::shared_ptr< cuCGSolver<_real,_real> > cg(new cuCGSolver<_real,_real>());
  cg->add_matrix_operator( E );   // encoding matrix
  cg->add_matrix_operator( Rx );  // regularization matrix
  cg->add_matrix_operator( Ry );  // regularization matrix
  cg->set_iterations( num_cg_iterations );
  cg->set_limit( 1e-4 );
  cg->set_output_mode( cuCGSolver<_real,_real>::OUTPUT_WARNINGS );  
  
  // Setup split-Bregman solver
  cuSBSolver<_real,_real> sb;
  sb.set_inner_solver( cg );
  sb.set_encoding_operator( E );
  //sb.add_regularization_operator( Rx ); // Anisotropic denoising
  //sb.add_regularization_operator( Ry ); // Anisotropic denoising
  sb.add_regularization_group_operator( Rx ); // Isotropic denoising
  sb.add_regularization_group_operator( Ry); // Isotropic denoising
  sb.add_group();
  sb.set_outer_iterations(num_outer_iterations);
  sb.set_inner_iterations(num_inner_iterations);
  sb.set_image_dimensions(data.get_dimensions());
  sb.set_output_mode( cuCGSolver<_real,_real>::OUTPUT_VERBOSE );
  
  // Run split-Bregman solver
  boost::shared_ptr< cuNDArray<_real> > sbresult = sb.solve(&data);

  // All done, write out the result
  boost::shared_ptr< hoNDArray<_real> > host_result = sbresult->to_host();
  write_nd_array<_real>(host_result.get(), (char*)parms.get_parameter('r')->get_string_value());
  
  return 0;
}
