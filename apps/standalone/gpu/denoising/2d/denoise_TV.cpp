/*
  Total variation denoising based on the paper 
  "The Split Bregman Method for L1-Regularized Problems" by Tom Goldstein and Stanley Osher. 
  Siam J. Imaging Sciences. Vol. 2, No. 2, pp. 323-343.
*/

// Gadgetron includes
#include "cuNDArray.h"
#include "hoNDArray_fileio.h"
#include "cuSbCgSolver.h"
#include "cuCgSolver.h"
#include "identityOperator.h"
#include "cuPartialDerivativeOperator.h"
#include "parameterparser.h"
#include "cuNDDWT.h"
#include "cuDWTOperator.h"
#include <boost/make_shared.hpp>
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
  parms.add_parameter( 'l', COMMAND_LINE_FLOAT,  1, "Total variation weight (lambda)", true, "50.0" );
  parms.add_parameter( 'm', COMMAND_LINE_FLOAT,  1, "Regularization weight (mu)", true, "25.0" );
  parms.add_parameter('w', COMMAND_LINE_FLOAT, 1, "Wavelet weight" ,true, "0");

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
  boost::shared_ptr< hoNDArray<_real> > host_data = 
    read_nd_array<_real>((char*)parms.get_parameter('d')->get_string_value());

  if( !host_data.get() ){
    GINFO_STREAM(endl << "Input image not found. Quitting!\n" << endl);
    return 1;
  }
  
  if( host_data->get_number_of_dimensions() != 2 ){
    GINFO_STREAM(endl << "Input image is not two-dimensional. Quitting!\n" << endl);
    return 1;
  }
  
  // Upload host data to device
  cuNDArray<_real> data(*host_data);
  
  _real mu = (_real) parms.get_parameter('m')->get_float_value();
  _real lambda = (_real)parms.get_parameter('l')->get_float_value();

  if( mu <= (_real) 0.0 ) {
    GINFO_STREAM(endl << "Regularization parameter mu should be strictly positive. Quitting!\n" << endl);
    return 1;
  }

  unsigned int num_cg_iterations = parms.get_parameter('i')->get_int_value();
  unsigned int num_inner_iterations = parms.get_parameter('I')->get_int_value();
  unsigned int num_outer_iterations = parms.get_parameter('O')->get_int_value();
  
 // Define encoding operator (identity)
  boost::shared_ptr< identityOperator<cuNDArray<_real> > > E( new identityOperator<cuNDArray<_real> >() );
  E->set_weight( mu );
  E->set_domain_dimensions(data.get_dimensions());
  E->set_codomain_dimensions(data.get_dimensions());

  // Setup split-Bregman solver
  cuSbCgSolver<_real> sb;
  sb.set_encoding_operator( E );
  sb.set_max_outer_iterations(num_outer_iterations);
  sb.set_max_inner_iterations(num_inner_iterations);
  sb.set_output_mode( cuCgSolver<_real>::OUTPUT_VERBOSE );
   // Setup regularization operators

  if (lambda > 0){
  boost::shared_ptr< cuPartialDerivativeOperator<_real,2> > Rx( new cuPartialDerivativeOperator<_real,2>(0) );
  Rx->set_weight( lambda );
  Rx->set_domain_dimensions(data.get_dimensions());
  Rx->set_codomain_dimensions(data.get_dimensions());

  boost::shared_ptr< cuPartialDerivativeOperator<_real,2> > Ry( new cuPartialDerivativeOperator<_real,2>(1) );
  Ry->set_weight( lambda );
  Ry->set_domain_dimensions(data.get_dimensions());
  Ry->set_codomain_dimensions(data.get_dimensions());
  //sb.add_regularization_operator( Rx ); // Anisotropic denoising
  //sb.add_regularization_operator( Ry ); // Anisotropic denoising
  sb.add_regularization_group_operator( Rx ); // Isotropic denoising
  sb.add_regularization_group_operator( Ry); // Isotropic denoising
  sb.add_group();
  }
  
  _real wavelet = parms.get_parameter('w')->get_float_value();
  if (wavelet > 0){
	  auto dwt = boost::make_shared<cuDWTOperator<_real,2>>();
	  dwt->set_levels(3);
	  dwt->set_weight(wavelet);
	  sb.add_regularization_operator(dwt);
	  dwt->set_domain_dimensions(data.get_dimensions());
	  dwt->set_codomain_dimensions(data.get_dimensions());
	  dwt->use_random(true);
  }

  // Setup inner conjugate gradient solver
  sb.get_inner_solver()->set_max_iterations( num_cg_iterations );
  sb.get_inner_solver()->set_tc_tolerance( 1e-4 );
  sb.get_inner_solver()->set_output_mode( cuCgSolver<_real>::OUTPUT_WARNINGS );

  // Run split-Bregman solver
  boost::shared_ptr< cuNDArray<_real> > sbresult = sb.solve(&data);

  /*
  boost::shared_ptr< cuNDArray<_real> > sbresult(new cuNDArray<_real>(data.get_dimensions()));
  clear(sbresult.get());

  vector_td<float,4> daubechies4({0.6830127f,1.1830127f,0.3169873f,-0.1830127f});
  vector_td<float,2> haahr(1.0f,1.0f);
  vector_td<float,6> daubechies6{0.47046721f,1.14111692f,0.650365f,-0.19093442f, -0.12083221f,0.0498175f};

  cuDWTOperator<float,2> dwt;
  dwt.set_levels(3);
  dwt.mult_M(&data,sbresult.get());
  //data = *sbresult;
  shrink1(sbresult.get(),30.0f,&data);
  dwt.mult_MH(&data,sbresult.get());*/
  //clear(sbresult.get());
 // All done, write out the result
  boost::shared_ptr< hoNDArray<_real> > host_result = sbresult->to_host();
  write_nd_array<_real>(host_result.get(), (char*)parms.get_parameter('r')->get_string_value());
  
  return 0;
}
