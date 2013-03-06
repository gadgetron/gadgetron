/*
 * protonReconstruct.cpp
 *
 *  Created on: Dec 22, 2011
 *      Author: u051747
 */
#include <iostream>
#include "parameterparser.h"
#include "cuNDArray.h"
#include "hoCuNDArray.h"

#include "cuEncodedImageOperator.h"
#include "ndarray_vector_td_utilities.h"
#include "hoCuOperatorPathBackprojection.h"
#include "hoNDArray_fileio.h"
#include "check_CUDA.h"

#include "hoCuGPBBSolver.h"
#include "hoCuPartialDerivativeOperator.h"
#include "hdf5_utils.h"

#include "encodingOperatorContainer.h"
#include "ndarray_vector_td_utilities.h"
#include "hoCuNDArray_utils.h"
#include "imageOperator.h"
using namespace std;
using namespace Gadgetron;
typedef float _real;


int main( int argc, char** argv)
{

  //
  // Parse command line
  //
	cudaSetDevice(0);
	cudaDeviceReset();
  _real background =  0.00106;
  ParameterParser parms;
  parms.add_parameter( 'p', COMMAND_LINE_STRING, 1, "Input projection file name (.real)", true,"projections.real" );
  parms.add_parameter( 'w', COMMAND_LINE_STRING, 1, "Input uncertainties file name (.real)", false);
  parms.add_parameter( 't', COMMAND_LINE_STRING, 1, "Prior image file name (.real)", false);
  parms.add_parameter( 'v', COMMAND_LINE_FLOAT, 1, "Total Variation weight",false);
  parms.add_parameter( 'k', COMMAND_LINE_FLOAT,    1, "Prior Image TV ", false);
  parms.add_parameter( 'c', COMMAND_LINE_FLOAT,    1, "Prior image kappa ", true, "0" );
  parms.add_parameter( 's', COMMAND_LINE_STRING, 1, "Splines projection file name (.real)", true,"splines.real" );
  parms.add_parameter( 'f', COMMAND_LINE_STRING, 1, "Output image file name ", true,"image.hdf5" );


  parms.add_parameter( 'm', COMMAND_LINE_INT,    3, "x y and z image size in voxels ", true, "256 256 1" );
  parms.add_parameter( 'R', COMMAND_LINE_INT,    1, "Rescale directions ", false);

  parms.add_parameter( 'i', COMMAND_LINE_INT,    1, "Number of iterations", true, "10" );
  parms.add_parameter( 'e', COMMAND_LINE_FLOAT,    1, "Residual ", true, "1e-8" );

  parms.add_parameter( 'x', COMMAND_LINE_FLOAT,  1, "X size  (cm)", true, "1.0" );
  parms.add_parameter( 'y', COMMAND_LINE_FLOAT,  1, "Y size (cm)", true, "1.0" );
  parms.add_parameter( 'z', COMMAND_LINE_FLOAT,  1, "Z size (cm)", true, "1.0" );
  parms.add_parameter('O', COMMAND_LINE_FLOAT, 3, "X, Y and Z origin (cm)",true, "0 0 0 ");


  parms.parse_parameter_list(argc, argv);
  if( parms.all_required_parameters_set() ){
    cout << " Running reconstruction with the following parameters: " << endl;
    parms.print_parameter_list();
  }
  else{
    cout << " Some required parameters are missing: " << endl;
    parms.print_parameter_list();
    parms.print_usage();
    return 1;
  }


  boost::shared_ptr<hoCuNDArray<vector_td<_real,3> > > splines = boost::static_pointer_cast<hoCuNDArray<vector_td<_real,3> > >(read_nd_array< vector_td<_real,3> >((char*)parms.get_parameter('s')->get_string_value()));
  cout << "Number of spline elements: " << splines->get_number_of_elements() << endl;

  boost::shared_ptr< hoCuNDArray<_real> > projections = boost::static_pointer_cast<hoCuNDArray<_real> > (read_nd_array<_real >((char*)parms.get_parameter('p')->get_string_value()));

  std::cout << "Number of elements " << projections->get_number_of_elements() << std::endl;


  if (projections->get_number_of_elements() != splines->get_number_of_elements()/4){
	  cout << "Critical error: Splines and projections do not match dimensions" << endl;
	  return 0;
  }
  vector_td<_real,3> physical_dims;
  vector<unsigned int> ndims;
  ndims.push_back(3);


  physical_dims[0]= (_real) parms.get_parameter('x')->get_float_value();
  physical_dims[1]= (_real) parms.get_parameter('y')->get_float_value();
  physical_dims[2]= (_real) parms.get_parameter('z')->get_float_value();

  vector_td<_real,3> origin;


  origin[0]= (_real) parms.get_parameter('O')->get_float_value(0);
  origin[1]= (_real) parms.get_parameter('O')->get_float_value(1);
  origin[2]= (_real) parms.get_parameter('O')->get_float_value(2);


  hoCuGPBBSolver< _real> solver;

  solver.set_max_iterations( parms.get_parameter('i')->get_int_value() );
  //solver.set_tc_tolerance( (_real) parms.get_parameter('e')->get_float_value());
  //solver.set_alpha(1e-7);
  solver.set_output_mode( hoCuGPBBSolver< _real>::OUTPUT_VERBOSE );
   solver.set_non_negativity_constraint(true);
  boost::shared_ptr< hoCuOperatorPathBackprojection<_real> > E (new hoCuOperatorPathBackprojection<_real> );


  if (parms.get_parameter('R')->get_is_set()){
  	E->rescale_directions(parms.get_parameter('R')->get_int_value());
  }
  vector<unsigned int> rhs_dims;
  rhs_dims.push_back(parms.get_parameter('m')->get_int_value(0));
  rhs_dims.push_back(parms.get_parameter('m')->get_int_value(1));
  rhs_dims.push_back(parms.get_parameter('m')->get_int_value(2));
  cout << "RHS dims " << rhs_dims[0] << " " << rhs_dims[1] << " " << rhs_dims[2] << endl;


  boost::shared_ptr<hoCuNDArray<_real > > weights;
  boost::shared_ptr<hoCuNDArray<_real > > uncertainties;
  if (parms.get_parameter('w')->get_is_set()){
  	cout << "Weights currently not supported" << endl;
  	/*
	  boost::shared_ptr< hoCuNDArray<_real> > uncertainties = read_nd_array<_real >((char*)parms.get_parameter('w')->get_string_value());

	  _real* uncertainties_ptr = uncertainties->get_data_ptr();
	  for (int i =0; i < uncertainties->get_number_of_elements(); i++){
		  uncertainties_ptr[i] = 1/uncertainties_ptr[i];
	  }


	  //E->setup(splines,physical_dims,projections,weights,background);
	  //cuNDA_scale(weights.get(),projections.get());
	  
    	E->setup(splines,physical_dims,projections,background);
    	*/

  } else{

    	E->setup(splines,physical_dims,projections,origin,background);
   }



  E->set_domain_dimensions(&rhs_dims);
  E->set_codomain_dimensions(projections->get_dimensions().get());

  boost::shared_ptr<encodingOperatorContainer< hoCuNDArray<float> > > enc (new encodingOperatorContainer<hoCuNDArray<float> >());
  enc->set_domain_dimensions(&rhs_dims);
  enc->add_operator(E);

  boost::shared_ptr<hoCuNDArray<_real > > prior;
  if (parms.get_parameter('t')->get_is_set()){
 	  std::cout << "Prior image regularization in use" << std::endl;
   	  prior = boost::static_pointer_cast<hoCuNDArray<_real > >(read_nd_array<_real >((char*)parms.get_parameter('t')->get_string_value()));

   	  prior->reshape(&rhs_dims);
   	  _real offset = _real(0.01);
   	  //cuNDA_add(offset,prior.get());


   	  if (parms.get_parameter('c')->get_float_value()>0){

 		  boost::shared_ptr<imageOperator<hoCuNDArray<_real>,hoCuNDArray<_real> > > I (new imageOperator<hoCuNDArray<_real>,hoCuNDArray<_real> >());
		  I->compute(prior.get());

		  I->set_weight((_real) parms.get_parameter('c')->get_float_value());

		  I->set_codomain_dimensions(&rhs_dims);
		  I->set_domain_dimensions(&rhs_dims);
		  hoCuNDArray<_real> tmp = *prior;
		  

		  I->mult_M(prior.get(),&tmp);

		  //cuNDA_scal(I->get_weight(),&tmp);
		  std::vector<hoCuNDArray<_real>* > proj;
		  proj.push_back(projections.get());
		  proj.push_back(&tmp);
		  enc->add_operator(I);
		  projections = enc->create_codomain(proj);
   	  }
  }

  if (parms.get_parameter('v')->get_is_set()){



    boost::shared_ptr< hoCuPartialDerivativeOperator<_real,3> > Rx( new hoCuPartialDerivativeOperator<_real,3>(0) );
    boost::shared_ptr< hoCuPartialDerivativeOperator<_real,3> > Ry( new hoCuPartialDerivativeOperator<_real,3>(1) );
    Rx->set_codomain_dimensions(&rhs_dims);
    Ry->set_codomain_dimensions(&rhs_dims);
    Rx->set_domain_dimensions(&rhs_dims);
    Ry->set_domain_dimensions(&rhs_dims);
    Rx->set_weight(parms.get_parameter('v')->get_float_value());
    Ry->set_weight(parms.get_parameter('v')->get_float_value());

    solver.add_regularization_group_operator(Rx);
    solver.add_regularization_group_operator(Ry);
    solver.add_group(1);

		  
  }


  
  solver.set_encoding_operator(enc);
	if (parms.get_parameter('t')->get_is_set()){
	 solver.set_x0(prior);
	}
	//hoCuNDA_clear(projections.get());
	//CHECK_FOR_CUDA_ERROR();

	float res = dot(projections.get(),projections.get());

	boost::shared_ptr< hoCuNDArray<_real> > result = solver.solve(projections.get());

	//write_nd_array<_real>(result.get(), (char*)parms.get_parameter('f')->get_string_value());
	std::stringstream ss;
	for (int i = 0; i < argc; i++){
		ss << argv[i] << " ";
	}
	saveNDArray2HDF5<3>(result.get(),parms.get_parameter('f')->get_string_value(),physical_dims,origin,ss.str(), solver.get_max_iterations());
}



