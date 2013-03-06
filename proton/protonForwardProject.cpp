/*
 * protonForwardProject.cpp
 *
 *  Created on: Jan 3, 2012
 *      Author: u051747
 */


#include <iostream>
#include "parameterparser.h"
#include "cuNDArray.h"
#include "cuCgSolver.h"
#include "ndarray_vector_td_utilities.h"
#include "cuOperatorPathBackprojection.h"
#include "hoNDArray_fileio.h"
#include <vector>
using namespace std;
typedef float _real;

using namespace Gadgetron;


int main( int argc, char** argv)
{
	_real background =  0.00106;
  //
  // Parse command line
  //

  ParameterParser parms;
  parms.add_parameter( 'p', COMMAND_LINE_STRING, 1, "Input image file name (.real)", true,"phantom.real" );
  parms.add_parameter( 's', COMMAND_LINE_STRING, 1, "Splines projection file name (.real)", true,"splines.real" );
  parms.add_parameter( 'f', COMMAND_LINE_STRING, 1, "Output projection file name ", true,"projections.real" );


  parms.add_parameter( 'x', COMMAND_LINE_FLOAT,  1, "X size  (cm)", true, "1.0" );
  parms.add_parameter( 'y', COMMAND_LINE_FLOAT,  1, "Y size (cm)", true, "1.0" );
  parms.add_parameter( 'z', COMMAND_LINE_FLOAT,  1, "Z size (cm)", true, "1.0" );

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


  boost::shared_ptr<hoNDArray<vector_td<_real,3> > > host_splines = read_nd_array< vector_td<_real,3> >((char*)parms.get_parameter('s')->get_string_value());
  boost::shared_ptr<cuNDArray<vector_td<_real,3> > > splines (new cuNDArray< vector_td<_real,3> >(host_splines.get()));

  boost::shared_ptr< hoNDArray<_real> > host_phantom = read_nd_array<_real >((char*)parms.get_parameter('p')->get_string_value());
  boost::shared_ptr<cuNDArray<_real > > phantom (new cuNDArray<_real >(host_phantom.get()));


  vector_td<_real,3> physical_dims;
  vector<unsigned int> ndims;
  ndims.push_back(3);


  physical_dims.vec[0]= (_real) parms.get_parameter('x')->get_float_value();
  physical_dims.vec[1]= (_real) parms.get_parameter('y')->get_float_value();
  physical_dims.vec[2]= (_real) parms.get_parameter('z')->get_float_value();




  boost::shared_ptr< cuOperatorPathBackprojection<_real> > E (new cuOperatorPathBackprojection<_real> );
  cout << "Performing setup" << endl;

  //

  boost::shared_ptr<cuNDArray<_real> > projections = boost::shared_ptr<cuNDArray<_real> >(new cuNDArray<_real>);
  vector<unsigned int> projection_dims;
  projection_dims.push_back(splines->get_dimensions()->at(0)/4);

  projections->create(&projection_dims);
  projections->clear();
  E->setup(splines,physical_dims,projections,background);
  projections->abs();
  cout << "Starting forward projection" << endl;
  E->mult_M(phantom.get(),projections.get());


   boost::shared_ptr< hoNDArray<_real> > host_result = projections->to_host();
   write_nd_array<_real>(host_result.get(), (char*)parms.get_parameter('f')->get_string_value());
   cout <<" AAAAND, done" << endl;
}






