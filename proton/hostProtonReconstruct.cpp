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


#include "hoCuOperatorPathBackprojection.h"
#include "hoNDArray_fileio.h"
#include "check_CUDA.h"

#include "hoCuGPBBSolver.h"
#include "hoCuPartialDerivativeOperator.h"
#include "hdf5_utils.h"

#include "encodingOperatorContainer.h"
#include "hoCuOperator.h"
#include "hoImageOperator.h"
#include "identityOperator.h"
#include <boost/program_options.hpp>
#include "vector_td_io.h"

#include "hoCuTVOperator.h"
#include "hoCuTVpicsOperator.h"
#include "projectionSpaceOperator.h"
#include "weightingOperator.h"
#include "ABOCSSolver.h"
using namespace std;
using namespace Gadgetron;
typedef float _real;


namespace po = boost::program_options;
int main( int argc, char** argv)
{

  //
  // Parse command line
  //

  _real background =  0.00106;
  std::string projectionsName;
	std::string splinesName;
	std::string outputFile;
	vector_td<int,3> dimensions;
	vector_td<float,3> physical_dims;
	vector_td<float,3> origin;
	int iterations;
	int device;
	po::options_description desc("Allowed options");
	desc.add_options()
			("help", "produce help message")
			("projections,p", po::value<std::string>(&projectionsName)->default_value("projections.real"), "File containing the projection data")
			("splines,s", po::value<std::string>(&splinesName)->default_value("splines.real"), "File containing the spline trajectories")
			("dimensions,d", po::value<vector_td<int,3> >(&dimensions)->default_value(vector_td<int,3>(512,512,1)), "Pixel dimensions of the image")
			("size,S", po::value<vector_td<float,3> >(&physical_dims)->default_value(vector_td<float,3>(20,20,5)), "Dimensions of the image")
			("center,c", po::value<vector_td<float,3> >(&origin)->default_value(vector_td<float,3>(0,0,0)), "Center of the reconstruction")
			("iterations,i", po::value<int>(&iterations)->default_value(10), "Dimensions of the image")
			("output,f", po::value<std::string>(&outputFile)->default_value("image.hdf5"), "Output filename")
			("prior,P", po::value<std::string>(),"Prior image filename")
			("prior-weight,k",po::value<float>(),"Weight of the prior image")
			("wprior-weight,w",po::value<float>(),"Weight of the weighted prior image")
			("PICS-weight,C",po::value<float>(),"TV Weight of the prior image (Prior image compressed sensing)")
			("TV,T",po::value<float>(),"TV Weight ")
			("device",po::value<int>(&device)->default_value(0),"Number of the device to use (0 indexed)")

	;

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);

	if (vm.count("help")) {
		cout << desc << "\n";
		return 1;
	}
	std::cout << "Command line options:" << std::endl;
	for (po::variables_map::iterator it = vm.begin(); it != vm.end(); ++it){
		boost::any a = it->second.value();
		std::cout << it->first << ": ";
		if (a.type() == typeid(std::string)) std::cout << it->second.as<std::string>();
		else if (a.type() == typeid(int)) std::cout << it->second.as<int>();
		else if (a.type() == typeid(float)) std::cout << it->second.as<float>();
		else if (a.type() == typeid(vector_td<float,3>)) std::cout << it->second.as<vector_td<float,3> >();
		else if (a.type() == typeid(vector_td<int,3>)) std::cout << it->second.as<vector_td<int,3> >();
		else std::cout << "Unknown type" << std::endl;
		std::cout << std::endl;
	}
	cudaSetDevice(device);
	cudaDeviceReset();


  boost::shared_ptr<hoCuNDArray<vector_td<_real,3> > > splines = boost::static_pointer_cast<hoCuNDArray<vector_td<_real,3> > >(read_nd_array< vector_td<_real,3> >(splinesName.c_str()));
  cout << "Number of spline elements: " << splines->get_number_of_elements() << endl;

  boost::shared_ptr< hoCuNDArray<_real> > projections = boost::static_pointer_cast<hoCuNDArray<_real> > (read_nd_array<_real >(projectionsName.c_str()));

  std::cout << "Number of elements " << projections->get_number_of_elements() << std::endl;


  if (projections->get_number_of_elements() != splines->get_number_of_elements()/4){
	  cout << "Critical error: Splines and projections do not match dimensions" << endl;
	  return 0;
  }

  ABOCSSolver< hoCuGPBBSolver< _real> > solver;
  solver.set_eps(_real(7.6e6));

  solver.set_max_iterations( iterations);
  //solver.set_tc_tolerance( (_real) parms.get_parameter('e')->get_float_value());
  //solver.set_alpha(1e-7);
  solver.set_output_mode( hoCuGPBBSolver< _real>::OUTPUT_VERBOSE );
  solver.set_non_negativity_constraint(true);
  boost::shared_ptr< hoCuOperatorPathBackprojection<_real> > E (new hoCuOperatorPathBackprojection<_real> );
  E->setup(splines,physical_dims,projections,origin,background);

  std::vector<unsigned int> rhs_dims(&dimensions[0],&dimensions[3]); //Quick and dirty vector_td to vector
  E->set_domain_dimensions(&rhs_dims);
  E->set_codomain_dimensions(projections->get_dimensions().get());



  boost::shared_ptr<hoCuNDArray<_real > > prior;
  if (vm.count("prior")){
 	  std::cout << "Prior image used as initial guess" << std::endl;
		prior = boost::static_pointer_cast<hoCuNDArray<_real > >(read_nd_array<_real >(vm["prior"].as<std::string>().c_str()));

		prior->reshape(&rhs_dims);
		_real offset = _real(0.01);
		//cuNDA_add(offset,prior.get());


		if (vm.count("prior-weight")){
			std::cout << "Prior image difference regularization in use" << std::endl;
			boost::shared_ptr<identityOperator<hoCuNDArray<_real> > > Itmp ( new identityOperator<hoCuNDArray<_real> >);

			boost::shared_ptr<projectionSpaceOperator<hoCuNDArray<_real> > > I (new projectionSpaceOperator<hoCuNDArray<_real> >(Itmp));

			I->set_weight(vm["prior-weight"].as<float>());

			I->set_codomain_dimensions(&rhs_dims);
			I->set_domain_dimensions(&rhs_dims);
			I->set_projections(prior);
			//hoCuNDArray<_real> tmp = *prior;


			//I->mult_M(prior.get(),&tmp);
			solver.add_regularization_operator(I);
		}

		if (vm.count("wprior-weight")){
			std::cout << "Prior image difference regularization in use" << std::endl;

			boost::shared_ptr<weightingOperator<hoCuNDArray<_real> > > Wtmp ( new weightingOperator<hoCuNDArray<_real> >());



			boost::shared_ptr<hoCuNDArray<_real> > weights(new hoCuNDArray<_real>(*prior));

			abs_inplace(weights.get());
			*weights += _real(0.01);

			reciprocal_inplace(weights.get());

			Wtmp->set_weights(weights);

			boost::shared_ptr<projectionSpaceOperator<hoCuNDArray<_real> > > W (new projectionSpaceOperator<hoCuNDArray<_real> >(Wtmp));
			W->set_weight(vm["wprior-weight"].as<float>());

			W->set_codomain_dimensions(&rhs_dims);
			W->set_domain_dimensions(&rhs_dims);

			boost::shared_ptr<hoCuNDArray<_real> > tmp_array(new hoCuNDArray<_real>(*prior));
			Wtmp->mult_M(prior.get(),tmp_array.get());

			W->set_projections(tmp_array);
			//hoCuNDArray<_real> tmp = *prior;


			//I->mult_M(prior.get(),&tmp);
			solver.add_regularization_operator(W);
		}
		if (vm.count("PICS-weight")){
			std::cout << "PICS in used" << std::endl;
			boost::shared_ptr<hoCuTVpicsOperator<float,3> > pics (new hoCuTVpicsOperator<float,3>);

			pics->set_prior(prior);
			pics->set_weight(vm["PICS-weight"].as<float>());
			solver.add_nonlinear_operator(pics);



		}
		solver.set_x0(prior);
  }

  if (vm.count("TV")){
	  std::cout << "Total variation regularization in use" << std::endl;
	  boost::shared_ptr<hoCuTVOperator<float,3> > tv(new hoCuTVOperator<float,3>);
	  tv->set_weight(vm["TV"].as<float>());
	  solver.add_nonlinear_operator(tv);
  }

  solver.set_encoding_operator(E);

	//hoCuNDA_clear(projections.get());
	//CHECK_FOR_CUDA_ERROR();

	//float res = dot(projections.get(),projections.get());

	boost::shared_ptr< hoCuNDArray<_real> > result = solver.solve(projections.get());

	//write_nd_array<_real>(result.get(), (char*)parms.get_parameter('f')->get_string_value());
	std::stringstream ss;
	for (int i = 0; i < argc; i++){
		ss << argv[i] << " ";
	}
	saveNDArray2HDF5<3>(result.get(),outputFile,physical_dims,origin,ss.str(), solver.get_max_iterations());
}



