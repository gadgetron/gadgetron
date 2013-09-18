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


#include "hoCuProtonSubsetOperator.h"
#include "hoNDArray_fileio.h"
#include "check_CUDA.h"

#include "hoCuGPBBSolver.h"
#include "hoCuPartialDerivativeOperator.h"
#include "hoCuNDArray_blas.h"
#include "hoCuNDArray_operators.h"

#include "osSARTSolver.h"
#include "hoOSGPBBSolver.h"
#include "hoOSCGSolver.h"
#include "hoOSCGPBBSolver.h"
//#include "hoABIBBSolver.h"
#include "hdf5_utils.h"

#include "encodingOperatorContainer.h"
#include "hoCuOperator.h"
#include "hoImageOperator.h"
#include "identityOperator.h"
#include <boost/program_options.hpp>
#include "vector_td_io.h"

#include "hoCuTvOperator.h"
#include "hoCuTvPicsOperator.h"
#include "projectionSpaceOperator.h"

using namespace std;
using namespace Gadgetron;
typedef float _real;

typedef solver<hoCuNDArray<_real>, hoCuNDArray<_real> > baseSolver;
namespace po = boost::program_options;
int main( int argc, char** argv)
{

  //
  // Parse command line
  //

  _real background =  0.00106;
  std::string dataName;
	std::string outputFile;
	vector_td<int,3> dimensions;
	vector_td<float,3> physical_dims;
	vector_td<float,3> origin;
	int iterations;
	int device;
	int subsets;
	po::options_description desc("Allowed options");
	desc.add_options()
			("help", "produce help message")
			("data,D", po::value<std::string>(&dataName)->default_value("data.hdf5"), "HDF5 file containing projections and splines")
			("dimensions,d", po::value<vector_td<int,3> >(&dimensions)->default_value(vector_td<int,3>(512,512,1)), "Pixel dimensions of the image")
			("size,S", po::value<vector_td<float,3> >(&physical_dims)->default_value(vector_td<float,3>(20,20,5)), "Dimensions of the image")
			("center,c", po::value<vector_td<float,3> >(&origin)->default_value(vector_td<float,3>(0,0,0)), "Center of the reconstruction")
			("iterations,i", po::value<int>(&iterations)->default_value(10), "Dimensions of the image")
			("output,f", po::value<std::string>(&outputFile)->default_value("image.hdf5"), "Output filename")
			("prior,P", po::value<std::string>(),"Prior image filename")
			("prior-weight,k",po::value<float>(),"Weight of the prior image")
			("device",po::value<int>(&device)->default_value(0),"Number of the device to use (0 indexed)")
			("TV,T",po::value<float>(),"TV Weight ")
			("subsets,n", po::value<int>(&subsets)->default_value(10), "Number of subsets to use")
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
	//cudaDeviceReset();


/*
  hoCuGPBBSolver< _real> solver;

  solver.set_max_iterations( iterations);

  solver.set_output_mode( hoCuGPBBSolver< _real>::OUTPUT_VERBOSE );*/
	//osSARTSolver<hoCuNDArray<_real> > solver;
	//hoOSCGPBBSolver<hoCuNDArray<_real> > solver;
	//hoABIBBSolver<hoCuNDArray<_real> > solver;
	hoOSCGSolver<hoCuNDArray<_real> > solver;
  solver.set_non_negativity_constraint(true);
  solver.set_max_iterations(iterations);
  boost::shared_ptr< hoCuProtonSubsetOperator<_real> > E (new hoCuProtonSubsetOperator<_real>(subsets) );
  boost::shared_ptr<hoCuNDArray<_real> >  projections = E->load_data(dataName,physical_dims,origin,background);

  std::vector<unsigned int> rhs_dims(&dimensions[0],&dimensions[3]); //Quick and dirty vector_td to vector
  E->set_domain_dimensions(&rhs_dims);
  E->set_codomain_dimensions(projections->get_dimensions().get());

  boost::shared_ptr<hoCuNDArray<_real > > prior;
  if (vm.count("prior")){
 	  std::cout << "Prior image regularization in use" << std::endl;
		prior = boost::static_pointer_cast<hoCuNDArray<_real > >(read_nd_array<_real >(vm["prior"].as<std::string>().c_str()));

		prior->reshape(&rhs_dims);
		_real offset = _real(0.01);
		//cuNDA_add(offset,prior.get());


		if (vm.count("prior-weight")){

		//boost::shared_ptr<hoImageOperator<_real> > I (new hoImageOperator<_real > ());
		//I->compute(prior.get());
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
		solver.set_x0(prior);
  }
  if (vm.count("TV")){
	  std::cout << "Total variation regularization in use" << std::endl;
	  boost::shared_ptr<hoCuTvOperator<float,3> > tv(new hoCuTvOperator<float,3>);
	  tv->set_weight(vm["TV"].as<float>());
	  solver.add_nonlinear_operator(tv);
  }

  solver.set_encoding_operator(E);
  solver.set_output_mode(baseSolver::OUTPUT_VERBOSE);
	//hoCuNDA_clear(projections.get());
	//CHECK_FOR_CUDA_ERROR();


	//float res = dot(projections.get(),projections.get());

	boost::shared_ptr< hoCuNDArray<_real> > result = solver.solve(projections.get());

	hoCuNDArray<_real> tmp_proj(projections->get_dimensions());
	E->mult_M(result.get(),&tmp_proj,false);
	tmp_proj -= *projections;

	std::cout << "L2 norm of residual: " << nrm2(&tmp_proj) << std::endl;
	//write_nd_array<_real>(result.get(), (char*)parms.get_parameter('f')->get_string_value());
	std::stringstream ss;
	for (int i = 0; i < argc; i++){
		ss << argv[i] << " ";
	}
	saveNDArray2HDF5<3>(result.get(),outputFile,physical_dims,origin,ss.str(), solver.get_max_iterations());
}



