/*
 * protonReconstruct.cpp
 *
 *  Created on: Dec 22, 2011
 *      Author: u051747
 */
#include <iostream>
#include <boost/program_options.hpp>
#include "cuNDArray.h"
#include "cuCgSolver.h"
#include "cuImageOperator.h"
#include "cuEncodedImageOperator.h"
#include "ndarray_vector_td_utilities.h"
#include "cuOperatorPathBackprojection.h"
#include "cuPartialDerivativeOperator.h"
#include "cuLaplaceOperator.h"
#include "hoNDArray_fileio.h"
#include "check_CUDA.h"
#include "cuLwSolver.h"
#include "cuGPBBSolver.h"
#include "identityOperator.h"


#include <sstream>
#include "hdf5_utils.h"

#include "encodingOperatorContainer.h"
#include "cuSARTSolver.h"
#include "cuMLSolver.h"



#include "ndarray_vector_td_utilities.h"
using namespace std;
using namespace Gadgetron;
typedef float _real;
/*
boost::shared_ptr< cuNDArray<_real> >  recursiveSolver(cuNDArray<_real> * rhs,cuCGSolver<_real, _real> * cg,int depth){
	std::cout << "Recursion depth " << depth << " reached" << std::endl;
	if (depth == 0){
		std::cout << "Running solver for depth " << depth  << std::endl;
		return cg->solve(rhs);
	} else {
		boost::shared_ptr< cuNDArray<_real> > rhs_temp = cuNDA_downsample<_real,2>(rhs);

		boost::shared_ptr< cuNDArray<_real> > guess = recursiveSolver(rhs_temp.get(),cg,depth-1);
		guess = cuNDA_upsample<_real,2>(guess.get());
		std::cout << "Running solver for depth " << depth  << std::endl;
		return cg->solve(rhs,guess.get());

	}

}
*/

/*
template<class T> void notify(T val){
	std::cout << val <<std::endl;
}
*/

template<class T> void notify(std::string val){
	std::cout << val << std::endl;
}

namespace po = boost::program_options;
int main( int argc, char** argv)
{

  //
  // Parse command line
  //
  _real background =  0.00106;
  /*
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
  parms.add_parameter( 'n', COMMAND_LINE_INT,    1, "z image size ", true, "1" );

  parms.add_parameter( 'i', COMMAND_LINE_INT,    1, "Number of iterations", true, "10" );
  parms.add_parameter( 'e', COMMAND_LINE_FLOAT,    1, "Residual ", true, "1e-8" );

  parms.add_parameter( 'x', COMMAND_LINE_FLOAT,  1, "X size  (cm)", true, "1.0" );
  parms.add_parameter( 'y', COMMAND_LINE_FLOAT,  1, "Y size (cm)", true, "1.0" );
  parms.add_parameter( 'z', COMMAND_LINE_FLOAT,  1, "Z size (cm)", true, "1.0" );

  parms.add_parameter( 'd', COMMAND_LINE_INT,    1, "Multiscale depth ", true, "1" );
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
*/

  std::string projectionsName;
  std::string splinesName;
  std::string outputFile;
  vector_td<int,3> dimensions;
  vector_td<float,3> physical_dims;
  vector_td<float,3> origin;
  int iterations;
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

  ;


  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  std::cout << "Command line options:" << std::endl;
  for (po::variables_map::iterator it = vm.begin(); it != vm.end(); ++it){
  	boost::any a = it->second.value();
  	std::cout << it->first << " ";
  	if (a.type() == typeid(std::string)) std::cout << it->second.as<std::string>();
  	if (a.type() == typeid(int)) std::cout << it->second.as<int>();
  	if (a.type() == typeid(vector_td<float,3>)) std::cout << it->second.as<vector_td<float,3> >();
  	if (a.type() == typeid(vector_td<int,3>)) std::cout << it->second.as<vector_td<int,3> >();
  	std::cout << std::endl;

  }
  std::cout <<  std::endl;
	boost::shared_ptr<hoNDArray<vector_td<_real,3> > > host_splines = read_nd_array< vector_td<_real,3> >(splinesName.c_str());
  boost::shared_ptr<cuNDArray<vector_td<_real,3> > > splines (new cuNDArray< vector_td<_real,3> >(host_splines.get()));
  cout << "Number of spline elements: " << splines->get_number_of_elements() << endl;

  boost::shared_ptr< hoNDArray<_real> > host_projections = read_nd_array<_real >(projectionsName.c_str());
  boost::shared_ptr<cuNDArray<_real > > projections (new cuNDArray<_real >(host_projections.get()));
  boost::shared_ptr<cuNDArray<_real > > projections_old = projections;

  std::cout << "Number of elements " << projections->get_number_of_elements() << std::endl;

  cout << "Number of projection elements: " << projections->get_number_of_elements() << endl;
  if (projections->get_number_of_elements() != splines->get_number_of_elements()/4){
	  cout << "Critical error: Splines and projections do not match dimensions" << endl;
	  return 0;
  }

  vector<unsigned int> ndims;
  ndims.push_back(3);



  cuGPBBSolver<_real> cg;

  cg.set_max_iterations( iterations);

  cg.set_output_mode( cuCgSolver<_real>::OUTPUT_VERBOSE );
   cg.set_non_negativity_constraint(true);
  boost::shared_ptr< cuOperatorPathBackprojection<_real> > E (new cuOperatorPathBackprojection<_real> );


  vector<unsigned int> rhs_dims;


	rhs_dims.push_back(dimensions[0]);
	rhs_dims.push_back(dimensions[1]);
	rhs_dims.push_back(dimensions[2]);
  cout << "RHS dims " << rhs_dims[0] << " " << rhs_dims[1] << " " << rhs_dims[2] << endl;


  /*boost::shared_ptr<cuNDArray<_real> > guess(new cuNDArray<_real>);
  vector<unsigned int> guess_dims  = rhs_dims;
  guess_dims[0] /=multiScale_depth;
  guess_dims[1] /=multiScale_depth;
  guess_dims[2] /=multiScale_depth;
  guess->create(&guess_dims);
  cuNDA_clear(guess.get(),_real(1));
   */
  boost::shared_ptr<cuNDArray<_real > > weights;
  boost::shared_ptr<cuNDArray<_real > > uncertainties;
  E->setup(splines,physical_dims,projections,origin,background);




  E->set_domain_dimensions(&rhs_dims);
  E->set_codomain_dimensions(projections->get_dimensions().get());
  boost::shared_ptr<encodingOperatorContainer<cuNDArray<_real> > > enc (new encodingOperatorContainer<cuNDArray<_real> >());
  enc->set_domain_dimensions(&rhs_dims);
  enc->add_operator(E);

  boost::shared_ptr<cuNDArray<_real > > prior;
  /*
  if (parms.get_parameter('t')->get_is_set()){
 	  std::cout << "Prior image regularization in use" << std::endl;
   	  boost::shared_ptr< hoNDArray<_real> > host_prior = read_nd_array<_real >((char*)parms.get_parameter('t')->get_string_value());
   	  prior = boost::shared_ptr<cuNDArray<_real > > (new cuNDArray<_real >(host_prior.get()));
	  prior->reshape(&rhs_dims);
   	  _real offset = _real(0.01);
   	  //cuNDA_add(offset,prior.get());


   	  if (parms.get_parameter('c')->get_float_value()>0){

 		  boost::shared_ptr<cuImageOperator<_real> > I (new cuImageOperator<_real>());
		  I->compute(prior.get());
	    //	    boost::shared_ptr<cuIdentityOperator<_real> > I (new cuIdentityOperator<_real>());
		  I->set_weight((_real) parms.get_parameter('c')->get_float_value());
		  //cuNDA_threshold_min(0.1f,I->get());
		  //boost::shared_ptr< hoNDArray<_real> > host_dxy = I->get()->to_host();
		  //write_nd_array<_real>(host_dxy.get(), "priorInverse.real");
		  I->set_codomain_dimensions(&rhs_dims);
		  I->set_domain_dimensions(&rhs_dims);
		  cuNDArray<_real> tmp = *prior;
		  

		  I->mult_M(prior.get(),&tmp);

		  //cuNDA_scal(I->get_weight(),&tmp);
		  std::vector<cuNDArray<_real>* > proj;
		  proj.push_back(projections.get());
		  proj.push_back(&tmp);
		  enc->add_operator(I);
		  projections = enc->create_codomain(proj);
   	  }
  }

  if (parms.get_parameter('v')->get_is_set()){



    boost::shared_ptr< cuPartialDerivativeOperator<_real,3> > Rx( new cuPartialDerivativeOperator<_real,3>(0) );
    boost::shared_ptr< cuPartialDerivativeOperator<_real,3> > Ry( new cuPartialDerivativeOperator<_real,3>(1) );
    Rx->set_codomain_dimensions(&rhs_dims);
    Ry->set_codomain_dimensions(&rhs_dims);
    Rx->set_domain_dimensions(&rhs_dims);
    Ry->set_domain_dimensions(&rhs_dims);
    Rx->set_weight(parms.get_parameter('v')->get_float_value());
    Ry->set_weight(parms.get_parameter('v')->get_float_value());

    cg.add_regularization_group_operator(Rx);
    cg.add_regularization_group_operator(Ry);
    cg.add_group(1);


		  
  }

*/



  boost::shared_ptr< cuNDArray<_real> > cgresult;

  

  cg.set_encoding_operator(enc);
  /*
     if (parms.get_parameter('t')->get_is_set()){
    	 prior->reshape(&rhs_dims);
	 cgresult=prior;
	 cg.set_x0(cgresult);
     } else {
       cgresult = boost::shared_ptr< cuNDArray<_real> >(new cuNDArray<_real>());
       cgresult->create(&rhs_dims);
       cgresult->clear();
     }
     */
     //cgresult = sb.solve(projections.get());

     //cgresult=sart.solve(projections.get());
     //cgresult=ml.solve(projections.get());
     cgresult = cg.solve(projections.get());

     cuNDArray<_real> tp = *projections_old;

     E->mult_M(cgresult.get(),&tp);
     axpy(-1.0f,projections_old.get(),&tp);
     std::cout << "Total residual " << nrm2(&tp) << std::endl;
     //cgresult = recursiveSolver(&rhs,&cg,multiScale_depth-1);
    /* boost::shared_ptr< vector<cuNDArray<_real> *> > projVec = boost::shared_ptr< vector<cuNDArray<_real> *> >(new vector<cuNDArray<_real> *>);
     projVec->push_back(projections.get());
     projVec->push_back(0);
     cgresult = cg.solve(projVec,&rhs_dims);*/
/*
     for (int i = 0; i < 30; i++) {
       cg.set_x0(cgresult);
	 cgresult = cg.solve(projections.get());

	 cuNDArray<_real> tmp_proj = *projections_old;

	 E->mult_M(cgresult.get(),&tmp_proj);
	 cuNDA_axpy(-1.0f,projections_old.get(),&tmp_proj);

	 cuNDA_threshold_amin(uncertainties.get(),&tmp_proj);
	 std::cout << "Iteration " << i << " residual: " << cuNDA_dot(&tmp_proj,&tmp_proj) << std::endl;

	 cuNDArray<_real> proj_view;
	 boost::shared_ptr<std::vector<unsigned int> > projection_dims = projections_old->get_dimensions();
	 proj_view.create(projection_dims.get(),projections->get_data_ptr());

	 cuNDA_axpy(-1.0f,&tmp_proj,&proj_view);

	 printf("Dumping frame\n");
	 char filename_real[256];
	 sprintf(&(filename_real[0]), "out_u_%04i.real", i);
	 boost::shared_ptr<hoNDArray<_real> > tmp = cgresult->to_host();
	 write_nd_array<_real>(tmp.get(), filename_real);
       }
	  
     */

   boost::shared_ptr< hoNDArray<_real> > host_result = cgresult->to_host();
   //write_nd_array<_real>(host_result.get(), (char*)parms.get_parameter('f')->get_string_value());
   std::stringstream ss;
   	for (int i = 0; i < argc; i++){
   		ss << argv[i] << " ";
   	}

   	saveNDArray2HDF5<3>(host_result.get(),outputFile,physical_dims,origin,ss.str(), cg.get_max_iterations());

}



