#include "hoCuNDArray_utils.h"
#include "radial_utilities.h"
#include "hoNDArray_fileio.h"
#include "cuNDArray.h"
#include "imageOperator.h"
#include "identityOperator.h"
#include "hoPartialDerivativeOperator.h"
#include "hoCudaConebeamProjectionOperator.h"
#include "cuConvolutionOperator.h"
#include "hoCuNDArray_blas.h"
#include "hoCuNDArray_operators.h"
#include "hoCuNDArray_blas.h"
#include "cgSolver.h"
#include "CBCT_acquisition.h"
#include "complext.h"
#include "encodingOperatorContainer.h"
#include "vector_td_io.h"
#include "hoPartialDerivativeOperator.h"
#include "hoGpBbSolver.h"
#include "hoCuTvOperator.h"
#include "hoCuTvPicsOperator.h"
#include "hoCuNCGSolver.h"
#include <iostream>
#include <algorithm>
#include <sstream>
#include <math_constants.h>
#include <boost/program_options.hpp>

using namespace std;
using namespace Gadgetron;

namespace po = boost::program_options;

int main(int argc, char** argv) 
{  
  string acquisition_filename;
  string outputFile;  
  uintd3 imageSize;
  floatd3 voxelSize;
  int device;
  unsigned int iterations;
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help", "produce help message")
    ("acquisition,a", po::value<string>(&acquisition_filename)->default_value("acquisition.hdf5"), "Acquisition data")
    ("samples,n",po::value<unsigned int>(),"Number of samples per ray")
    ("output,f", po::value<string>(&outputFile)->default_value("reconstruction.real"), "Output filename")
    ("size,s",po::value<uintd3>(&imageSize)->default_value(uintd3(512,512,1)),"Image size in pixels")
    ("binning,b",po::value<string>(),"Binning file for 4d reconstruction")
    ("SAG","Use exact SAG correction if present")
    ("voxelSize,v",po::value<floatd3>(&voxelSize)->default_value(floatd3(0.488f,0.488f,1.0f)),"Voxel size in mm")
    ("dimensions,d",po::value<floatd3>(),"Image dimensions in mm. Overwrites voxelSize.")
    ("iterations,i",po::value<unsigned int>(&iterations)->default_value(10),"Number of iterations")
    ("TV,T",po::value<float>(),"TV Weight ")
    ("prior", po::value<std::string>(),"Prior image filename")
    ("PICS",po::value<float>(),"TV Weight of the prior image (Prior image compressed sensing)")
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
    else if (a.type() == typeid(unsigned int)) std::cout << it->second.as<unsigned int>();
    else if (a.type() == typeid(float)) std::cout << it->second.as<float>();
    else if (a.type() == typeid(vector_td<float,3>)) std::cout << it->second.as<vector_td<float,3> >();
    else if (a.type() == typeid(vector_td<int,3>)) std::cout << it->second.as<vector_td<int,3> >();
    else if (a.type() == typeid(vector_td<unsigned int,3>)) std::cout << it->second.as<vector_td<unsigned int,3> >();
    else std::cout << "Unknown type" << std::endl;
    std::cout << std::endl;
  }
  cudaSetDevice(device);
  cudaDeviceReset();

  //Really weird stuff. Needed to initialize the device?? Should find real bug.
  cudaDeviceManager::Instance()->lockHandle();
  cudaDeviceManager::Instance()->unlockHandle();
       
  boost::shared_ptr<CBCT_acquisition> ps(new CBCT_acquisition());
  ps->load(acquisition_filename);
  unsigned int ip_w = ps->get_projections()->get_size(0);
  unsigned int ip_h = ps->get_projections()->get_size(1);
  uintd2 ps_dims_in_pixels(ip_w, ip_h);

  ps->get_geometry()->print(std::cout);

  //ps->get_geometry()->set_spacing(floatd2(ip_sx,ip_sy));
  float SDD = ps->get_geometry()->get_SDD();
  float SAD = ps->get_geometry()->get_SAD();
  // Projection space physical dimensions 40cm x 30cm

  boost::shared_ptr<CBCT_binning> binning(new CBCT_binning());
  if (vm.count("binning")){
    std::cout << "Loading binning data" << std::endl;
    binning->load(vm["binning"].as<string>());
	  
  } else binning->set_as_default_3d_bin(ps->get_projections()->get_size(2));
	
  binning->print(std::cout);

  floatd3 imageDimensions;
  if (vm.count("dimensions")){
    imageDimensions = vm["dimensions"].as<floatd3>();
    voxelSize = imageDimensions/imageSize;
  }
  else imageDimensions = voxelSize*imageSize;


  //floatd3 is_dims_in_mm = is_dims_in_pixels * is_spacing_in_mm;

  float lengthOfRay_in_mm = norm(imageDimensions);
  unsigned int numSamplesPerPixel = 3;
  float minSpacing = min(voxelSize)/numSamplesPerPixel;

  unsigned int numSamplesPerRay;
  if (vm.count("samples")) numSamplesPerRay = vm["samples"].as<unsigned int>();
  else numSamplesPerRay = ceil( lengthOfRay_in_mm / minSpacing );

  float step_size_in_mm = lengthOfRay_in_mm / numSamplesPerRay;

  size_t numProjs = ps->get_projections()->get_size(2);
  size_t height = ip_h;
  size_t width = ip_w;

  size_t needed_bytes = 2 * prod(imageSize) * sizeof(float);

  std::vector<unsigned int> is_dims = to_std_vector(imageSize);

  std::cout << "IS dimensions " << is_dims[0] << " " << is_dims[1] << " " << is_dims[2] << std::endl;

  std::cout << "Image size " << imageDimensions << std::endl;

  is_dims.push_back(binning->get_number_of_bins());
  std::vector<unsigned int> ps_dims;
  ps_dims.push_back(ps_dims_in_pixels[0]);
  ps_dims.push_back(ps_dims_in_pixels[1]);
  ps_dims.push_back(numProjs);

  boost::shared_ptr< hoCuNDArray<float> > projections (new hoCuNDArray<float>(*ps->get_projections()));

  clamp_min(projections.get(),0.0f);
  //Standard 3d FDK
  // Define encoding matrix
  boost::shared_ptr< hoCudaConebeamProjectionOperator >
    E( new hoCudaConebeamProjectionOperator() );

  /*E->setup( ps_g, ps_bd, ps_g->getAnglesArray(),ps_g->getOffsetXArray(), ps_g->getOffsetYArray(), 1u,
    voxelSize, ps_dims_in_pixels,
    numSamplesPerRay, false);
    E->set_codomain_dimensions(projections->get_dimensions().get()); */
  // Form right hand side

  E->setup(ps,binning,imageDimensions);
  E->set_domain_dimensions(&is_dims);
  E->set_num_projections_per_batch(1);



  //hoGpBbSolver<float> solver;
  hoCuNCGSolver<float> solver;
  solver.set_domain_dimensions(&is_dims);
  solver.set_encoding_operator(E);
  solver.set_output_mode(hoCuNCGSolver<float>::OUTPUT_VERBOSE);
  solver.set_max_iterations(iterations);

  solver.set_non_negativity_constraint(true);

  if (vm.count("TV")){
    std::cout << "Total variation regularization in use" << std::endl;
    boost::shared_ptr<hoCuTvOperator<float,3> > tv(new hoCuTvOperator<float,3>);
    tv->set_weight(vm["TV"].as<float>());
    solver.add_nonlinear_operator(tv);
  }

  if (vm.count("PICS")){
    std::cout << "PICS in used" << std::endl;
    CBCT_binning *binning_pics = new CBCT_binning();
    binning_pics->set_as_default_3d_bin(projections->get_size(2));
    std::vector<unsigned int> is_dims3d = to_std_vector(imageSize);
    boost::shared_ptr< hoCudaConebeamProjectionOperator >
      Ep( new hoCudaConebeamProjectionOperator() );
    /*Ep->setup( ps_g, ps_bd_pics, ps_g->getAnglesArray(), ps_g->getOffsetXArray(), ps_g->getOffsetYArray(), 1u,
      voxelSize, ps_dims_in_pixels,
      numSamplesPerRay, true);*/
    Ep->set_codomain_dimensions(projections->get_dimensions().get());
    // Form right hand side
    Ep->set_domain_dimensions(&is_dims3d);

    boost::shared_ptr<hoCuNDArray<float> > prior3d(new hoCuNDArray<float>(&is_dims3d));
    Ep->mult_MH(projections.get(),prior3d.get());

    hoCuNDArray<float> tmp_proj(*projections);
    Ep->mult_M(prior3d.get(),&tmp_proj);
    float s = dot(projections.get(),&tmp_proj)/dot(&tmp_proj,&tmp_proj);
    *prior3d *= s;
    boost::shared_ptr<hoCuNDArray<float> > prior(new hoCuNDArray<float>(*expand( prior3d.get(), is_dims.back() )));
    boost::shared_ptr<hoCuTvPicsOperator<float,3> > pics (new hoCuTvPicsOperator<float,3>);
    pics->set_prior(prior);
    pics->set_weight(vm["PICS"].as<float>());
    solver.add_nonlinear_operator(pics);
    solver.set_x0(prior);
    delete binning_pics;
  }

  boost::shared_ptr< hoCuNDArray<float> > result = solver.solve(projections.get());
  //boost::shared_ptr< hoCuNDArray<float> > result(new hoCuNDArray<float>(&is_dims));
  //E->mult_MH(projections.get(),result.get());
  write_nd_array<float>( result.get(), outputFile.c_str());


}
