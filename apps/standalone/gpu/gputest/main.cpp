#include <iostream>

#include "cuNDArray.h"
#include "hoNDArray_fileio.h"
#include "cuNDFFT.h"
#include "cuCartesianSenseOperator.h"
#include "cuNonCartesianSenseOperator.h"
#include "cuCGPrecondWeights.h"
#include "cuCGSolver.h"
#include "GPUTimer.h"
#include "vector_td.h"

#include <memory>

int main(int argc, char** argv)
{
  std::cout << "Simple GPU Test program" << std::endl;
  
  boost::shared_ptr< hoNDArray<float_complext> > phantom = read_nd_array<float_complext>("phantom.cplx");
  boost::shared_ptr< hoNDArray<float_complext> > csm = read_nd_array<float_complext>("csm.cplx");
  boost::shared_ptr< hoNDArray<float_complext> > D = read_nd_array<float_complext>("D.cplx");
  boost::shared_ptr< hoNDArray<float> > idxf = read_nd_array<float>("idx.real");  
  boost::shared_ptr< hoNDArray<floatd2::Type> > co = read_nd_array<floatd2::Type>("co.cplx");
  boost::shared_ptr< hoNDArray<float> > w = read_nd_array<float>("w.real");
  
  if (csm->get_number_of_dimensions() == 2) {
    std::vector<unsigned int> tmp_reshape_dim = *csm->get_dimensions();
    tmp_reshape_dim.push_back(1);
    csm->reshape(&tmp_reshape_dim);
  }

  std::cout << "Done reading input data" << std::endl;

  hoNDArray<unsigned int> idx;
  idx.create(idxf->get_dimensions().get());

  for (unsigned int i = 0; i < idxf->get_number_of_elements(); i++) {
    idx.get_data_ptr()[i] = static_cast<unsigned int>(idxf->get_data_ptr()[i]);
  }

  cuNDArray<float_complext> phantom_dev(phantom.get());
  cuNDArray<float_complext> *csm_dev = new cuNDArray<float_complext>(csm.get());
  cuNDArray<unsigned int> *idx_dev = new cuNDArray<unsigned int>(&idx);
  cuCartesianSenseOperator<float,2> *E = new cuCartesianSenseOperator<float,2>();
  cuNDArray<floatd2::Type> co_dev(co.get()); co_dev.squeeze();
  cuNDArray<float> *w_dev = new cuNDArray<float>(w.get());

  E->set_csm(boost::shared_ptr< cuNDArray<float_complext> >(csm_dev));
  E->set_sampling_indices( boost::shared_ptr< cuNDArray<unsigned int> >(idx_dev));

  cuNonCartesianSenseOperator<float,2> *E_noncart = new cuNonCartesianSenseOperator<float,2>();
  uintd2 matrix_size(128,128);
  uintd2 matrix_size_os(256,256);
  float kernel_width = 5.5f;
  E_noncart->setup( matrix_size, matrix_size_os, kernel_width );
  E_noncart->set_csm(boost::shared_ptr<cuNDArray<float_complext> >(csm_dev));
  if (E_noncart->preprocess(&co_dev) < 0) {
    std::cout << "Failed to set trajectory on encoding matrix" << std::endl;
  }
  if (E_noncart->set_dcw(boost::shared_ptr< cuNDArray<float> >(w_dev)) < 0) {
    std::cout << "Failed to set weights on encoding matrix" << std::endl;
  }

  std::vector<unsigned int> dims_out;
  dims_out.push_back(idx.get_number_of_elements());
  dims_out.push_back(csm->get_size(csm->get_number_of_dimensions()-1));

  cuNDArray<float_complext> tmp_out_dev;
  tmp_out_dev.create(&dims_out);

  if (E->mult_M(&phantom_dev,&tmp_out_dev,false) < 0) {
    std::cerr << "Failed to multiply with system matrix E" << std::endl;
  }

  boost::shared_ptr< hoNDArray<float_complext> > tmp_out = tmp_out_dev.to_host();
  write_nd_array<float_complext>(tmp_out.get(),"tmp_out.cplx");

  cuNDArray<float_complext> tmp2_out_dev;
  tmp2_out_dev.create(phantom->get_dimensions().get());
  
  if (E->mult_MH(&tmp_out_dev,&tmp2_out_dev,false) < 0) {
    std::cerr << "Failed to multiply with system matrix EH" << std::endl;
  }

  boost::shared_ptr< hoNDArray<float_complext> > tmp2_out = tmp2_out_dev.to_host();
  //write_nd_array<float_complext>(tmp2_out,"tmp2_out.cplx");
  
  cuNDArray<float_complext> *D_dev = new cuNDArray<float_complext>(D.get());
  
  cuCGPrecondWeights<float_complext> *Dm = new cuCGPrecondWeights<float_complext>();
  Dm->set_weights(boost::shared_ptr< cuNDArray<float_complext> >(D_dev));

  cuCGSolver<float,float_complext> cg;
  cg.add_matrix_operator( boost::shared_ptr< matrixOperator<float, cuNDArray<float_complext> > >(E) );
  cg.set_preconditioner( boost::shared_ptr< cuCGPreconditioner<float_complext> >(Dm) );
  cg.set_iterations(10);
  cg.set_limit(1e-5);
  cg.set_output_mode(cuCGSolver<float, float_complext>::OUTPUT_VERBOSE);

  boost::shared_ptr< cuNDArray<float_complext> > cgresult;
  {
    GPUTimer timer("GPU Conjugate Gradient solve");
    cgresult = cg.solve(&tmp2_out_dev);
  }
  
  boost::shared_ptr< hoNDArray<float_complext> > rho_out = cgresult->to_host();
  write_nd_array<float_complext>(rho_out.get(),"rho_out.cplx");
  
  std::vector<unsigned int> dims_out_nc;
  dims_out_nc.push_back(co->get_number_of_elements());
  dims_out_nc.push_back(csm->get_size(csm->get_number_of_dimensions()-1));

  cuNDArray<float_complext> tmp_out_nc_dev;
  tmp_out_nc_dev.create(&dims_out_nc);

  if (E_noncart->mult_M(&phantom_dev,&tmp_out_nc_dev,false) < 0) {
    std::cerr << "Failed to multiply with system matrix non-cartE" << std::endl;
  }

  boost::shared_ptr< hoNDArray<float_complext> > tmp_out_nc = tmp_out_nc_dev.to_host();
  write_nd_array<float_complext>(tmp_out_nc.get(),"tmp_out_nc.cplx");

  cuNDArray<float_complext> tmp2_out_nc_dev;
  tmp2_out_nc_dev.create(phantom->get_dimensions().get());
  

  if (E_noncart->mult_MH(&tmp_out_nc_dev,&tmp2_out_nc_dev,false) < 0) {
    std::cerr << "Failed to multiply with system matrix EH (non cartesian)" << std::endl;
  }
  
  boost::shared_ptr< hoNDArray<float_complext> > tmp2_out_nc = tmp2_out_nc_dev.to_host();
  write_nd_array<float_complext>(tmp2_out_nc.get(),"tmp2_out_nc.cplx");

  cuCGSolver<float, float_complext> cg_nc;
  cg_nc.add_matrix_operator( boost::shared_ptr< matrixOperator< float, cuNDArray<float_complext> > >(E_noncart) );
  cg_nc.set_preconditioner(  boost::shared_ptr< cuCGPreconditioner<float_complext> >(Dm) );
  cg_nc.set_iterations(5);
  cg_nc.set_limit(1e-5);
  cg_nc.set_output_mode(cuCGSolver<float, float_complext>::OUTPUT_VERBOSE);

  boost::shared_ptr< cuNDArray<float_complext> > cgresult_nc;
  {
    GPUTimer timer("GPU Conjugate Gradient solve");
    cgresult_nc = cg_nc.solve(&tmp2_out_nc_dev);
  }
  
  boost::shared_ptr< hoNDArray<float_complext> > rho_out_nc = cgresult_nc->to_host();
  write_nd_array<float_complext>(rho_out_nc.get(),"rho_out_nc.cplx");


  std::cout << "Reconstruction done" << std::endl;

  return 0;
}
