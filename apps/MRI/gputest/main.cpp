#include <iostream>

#include "cuNDArray.h"
#include "hoNDArray_fileio.h"
#include "cuNDFFT.h"
#include "cgOperatorCartesianSense.h"
#include "cgOperatorNonCartesianSense.h"
#include "cuCG.h"
#include "GPUTimer.h"
#include "vector_td.h"

#include <memory>

int main(int argc, char** argv)
{
  std::cout << "Simple GPU Test program" << std::endl;
  
  hoNDArray<float_complext::Type> phantom = read_nd_array<float_complext::Type>("phantom.cplx");
  hoNDArray<float_complext::Type> csm = read_nd_array<float_complext::Type>("csm.cplx");
  hoNDArray<float_complext::Type> D = read_nd_array<float_complext::Type>("D.cplx");
  hoNDArray<float>  idxf = read_nd_array<float>("idx.real");  
  hoNDArray<floatd2::Type>  co = read_nd_array<floatd2::Type>("co.cplx");
  hoNDArray<float>   w = read_nd_array<float>("w.real");
  
  if (csm.get_number_of_dimensions() == 2) {
    std::vector<unsigned int> tmp_reshape_dim = csm.get_dimensions();
    tmp_reshape_dim.push_back(1);
    csm.reshape(tmp_reshape_dim);
  }

  std::cout << "Done reading input data" << std::endl;

  hoNDArray<unsigned int> idx;
  idx.create(idxf.get_dimensions());

  for (unsigned int i = 0; i < idxf.get_number_of_elements(); i++) {
    idx.get_data_ptr()[i] = static_cast<unsigned int>(idxf.get_data_ptr()[i]);
  }

  cuNDArray<float_complext::Type> phantom_dev(phantom);
  cuNDArray<float_complext::Type> *csm_dev = new cuNDArray<float_complext::Type>(csm);
  cuNDArray<unsigned int> *idx_dev = new cuNDArray<unsigned int>(idx);
  cgOperatorCartesianSense<float,2> *E = new cgOperatorCartesianSense<float,2>();
  cuNDArray<floatd2::Type> co_dev(co); co_dev.squeeze();
  cuNDArray<float> *w_dev = new cuNDArray<float>(w);

  E->set_csm(boost::shared_ptr< cuNDArray<float_complext::Type> >(csm_dev));
  E->set_sampling_indices( boost::shared_ptr< cuNDArray<unsigned int> >(idx_dev));

  cgOperatorNonCartesianSense<float,2> *E_noncart = new cgOperatorNonCartesianSense<float,2>();
  uintd2 matrix_size(128,128);
  uintd2 matrix_size_os(256,256);
  float kernel_width = 5.5f;
  E_noncart->setup( matrix_size, matrix_size_os, kernel_width );
  E_noncart->set_csm(boost::shared_ptr<cuNDArray<float_complext::Type> >(csm_dev));
  if (E_noncart->preprocess(&co_dev) < 0) {
    std::cout << "Failed to set trajectory on encoding matrix" << std::endl;
  }
  if (E_noncart->set_dcw(boost::shared_ptr< cuNDArray<float> >(w_dev)) < 0) {
    std::cout << "Failed to set weights on encoding matrix" << std::endl;
  }

  std::vector<unsigned int> dims_out;
  dims_out.push_back(idx.get_number_of_elements());
  dims_out.push_back(csm.get_size(csm.get_number_of_dimensions()-1));

  cuNDArray<float_complext::Type> tmp_out_dev;
  tmp_out_dev.create(dims_out);

  if (E->mult_M(&phantom_dev,&tmp_out_dev,false) < 0) {
    std::cerr << "Failed to multiply with system matrix E" << std::endl;
  }

  hoNDArray<float_complext::Type> tmp_out = tmp_out_dev.to_host();
  write_nd_array<float_complext::Type>(tmp_out,"tmp_out.cplx");

  cuNDArray<float_complext::Type> tmp2_out_dev;
  tmp2_out_dev.create(phantom.get_dimensions());
  
  if (E->mult_MH(&tmp_out_dev,&tmp2_out_dev,false) < 0) {
    std::cerr << "Failed to multiply with system matrix EH" << std::endl;
  }

  hoNDArray<float_complext::Type> tmp2_out = tmp2_out_dev.to_host();
  //write_nd_array<float_complext::Type>(tmp2_out,"tmp2_out.cplx");
  
  cuNDArray<float_complext::Type> *D_dev = new cuNDArray<float_complext::Type>(D);
  
  cuCGPrecondWeight<float_complext::Type> *Dm = new cuCGPrecondWeight<float_complext::Type>();
  Dm->set_weights(boost::shared_ptr< cuNDArray<float_complext::Type> >(D_dev));

  cuCG<float,float_complext::Type> cg;
  cg.add_matrix_operator( boost::shared_ptr< cuCGMatrixOperator<float,float_complext::Type> >(E) );
  cg.set_preconditioner( boost::shared_ptr< cuCGPreconditioner<float_complext::Type> >(Dm) );
  cg.set_iterations(10);
  cg.set_limit(1e-5);
  cg.set_output_mode(cuCG<float, float_complext::Type>::OUTPUT_VERBOSE);

  boost::shared_ptr< cuNDArray<float_complext::Type> > cgresult;
  {
    GPUTimer timer("GPU Conjugate Gradient solve");
    cgresult = cg.solve(&tmp2_out_dev);
  }
  
  hoNDArray<float_complext::Type> rho_out = cgresult->to_host();
  write_nd_array<float_complext::Type>(rho_out,"rho_out.cplx");
  
  std::vector<unsigned int> dims_out_nc;
  dims_out_nc.push_back(co.get_number_of_elements());
  dims_out_nc.push_back(csm.get_size(csm.get_number_of_dimensions()-1));

  cuNDArray<float_complext::Type> tmp_out_nc_dev;
  tmp_out_nc_dev.create(dims_out_nc);

  if (E_noncart->mult_M(&phantom_dev,&tmp_out_nc_dev,false) < 0) {
    std::cerr << "Failed to multiply with system matrix non-cartE" << std::endl;
  }

  hoNDArray<float_complext::Type> tmp_out_nc = tmp_out_nc_dev.to_host();
  write_nd_array<float_complext::Type>(tmp_out_nc,"tmp_out_nc.cplx");

  cuNDArray<float_complext::Type> tmp2_out_nc_dev;
  tmp2_out_nc_dev.create(phantom.get_dimensions());
  

  if (E_noncart->mult_MH(&tmp_out_nc_dev,&tmp2_out_nc_dev,false) < 0) {
    std::cerr << "Failed to multiply with system matrix EH (non cartesian)" << std::endl;
  }
  hoNDArray<float_complext::Type> tmp2_out_nc = tmp2_out_nc_dev.to_host();
  write_nd_array<float_complext::Type>(tmp2_out_nc,"tmp2_out_nc.cplx");

  cuCG<float, float_complext::Type> cg_nc;
  cg_nc.add_matrix_operator( boost::shared_ptr< cuCGMatrixOperator<float,float_complext::Type> >(E_noncart) );
  cg_nc.set_preconditioner(  boost::shared_ptr< cuCGPreconditioner<float_complext::Type> >(Dm) );
  cg_nc.set_iterations(5);
  cg_nc.set_limit(1e-5);
  cg_nc.set_output_mode(cuCG<float, float_complext::Type>::OUTPUT_VERBOSE);

  boost::shared_ptr< cuNDArray<float_complext::Type> > cgresult_nc;
  {
    GPUTimer timer("GPU Conjugate Gradient solve");
    cgresult_nc = cg_nc.solve(&tmp2_out_nc_dev);
  }
  hoNDArray<float_complext::Type> rho_out_nc = cgresult_nc->to_host();
  write_nd_array<float_complext::Type>(rho_out_nc,"rho_out_nc.cplx");


  std::cout << "Reconstruction done" << std::endl;

  return 0;
}
