#include <iostream>

#include "cuNDArray.h"
#include "hoNDArray_fileio.h"
#include "cuNDFFT.h"
#include "cgOperatorCartesianSense.h"
#include "cgOperatorNonCartesianSense.h"
#include "cuCG.h"
#include "GPUTimer.h"

int main(int argc, char** argv)
{
  std::cout << "Simple GPU Test program" << std::endl;
  
  hoNDArray<float2> phantom = read_nd_array<float2>("phantom.cplx");
  hoNDArray<float2> csm = read_nd_array<float2>("csm.cplx");
  hoNDArray<float2> D = read_nd_array<float2>("D.cplx");
  hoNDArray<float>  idxf = read_nd_array<float>("idx.real");  
  hoNDArray<float2>  co = read_nd_array<float2>("co.cplx");
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

  cuNDArray<float2> phantom_dev(phantom);
  cuNDArray<float2> csm_dev(csm);
  cuNDArray<unsigned int> idx_dev(idx);
  cgOperatorCartesianSense E;
  cuNDArray<float2> co_dev(co);
  cuNDArray<float> w_dev(w);

  E.set_csm(&csm_dev);
  E.set_sampling_indices(&idx_dev);


  cgOperatorNonCartesianSense E_noncart;
  E_noncart.set_csm(&csm_dev);
  if (E_noncart.set_trajectories(&co_dev) < 0) {
    std::cout << "Failed to set trajectory on encoding matrix" << std::endl;
  }
  if (E_noncart.set_weights(&w_dev) < 0) {
    std::cout << "Failed to set weights on encoding matrix" << std::endl;
  }

  std::vector<unsigned int> dims_out;
  dims_out.push_back(idx.get_number_of_elements());
  dims_out.push_back(csm.get_size(csm.get_number_of_dimensions()-1));

  cuNDArray<float2> tmp_out_dev;
  tmp_out_dev.create(dims_out);

  if (E.mult_M(&phantom_dev,&tmp_out_dev,false) < 0) {
    std::cerr << "Failed to multiply with system matrix E" << std::endl;
  }

  hoNDArray<float2> tmp_out = tmp_out_dev.to_host();
  write_nd_array<float2>(tmp_out,"tmp_out.cplx");

  cuNDArray<float2> tmp2_out_dev;
  tmp2_out_dev.create(phantom.get_dimensions());
  
  if (E.mult_MH(&tmp_out_dev,&tmp2_out_dev,false) < 0) {
    std::cerr << "Failed to multiply with system matrix EH" << std::endl;
  }

  hoNDArray<float2> tmp2_out = tmp2_out_dev.to_host();
  //write_nd_array<float2>(tmp2_out,"tmp2_out.cplx");
  

  cuNDArray<float2> D_dev(D);
  
  cuCGPrecondWeight<float2> Dm;
  Dm.set_weights(&D_dev);

  cuCG<float2> cg;
  cg.add_matrix_operator(&E, 1.0);
  cg.set_preconditioner(&Dm);
  cg.set_iterations(10);
  cg.set_limit(1e-5);
  cg.set_output_mode(cuCG<float2>::OUTPUT_VERBOSE);

  cuNDArray<float2> cgresult;
  {
    GPUTimer timer("GPU Conjugate Gradient solve");
    cgresult = cg.solve(&tmp2_out_dev);
  }

  hoNDArray<float2> rho_out = cgresult.to_host();
  write_nd_array<float2>(rho_out,"rho_out.cplx");
  
  std::vector<unsigned int> dims_out_nc;
  dims_out_nc.push_back(co.get_number_of_elements());
  dims_out_nc.push_back(csm.get_size(csm.get_number_of_dimensions()-1));

  cuNDArray<float2> tmp_out_nc_dev;
  tmp_out_nc_dev.create(dims_out_nc);

  if (E_noncart.mult_M(&phantom_dev,&tmp_out_nc_dev,false) < 0) {
    std::cerr << "Failed to multiply with system matrix non-cartE" << std::endl;
  }

  hoNDArray<float2> tmp_out_nc = tmp_out_nc_dev.to_host();
  write_nd_array<float2>(tmp_out_nc,"tmp_out_nc.cplx");

  cuNDArray<float2> tmp2_out_nc_dev;
  tmp2_out_nc_dev.create(phantom.get_dimensions());
  

  if (E_noncart.mult_MH(&tmp_out_nc_dev,&tmp2_out_nc_dev,false) < 0) {
    std::cerr << "Failed to multiply with system matrix EH (non cartesian)" << std::endl;
  }
  hoNDArray<float2> tmp2_out_nc = tmp2_out_nc_dev.to_host();
  write_nd_array<float2>(tmp2_out_nc,"tmp2_out_nc.cplx");

  cuCG<float2> cg_nc;
  cg_nc.add_matrix_operator(&E_noncart, 1.0);
  cg_nc.set_preconditioner(&Dm);
  cg_nc.set_iterations(5);
  cg_nc.set_limit(1e-5);
  cg_nc.set_output_mode(cuCG<float2>::OUTPUT_VERBOSE);

  cuNDArray<float2> cgresult_nc;
  {
    GPUTimer timer("GPU Conjugate Gradient solve");
    cgresult_nc = cg_nc.solve(&tmp2_out_nc_dev);
  }
  hoNDArray<float2> rho_out_nc = cgresult_nc.to_host();
  write_nd_array<float2>(rho_out_nc,"rho_out_nc.cplx");


  std::cout << "Reconstruction done" << std::endl;

  return 0;
}
