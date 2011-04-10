#include <iostream>

#include "cuNDArray.h"
#include "hoNDArray_fileio.h"
#include "cuNDFFT.h"
#include "cgOperatorCartesianSense.h"
#include "cuCG.h"

int main(int argc, char** argv)
{
  std::cout << "Simple GPU Test program" << std::endl;
  
  hoNDArray<float2> phantom = read_nd_array<float2>("phantom.cplx");
  hoNDArray<float2> csm = read_nd_array<float2>("csm.cplx");
  hoNDArray<float2> D = read_nd_array<float2>("D.cplx");
  hoNDArray<float>  idxf = read_nd_array<float>("idx.real");

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

  E.set_csm(&csm_dev);
  E.set_sampling_indices(&idx_dev);

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
  cg.set_iterations(20);
  cg.set_limit(1e-10);
  cg.set_output_mode(cuCG<float2>::OUTPUT_VERBOSE);

  cuNDArray<float2> cgresult = cg.solve(&tmp2_out_dev);
  hoNDArray<float2> rho_out = cgresult.to_host();
  write_nd_array<float2>(rho_out,"rho_out.cplx");

  std::cout << "Reconstruction done" << std::endl;

  return 0;
}
