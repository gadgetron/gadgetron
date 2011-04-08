#include <iostream>

#include "cuNDArray.h"
#include "hoNDArray_fileio.h"
#include "cuNDFFT.h"
#include "cgOperatorCartesianSense.h"

int main(int argc, char** argv)
{
  std::cout << "Simple GPU Test program" << std::endl;

  
  hoNDArray<float2> phantom = read_nd_array<float2>("phantom.cplx");
  hoNDArray<float2> csm = read_nd_array<float2>("csm.cplx");
  hoNDArray<float>  idxf = read_nd_array<float>("idx.real");

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
  write_nd_array<float2>(tmp2_out,"tmp2_out.cplx");


  if (E.mult_MH_M(&phantom_dev,&tmp2_out_dev,false) < 0) {
    std::cerr << "Failed to multiply with system matrix EHE" << std::endl;
  }

  hoNDArray<float2> tmp3_out = tmp2_out_dev.to_host();
  write_nd_array<float2>(tmp3_out,"tmp3_out.cplx");


  /*
  std::vector<unsigned int> dims;
  dims.push_back(128);
  dims.push_back(1);
  dims.push_back(128);

  hoNDArray<float2> hoND_tmp4;
  hoND_tmp4.create(dims);
  
  std::cout << "Dimensions before squeeze: ";
  for (unsigned int i = 0; i < hoND_tmp4.get_number_of_dimensions(); i++) {
    std::cout << hoND_tmp4.get_size(i) << " ";
  }
  std::cout << std::endl;
  hoND_tmp4.squeeze();

  std::cout << "Dimensions before squeeze: ";
  for (unsigned int i = 0; i < hoND_tmp4.get_number_of_dimensions(); i++) {
    std::cout << hoND_tmp4.get_size(i) << " ";
  }
  std::cout << std::endl;


  hoNDArray<float2> hoND_tmp = read_nd_array<float2>("test_data.cplx");

  
  hoNDArray<float2> hoND_tmp3;
  hoND_tmp3.create(dims);

  hoND_tmp3 = hoND_tmp;
  write_nd_array<float2>(hoND_tmp3,"test_data3.cplx");

  cuNDArray<float2> cuND_tmp(hoND_tmp);
  
  cuNDFFT ft;

  std::vector<unsigned int> ft_dims;
  ft_dims.push_back(1);
  ft_dims.push_back(2);

  ft.fft(&cuND_tmp);
  ft.ifft(&cuND_tmp);
  

  //cuND_tmp.permute(ft_dims,0,1);
  //cuND_tmp.permute(ft_dims,0,-1);
  hoNDArray<float2> hoND_tmp2 = cuND_tmp.to_host();

  write_nd_array<float2>(hoND_tmp2,"test_data_fft.cplx");

  */

  /*
  hoNDArray<float> hoND_tmp = read_nd_array<float>("test_data.real");
  
  cuNDArray<float> cuND_tmp(hoND_tmp);



  std::vector<unsigned int> permute_order;
  permute_order.push_back(1);
  permute_order.push_back(0);

  if (cuND_tmp.permute(permute_order, 0 ,-1) < 0) {
    std::cerr << "Unable to permute" << std::endl;
  }

  cuNDArray<float> cuND_tmp2;
  if (!cuND_tmp2.create(cuND_tmp.get_dimensions(), cuND_tmp.get_data_ptr())) {
    std::cerr << "Unable to create new array" << std::endl;
  }

  cuNDArray<float>* cuND_tmp22 = new cuNDArray<float>(cuND_tmp2);


  hoNDArray<float> hoND_tmp2 = cuND_tmp.to_host();

  write_nd_array<float>(hoND_tmp2,"test_data_perm.real");

  cuNDArray<float>* cuND_tmp3 = 0;
  cuND_tmp3 = cuNDArray<float>::allocate(cuND_tmp.get_dimensions());
  std::cout << "CUDA Array allocate with address: " << cuND_tmp3 << std::endl;
  delete cuND_tmp3;

  hoNDArray<float>* hoND_tmp3 = 0;
  hoND_tmp3 = hoNDArray<float>::allocate(hoND_tmp.get_dimensions());
  std::cout << "HOST Array allocate with address: " << hoND_tmp3 << std::endl;
  delete hoND_tmp3;
  */

  /*
  if (!cuND_tmp.create(dims)) {
    std::cout << "Failed to allocate array on GPU" << std::endl;
  }
  */

  return 0;
}
