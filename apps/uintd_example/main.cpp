#include <iostream>

#include "cuNDArray.h"
#include "hoNDArray_fileio.h"
#include "cuNDFFT.h"
#include "uintd.h"
#include "hoNDArray_fileio.h"
#include "cuNDArray_new_permute.h"

int main(int argc, char** argv)
{
  std::cout << "Simple UINTD Test program" << std::endl;

  hoNDArray<float2> phantom = read_nd_array<float2>("phantom.cplx");
  cuNDArray<float2> phantom_dev(phantom);

  std::vector<unsigned int> order;
  order.push_back(1);
  order.push_back(0);

  cuNDArray<float2> out;
  out.create(phantom.get_dimensions());

  cuNDArray_new_permute(&phantom_dev, &out, order);
  
  hoNDArray<float2> phantom_perm = out.to_host();
  write_nd_array<float2>(phantom_perm, "phantom_perm.cplx");

  return 0;
}
