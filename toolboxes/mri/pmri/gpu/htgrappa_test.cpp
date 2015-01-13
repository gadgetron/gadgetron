#include <iostream>
#include <memory>

#include "cuNDArray.h"
#include "hoNDArray_fileio.h"
#include "cuNDFFT.h"
#include "GPUTimer.h"
#include "htgrappa.h"

using namespace Gadgetron;
int main(int argc, char** argv)
{
  GDEBUG_STREAM("Simple HTGRAPPA program" << std::endl);
  {
    GPUTimer init_time("CUDA Initialization");
  }
  GPUTimer total_time("Total time elapsed");
  

  GPUTimer* timer_ptr = new GPUTimer("Loading data");
  hoNDArray<cuFloatComplex> time_average_host = 
    read_nd_array<cuFloatComplex>("time_average.cplx");

  hoNDArray<cuFloatComplex> b1_host = 
    read_nd_array<cuFloatComplex>("b1.cplx");

  cuNDArray<cuFloatComplex> time_average_dev(time_average_host);
  cuNDArray<cuFloatComplex> b1_dev(b1_host);
  delete timer_ptr;

  cuNDArray<cuFloatComplex> unmixing_dev;
  if (!unmixing_dev.create(b1_dev.get_dimensions())) {
    GDEBUG_STREAM("Unable to allocate memory for GRAPPA unmixing coefficients" << std::endl);
    return 0;
  }

  {
    GPUTimer unmix_timer("GRAPPA Unmixing");
    std::vector<unsigned int> kernel_size;
    kernel_size.push_back(5);
    kernel_size.push_back(4);
    if ( htgrappa_calculate_grappa_unmixing(&time_average_dev, 
					    &b1_dev,
					    4,
					    kernel_size,
					    &unmixing_dev) < 0) {
      GDEBUG_STREAM("Error calculating unmixing coefficients" << std::endl);
    }
  }

  /*
  std::auto_ptr< cuNDArray<float2> > b1 = 
    estimate_b1_map<uint2, float, float2>(&time_average_dev);
  */

  timer_ptr = new GPUTimer("Saving data");
  hoNDArray<cuFloatComplex> average_image = time_average_dev.to_host();
  write_nd_array<cuFloatComplex>(average_image, "average_image.cplx");
  delete timer_ptr;

  GDEBUG_STREAM("Reconstruction done" << std::endl);

  return 0;
}
