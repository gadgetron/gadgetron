/*
 * Copyright (C) 2011 Florian Rathgeber, florian.rathgeber@gmail.com
 *
 * This code is licensed under the MIT License.  See the FindCUDA.cmake script
 * for the text of the license.
 *
 * Based on code by Christopher Bruns published on Stack Overflow (CC-BY):
 * http://stackoverflow.com/questions/2285185
 */

#include <stdio.h>
#include <cuda_runtime.h>

int main() {
  int deviceCount, device, major = 9999, minor = 9999;
  int gpuDeviceCount = 0;
  struct cudaDeviceProp properties;

  if (cudaGetDeviceCount(&deviceCount) != cudaSuccess)
    return 1;
  /* machines with no GPUs can still report one emulation device */
  for (device = 0; device < deviceCount; ++device) {
    cudaGetDeviceProperties(&properties, device);
    if (properties.major != 9999 && properties.major > 1) {/* 9999 means emulation only and we do not support compute model 1.x*/
      ++gpuDeviceCount;
      if (gpuDeviceCount > 1)
      	printf(";");
      if (properties.major == 2) //Need a special case for Fermi. Compute capability 2.1 exists, but compute model 2.1 does not.
      	printf("%d%d",properties.major, 0);
      else
      	printf("%d%d",properties.major, properties.minor);
      /*  get minimum compute capability of all devices */
    }
  }
  /* don't just return the number of gpus, because other runtime cuda
     errors can also yield non-zero return values */
  if (gpuDeviceCount > 0) {
    /* this output will be parsed by FindCUDA.cmake */
    return 0; /* success */
  }
  return 1; /* failure */
}
