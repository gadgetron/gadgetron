#pragma once

#include "cuNDArray.h"
#include "gpucore_export.h"

namespace Gadgetron{

  template<class T> EXPORTGPUCORE boost::shared_ptr< cuNDArray<T> >
  shift_dim( cuNDArray<T> *in, int shift );

  template<class T> EXPORTGPUCORE void
  shift_dim( cuNDArray<T> *in, cuNDArray<T> *out, int shift );
  
  template<class T> EXPORTGPUCORE boost::shared_ptr< cuNDArray<T> >
  permute( cuNDArray<T> *in, std::vector<unsigned int> *dim_order, int shift_mode = 0 );
  
  template<class T> EXPORTGPUCORE void
  permute( cuNDArray<T> *in, cuNDArray<T> *out, std::vector<unsigned int> *dim_order, int shift_mode = 0 );

  // Expand (copy) array to new dimension (scalar and vector_td arrays)
  template<class T> EXPORTGPUCORE boost::shared_ptr<cuNDArray<T> > expand(cuNDArray<T> *data, unsigned int added_dim_size );
  
  // Sum over dimension (scalar and vector_td arrays)
  template<class T> EXPORTGPUCORE boost::shared_ptr<cuNDArray<T> > sum(cuNDArray<T> *data, unsigned int dim );

  
  /*

  // Correlation matrix over the last dimension in the input array (float/double/complext array)
  template<class T> EXPORTGPUCORE
  boost::shared_ptr<cuNDArray<T> >
  correlation(cuNDArray<T> *data,
  cuNDA_device alloc_device = CUNDA_CURRENT_DEVICE,
  cuNDA_device compute_device = CUNDA_CURRENT_DEVICE);

  // Downsample array to half size (Real arrays only)
  template<class REAL, unsigned int D> EXPORTGPUCORE
  boost::shared_ptr<cuNDArray<REAL> >
  downsample(cuNDArray<REAL> *data, cuNDA_device alloc_device =
  CUNDA_CURRENT_DEVICE,
  cuNDA_device compute_device = CUNDA_CURRENT_DEVICE);

  // Nearest neighbor upsampling of array to double size (Real arrays only)
  template<class REAL, unsigned int D> EXPORTGPUCORE
  boost::shared_ptr<cuNDArray<REAL> >
  upsample_nn(cuNDArray<REAL> *data, cuNDA_device alloc_device =
  CUNDA_CURRENT_DEVICE,
  cuNDA_device compute_device = CUNDA_CURRENT_DEVICE);

  // Linear interpolation upsampling of array to double size (Real arrays only)
  template<class REAL, unsigned int D> EXPORTGPUCORE
  boost::shared_ptr<cuNDArray<REAL> >
  upsample_lin(cuNDArray<REAL> *data, cuNDA_device alloc_device =
  CUNDA_CURRENT_DEVICE,
  cuNDA_device compute_device = CUNDA_CURRENT_DEVICE);
  */
  /**
   * Calculates the elementwise maximum of two arrays
   * @param[in] in1 First input array
   * @param[in] in2 Second input Array
   * @param[in] alloc_device Device on which to allocate the new array
   * @param[in] compute_device Device on which to do the computation
   * @return shared pointer to array containing the elementwise maximum of two arrays
   */
  //template<class T>
  //boost::shared_ptr< cuNDArray<T> >
  //maximum( cuNDArray<T> *in1,cuNDArray<T> *in2,
  //	    cuNDA_device alloc_device, cuNDA_device compute_device );
  /**
   * Calculates the elementwise minimum of two arrays
   * @param[in] in1 First input array
   * @param[in] in2 Second input Array
   * @param[in] alloc_device Device on which to allocate the new array
   * @param[in] compute_device Device on which to do the computation
   * @return shared pointer to array containing the elementwise minimum of two arrays
   */
  /*template<class T>
    boost::shared_ptr< cuNDArray<T> >
    minimum( cuNDArray<T> *in1,cuNDArray<T> *in2,
    cuNDA_device alloc_device, cuNDA_device compute_device );

  // Border fill (circular)
  template<class REAL, class T, unsigned int D> EXPORTGPUNFFT
  void zero_fill_border(REAL radius, cuNDArray<T> *image );


    // Mirror around the origin -- !! leaving the origin unchanged !!
    template<class T, unsigned int D> EXPORTGPUCORE
    void origin_mirror(cuNDArray<T> *in, cuNDArray<T> *out, bool zero_fill = true,
    cuNDA_device compute_device = CUNDA_CURRENT_DEVICE);
  */

  // Normalize by RSS (float/double/complext arrays)
}
