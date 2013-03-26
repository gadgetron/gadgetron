
  template <class T>  void cuNDArray_permute(cuNDArray<T> *in, cuNDArray<T> *out, std::vector<unsigned int> *order, int shift_mode);

// Sum over dimension (scalar and vector_td arrays)
template<class T> EXPORTGPUCORE
boost::shared_ptr<cuNDArray<T> >
sum(cuNDArray<T> *data, unsigned int dim, cuNDA_device alloc_device =
		CUNDA_CURRENT_DEVICE,
		cuNDA_device compute_device = CUNDA_CURRENT_DEVICE);

// Expand (copy) array to new dimension (scalar and vector_td arrays)
template<class T> EXPORTGPUCORE
boost::shared_ptr<cuNDArray<T> >
expand(cuNDArray<T> *data, unsigned int added_dim_size,
		cuNDA_device alloc_device = CUNDA_CURRENT_DEVICE,
		cuNDA_device compute_device = CUNDA_CURRENT_DEVICE);


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

/**
 * Calculates the elementwise maximum of two arrays
 * @param[in] in1 First input array
 * @param[in] in2 Second input Array
 * @param[in] alloc_device Device on which to allocate the new array
 * @param[in] compute_device Device on which to do the computation
 * @return shared pointer to array containing the elementwise maximum of two arrays
 */
template<class T>
boost::shared_ptr< cuNDArray<T> >
maximum( cuNDArray<T> *in1,cuNDArray<T> *in2,
	    cuNDA_device alloc_device, cuNDA_device compute_device );
/**
 * Calculates the elementwise minimum of two arrays
 * @param[in] in1 First input array
 * @param[in] in2 Second input Array
 * @param[in] alloc_device Device on which to allocate the new array
 * @param[in] compute_device Device on which to do the computation
 * @return shared pointer to array containing the elementwise minimum of two arrays
 */
template<class T>
boost::shared_ptr< cuNDArray<T> >
minimum( cuNDArray<T> *in1,cuNDArray<T> *in2,
	    cuNDA_device alloc_device, cuNDA_device compute_device );
// Crop (scalar and vector_td arrays)
template<class T, unsigned int D> EXPORTGPUCORE
void crop(typename uintd<D>::Type crop_offset, cuNDArray<T> *in,
		cuNDArray<T> *out, cuNDA_device compute_device = CUNDA_CURRENT_DEVICE);

// Expand with zero filling (real and complex types)
template<class T, unsigned int D> EXPORTGPUCORE
void expand_with_zero_fill(cuNDArray<T> *in, cuNDArray<T> *out,
		cuNDA_device compute_device = CUNDA_CURRENT_DEVICE);

// Zero fill border (rectangular) - (real and complex types)
template<class T, unsigned int D> EXPORTGPUCORE
void zero_fill_border(typename uintd<D>::Type matrix_size, cuNDArray<T> *image,
		cuNDA_device compute_device = CUNDA_NDARRAY_DEVICE);

// Border fill (circular) - (real and complex types)
template<class REAL, class T, unsigned int D> EXPORTGPUCORE
void zero_fill_border(REAL radius, cuNDArray<T> *image,
		cuNDA_device compute_device = CUNDA_NDARRAY_DEVICE);

// Mirror around the origin -- !! leaving the origin unchanged !!
template<class T, unsigned int D> EXPORTGPUCORE
void origin_mirror(cuNDArray<T> *in, cuNDArray<T> *out, bool zero_fill = true,
		cuNDA_device compute_device = CUNDA_CURRENT_DEVICE);

/**
   * @brief Normalize by the root sum of squares
   * @param[in] x Input array.
   * @return A new complex array containing the input array in the real component and zeros in the imaginary component.
   */

// Normalize by RSS (float/double/complext arrays)
template<class T> EXPORTCPUCOREMATH
void rss_normalize(hoNDArray<T> *in_out, unsigned int dim,
		hoNDA_device compute_device = HONDA_NDARRAY_DEVICE);
