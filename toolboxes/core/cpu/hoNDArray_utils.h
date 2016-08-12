#pragma once

#include "hoNDArray.h"
#include "vector_td_utilities.h"

#ifdef USE_OMP
#include <omp.h>
#endif

namespace Gadgetron {

  class ArrayIterator
  {
  public:

    ArrayIterator(std::vector<size_t> *dimensions, std::vector<size_t> *order)
    {
      dimensions_  = boost::shared_ptr< std::vector<size_t> > (new std::vector<size_t>);
      order_       = boost::shared_ptr< std::vector<size_t> > (new std::vector<size_t>);
      current_     = boost::shared_ptr< std::vector<size_t> > (new std::vector<size_t>);
      block_sizes_ = boost::shared_ptr< std::vector<size_t> > (new std::vector<size_t>);

      block_sizes_->push_back(1);
      for (size_t i = 0; i < order->size(); i++) {
        dimensions_->push_back((*dimensions)[i]);
        order_->push_back((*order)[i]);
        current_->push_back(0);
        if (i > 0) {
          block_sizes_->push_back((*block_sizes_)[i-1]*(*dimensions_)[i-1]);
        }
      }
      current_idx_ = 0;
    }

    inline size_t advance()
    {
      size_t order_index = 0;
      (*current_)[(*order_)[order_index]]++;
      while ((*current_)[(*order_)[order_index]] >= (*dimensions_)[(*order_)[order_index]]) {
        (*current_)[(*order_)[order_index]] = 0;
        order_index = (order_index+1)%dimensions_->size();
        (*current_)[(*order_)[order_index]]++;
      }

      current_idx_ = 0;
      for (size_t i = 0; i < dimensions_->size(); i++) {
        current_idx_ += (*current_)[i]*(*block_sizes_)[i];
      }	
      return current_idx_;
    }

    inline size_t get_current_idx() {
      return current_idx_;
    }

    boost::shared_ptr< std::vector<size_t> > get_current_sub() {
      return current_;
    }

  protected:
    boost::shared_ptr< std::vector<size_t> > dimensions_;
    boost::shared_ptr< std::vector<size_t> > order_;
    boost::shared_ptr< std::vector<size_t> > current_;
    boost::shared_ptr< std::vector<size_t> > block_sizes_;
    size_t current_idx_;
  };

  template<class T> boost::shared_ptr< hoNDArray<T> > shift_dim( hoNDArray<T> *in, int shift )  
  {
    if( in == 0x0 ) {
      throw std::runtime_error("shift_dim(): invalid input pointer provided");;
    }    
    std::vector<size_t> order;
    for (size_t i = 0; i < in->get_number_of_dimensions(); i++) {
      order.push_back(static_cast<size_t>((i+shift)%in->get_number_of_dimensions()));
    }
    return permute(in,&order);
  }

  template<class T> void shift_dim( hoNDArray<T> *in, hoNDArray<T> *out, int shift )
  {
    if( in == 0x0 || out == 0x0 ) {
      throw std::runtime_error("shift_dim(): invalid pointer provided");;
    }    
    std::vector<size_t> order;
    for (size_t i = 0; i < in->get_number_of_dimensions(); i++) {
      order.push_back(static_cast<size_t>((i+shift)%in->get_number_of_dimensions()));
    }
    permute(in,out,&order);
  }

  template<class T> boost::shared_ptr< hoNDArray<T> > 
  permute( hoNDArray<T> *in, std::vector<size_t> *dim_order, int shift_mode = 0) 
  {
    if( in == 0x0 || dim_order == 0x0 ) {
      throw std::runtime_error("permute(): invalid pointer provided");;
    }    

    std::vector<size_t> dims;
    for (size_t i = 0; i < dim_order->size(); i++)
      dims.push_back(in->get_dimensions()->at(dim_order->at(i)));
    boost::shared_ptr< hoNDArray<T> > out( new hoNDArray<T>() );    
    out->create(&dims);
    permute( in, out.get(), dim_order, shift_mode );
    return out;
  }

  template<class T> void 
  permute( hoNDArray<T> *in, hoNDArray<T> *out, std::vector<size_t> *dim_order, int shift_mode = 0) 
  {
    if( in == 0x0 || out == 0x0 || dim_order == 0x0 ) {
      throw std::runtime_error("permute(): invalid pointer provided");;
    }    

    if( in == out ){
      throw std::runtime_error("permute(): in-place permutation not supported");;
    }   

    // Check ordering array
    if (dim_order->size() > in->get_number_of_dimensions()) {
      throw std::runtime_error("hoNDArray::permute - Invalid length of dimension ordering array");;
    }

    std::vector<size_t> dim_count(in->get_number_of_dimensions(),0);
    for (size_t i = 0; i < dim_order->size(); i++) {
      if ((*dim_order)[i] >= in->get_number_of_dimensions()) {
        throw std::runtime_error("hoNDArray::permute - Invalid dimension order array");;
      }
      dim_count[(*dim_order)[i]]++;
    }

    // Create an internal array to store the dimensions
    std::vector<size_t> dim_order_int;

    // Check that there are no duplicate dimensions
    for (size_t i = 0; i < dim_order->size(); i++) {
      if (dim_count[(*dim_order)[i]] != 1) {
        throw std::runtime_error("hoNDArray::permute - Invalid dimension order array (duplicates)");;

      }
      dim_order_int.push_back((*dim_order)[i]);
    }

    for (size_t i = 0; i < dim_order_int.size(); i++) {
      if ((*in->get_dimensions())[dim_order_int[i]] != out->get_size(i)) {
        throw std::runtime_error("permute(): dimensions of output array do not match the input array");;
      }
    }

    // Pad dimension order array with dimension not mentioned in order array
    if (dim_order_int.size() < in->get_number_of_dimensions()) {
      for (size_t i = 0; i < dim_count.size(); i++) {
        if (dim_count[i] == 0) {
          dim_order_int.push_back(i);
        }
      }
    }

    T* o = out->get_data_ptr();

    // if memcpy can be used during permute
    size_t stride = 1;
    size_t num_dim_memcpy = 0;
    for (size_t i = 0; i < dim_order_int.size(); i++) {
        if (dim_order_int[i]==i){
            stride *= in->get_size(i);
            num_dim_memcpy = i;
        }
        else{
            break;
        }
    }

    if (stride == 1) {
        // point by point assignment is needed
        ArrayIterator it(in->get_dimensions().get(), &dim_order_int);
        for (size_t i = 0; i < in->get_number_of_elements(); i++) {
            o[i] = in->get_data_ptr()[it.get_current_idx()];
            it.advance();
        }
    }
    else {
        // memcpy can be used

        size_t nDim = in->get_number_of_dimensions();
        size_t num_memcpy = in->get_number_of_elements() / stride;

        if (num_dim_memcpy == nDim - 1){
            memcpy(out->begin(), in->begin(), in->get_number_of_bytes());
            return;
        }

        // for the array index calculation
        std::vector<size_t> dim_permute(nDim-num_dim_memcpy-1);
        for (size_t i = num_dim_memcpy+1; i < dim_order_int.size(); i++) {
            dim_permute[i - num_dim_memcpy - 1] = in->get_size(i);
        }

        size_t n;

        hoNDArray<T> permuteArray(dim_permute, in->begin(), false);

        // starting index for in and out array for every permute memcpy operation
        std::vector<size_t> ind_permute_in(dim_permute.size(), 0), ind_in(nDim, 0), ind_out(nDim, 0);

        for (n = 0; n < num_memcpy; n++) {
            permuteArray.calculate_index(n, ind_permute_in);
            memcpy(&ind_in[0] + num_dim_memcpy + 1, &ind_permute_in[0], sizeof(size_t)*ind_permute_in.size());

            // permute the indexes
            for (size_t i = 0; i < nDim; i++) {
                ind_out[i] = ind_in[dim_order_int[i]];
            }

            size_t offset_in = in->calculate_offset(ind_in);
            size_t offset_out = out->calculate_offset(ind_out);

            memcpy(o + offset_out, in->begin() + offset_in, sizeof(T)*stride);
        }
    }
  }

  // Expand array to new dimension
  template<class T> boost::shared_ptr<hoNDArray<T> > 
  expand(hoNDArray<T> *in, size_t new_dim_size )
  {
    if( in == 0x0 ){
      throw std::runtime_error("expand(): illegal input pointer.");;
    }

    const size_t number_of_elements_in = in->get_number_of_elements();    

    std::vector<size_t> dims = *in->get_dimensions(); 
    dims.push_back(new_dim_size);

    boost::shared_ptr< hoNDArray<T> > out(new hoNDArray<T>(&dims));

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for( long long int idx=0; idx<number_of_elements_in*new_dim_size; idx++ ){
      (*out)[idx] = in->at(idx%number_of_elements_in);
    }
    return out;
  }
  
  // Sum over dimension
  template<class T> boost::shared_ptr<hoNDArray<T> > 
  sum(hoNDArray<T> *in, size_t dim )
  {
    if( in == 0x0 ){
      throw std::runtime_error("sum(): illegal input pointer.");;
    }

    if( !(in->get_number_of_dimensions()>1) ){
      throw std::runtime_error("sum(): underdimensioned.");;
    }

    if( dim > in->get_number_of_dimensions()-1 ){
      throw std::runtime_error( "sum(): dimension out of range.");;
    }

    size_t number_of_batches = in->get_size(dim);
    size_t number_of_elements = in->get_number_of_elements()/number_of_batches;
    std::vector<size_t> dims = *in->get_dimensions(); dims.pop_back();

    boost::shared_ptr< hoNDArray<T> > out(new hoNDArray<T>());
    out->create(&dims);

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for( long long idx=0; idx<(long long)number_of_elements; idx++ ){
      T val(0);
      for( size_t j=0; j<number_of_batches; j++ ){
        size_t in_idx = j*number_of_elements+idx;
        val += in->get_data_ptr()[in_idx];      
      }
      out->get_data_ptr()[idx] = val;
    }
    return out;
  } 

  /**
  * @param[in] crop_offset starting position to crop
  * @param[in] crop_size Size of cropped array
  * @param[in] in input array
  * @param[out] out Output array after cropping
  */
  template<class T, unsigned int D> void
  crop(const vector_td<size_t, D>& crop_offset, const vector_td<size_t, D>& crop_size, hoNDArray<T> *in, hoNDArray<T> *out)
  {
      if (in == 0x0){
          throw std::runtime_error("crop: 0x0 array provided");;
      }

      if (in->get_number_of_dimensions() < D){
          std::stringstream ss;
          ss << "crop: number of image dimensions should be at least " << D;
          throw std::runtime_error(ss.str());;
      }

      if (out == 0x0){
          throw std::runtime_error("crop: 0x0 array provided");;
      }

      std::vector<size_t> dims = to_std_vector(crop_size);
      for (unsigned int d = D; d<in->get_number_of_dimensions(); d++){
          dims.push_back(in->get_size(d));
      }

      if (!out->dimensions_equal(&dims)){
          out->create(dims);
      }

      typename uint64d<D>::Type matrix_size_in = from_std_vector<size_t, D>(*in->get_dimensions());
      typename uint64d<D>::Type matrix_size_out = from_std_vector<size_t, D>(*out->get_dimensions());

      if (weak_greater(crop_offset + matrix_size_out, matrix_size_in)){
          throw std::runtime_error("crop: cropping size mismatch");;
      }

      size_t len = out->get_size(0);
      size_t num = out->get_number_of_elements() / len;

      long long k;

      T *in_ptr = in->get_data_ptr();
      T *out_ptr = out->get_data_ptr();

      #pragma omp parallel default(none) private(k) shared(in_ptr, out_ptr, num, len, in, out, crop_offset)
      {
          std::vector<size_t> ind;

      #pragma omp for 
          for (k = 0; k < (long long)num; k++){
              ind = out->calculate_index(k*len);
              for (unsigned int d = 0; d < D; d++){
                  ind[d] += crop_offset[d];
              }

              T* in_ptr_curr = in_ptr + in->calculate_offset(ind);
              memcpy(out_ptr + k*len, in_ptr_curr, sizeof(T)*len);
          }
      }
  }

  /**
  * @param[in] crop_size Size of cropped array
  * @param[in] in input array
  * @param[out] out Output array after cropping

  * Crop the input array around its center N/2; that is, the center pixel of in array is the center pixel of out array
  */
  template<class T, unsigned int D> void
  crop(const vector_td<size_t, D>& crop_size, hoNDArray<T> *in, hoNDArray<T> *out)
  {
    // compute crop offset, perserving the center
    if (in == 0x0){
        throw std::runtime_error("crop: 0x0 array provided");;
    }

    vector_td<size_t, D> crop_offset;

    unsigned int d;
    for (d = 0; d < D; d++)
    {
        crop_offset[d] = in->get_size(d) / 2 - crop_size[d] / 2;
    }

    crop(crop_offset, crop_size, in, out);
  }

  template<class T> void
  crop(size_t x, hoNDArray<T> *in, hoNDArray<T> *out)
  {
      vector_td<size_t, 1> crop_size(x);
      crop(crop_size, in, out);
  }

  template<class T> void
  crop(size_t x, size_t y, hoNDArray<T> *in, hoNDArray<T> *out)
  {
      vector_td<size_t, 2> crop_size(x, y);
      crop(crop_size, in, out);
  }

  template<class T> void
  crop(size_t x, size_t y, size_t z, hoNDArray<T> *in, hoNDArray<T> *out)
  {
      vector_td<size_t, 3> crop_size(x, y, z);
      crop(crop_size, in, out);
  }

  template<class T, unsigned int D> boost::shared_ptr< hoNDArray<T> >
  crop( const vector_td<size_t, D>& crop_offset, const vector_td<size_t, D>& crop_size, hoNDArray<T> *in )
  {
    boost::shared_ptr< hoNDArray<T> > out( new hoNDArray<T>() );
    crop(crop_offset, crop_size, in, out.get());
    return out;
  }

  /**
   * @param[in]     size    Size of the output array
   * @param[in]     in      Input array
   * @param[out]    out     Output array after padding
   * @param[in]     preset_out_with_val if true, out array will be filled with val before padding
   * @param[in]     val     Value to use for padding

   * The padding operations keep the center of array unchanged, e.g. the center is always N/2
   */
  template<class T, unsigned int D> void
  pad(const typename uint64d<D>::Type& size, hoNDArray<T> *in, hoNDArray<T>* out, bool preset_out_with_val = true, T val = T(0))
  {
      if (in == 0x0){
          throw std::runtime_error("pad: 0x0 array provided");;
      }

      if (out == 0x0){
          throw std::runtime_error("pad: 0x0 array provided");;
      }

      if (in->get_number_of_dimensions() < D){
          std::stringstream ss;
          ss << "pad: number of image dimensions should be at least " << D;
          throw std::runtime_error(ss.str());;
      }

      unsigned int d;

      std::vector<size_t> dims = to_std_vector(size);
      for (d = D; d<in->get_number_of_dimensions(); d++){
          dims.push_back(in->get_size(d));
      }

      if (!out->dimensions_equal(&dims)){
          out->create(dims);
      }

      if (in->dimensions_equal(&dims)){
          memcpy(out->begin(), in->begin(), in->get_number_of_bytes());
          return;
      }

      T *in_ptr = in->get_data_ptr();
      T *out_ptr = out->get_data_ptr();

      if (preset_out_with_val){
          if (val == T(0)){
              memset(out_ptr, 0, out->get_number_of_bytes());
          }
          else{
                size_t N = out->get_number_of_elements();
                long long n;
                #pragma omp parallel for default(none) private(n) shared(N, out_ptr, val)
                for (n = 0; n<(long long)N; n++)
                {
                    out_ptr[n] = val;
                }
          }
      }

      typename uint64d<D>::Type matrix_size_in = from_std_vector<size_t, D>(*in->get_dimensions());
      typename uint64d<D>::Type matrix_size_out = from_std_vector<size_t, D>(*out->get_dimensions());

      if (weak_greater(matrix_size_in, matrix_size_out)){
          throw std::runtime_error("pad: size mismatch, cannot expand");
      }

      typename uint64d<D>::Type offset(D);
      for (d = 0; d<D; d++){
          offset[d] = matrix_size_out[d]/2 - matrix_size_in[d]/2;
      }

      size_t len = in->get_size(0);
      size_t num = in->get_number_of_elements() / len;

      long long k;

#pragma omp parallel default(none) private(k, d) shared(in_ptr, out_ptr, num, len, in, out, offset)
      {
          std::vector<size_t> ind;

#pragma omp for 
          for (k = 0; k < (long long)num; k++){
              ind = in->calculate_index(k*len);
              for (d = 0; d < D; d++){
                  ind[d] += offset[d];
              }

              T* out_ptr_curr = out_ptr + out->calculate_offset(ind);
              memcpy(out_ptr_curr, in_ptr + k*len, sizeof(T)*len);
          }
      }
  }

  template<class T> void
  pad(size_t x, hoNDArray<T> *in, hoNDArray<T>* out, bool preset_out_with_val = true, T val = T(0))
  {
      typename uint64d<1>::Type padSize(x);
      pad<T, 1>(padSize, in, out, preset_out_with_val, val);
  }

  template<class T> void
  pad(size_t x, size_t y, hoNDArray<T> *in, hoNDArray<T>* out, bool preset_out_with_val = true, T val = T(0))
  {
      typename uint64d<2>::Type padSize(x, y);
      pad<T, 2>(padSize, in, out, preset_out_with_val, val);
  }

  template<class T> void
  pad(size_t x, size_t y, size_t z, hoNDArray<T> *in, hoNDArray<T>* out, bool preset_out_with_val = true, T val = T(0))
  {
      typename uint64d<3>::Type padSize(x, y, z);
      pad<T, 3>(padSize, in, out, preset_out_with_val, val);
  }

  /**
  * @param[in] size Size of the output array
  * @param[in] in Input array
  * @param[in] val Value to use for padding
  * @returns New array of the specified size, containing the original input array in the center and val outside.
  */
  template<class T, unsigned int D> boost::shared_ptr< hoNDArray<T> >
  pad(const typename uint64d<D>::Type& size, hoNDArray<T> *in, T val = T(0))
  {
    boost::shared_ptr< hoNDArray<T> > out(new hoNDArray<T>());
    pad<T,D>(size, in, out.get(), true, val);
    return out;
  }

  /// copy the sub array x(:, indLastDim) to all other places of the last dimensions
  template<typename T> 
  bool repmatLastDimension(hoNDArray<T>& x, size_t indLastDim)
  {
    try
      {
        size_t NDim = x.get_number_of_dimensions();
        size_t lastDim = x.get_size(NDim-1);
        GADGET_CHECK_RETURN_FALSE( indLastDim < lastDim );

        std::vector<size_t> ind(NDim, 0);
        ind[NDim-1] = indLastDim;
        size_t offsetIndLastDim = x.calculate_offset(ind);

        size_t N = x.get_number_of_elements() / lastDim;

        long long l;
#pragma omp parallel default(none) private(l) shared(lastDim, offsetIndLastDim, x, ind, indLastDim, N, NDim)
        {
            std::vector<size_t> indLocal(ind);

#pragma omp for
            for ( l=0; l<(long long)lastDim; l++ )
            {
                if ( l==indLastDim ) continue;
                indLocal[NDim-1] = l;
                size_t offsetInd = x.calculate_offset(indLocal);

                memcpy(x.begin()+offsetInd, x.begin()+offsetIndLastDim, sizeof(T)*N);
            }
        }
      }
    catch (...)
      {
        GERROR_STREAM("Errors in repmatLastDimension(hoNDArray<T>& x, size_t indLastDim) ... ");
        return false;
      }
    return true;
  }

  // Utility to check if all neighbors required for the linear interpolation exists
  // ... do not include dimensions of size 1

  template<class REAL, unsigned int D> inline bool
  is_border_pixel( vector_td<size_t,D> co, vector_td<size_t,D> dims )
  {
    for( size_t dim=0; dim<D; dim++ ){
      if( dims[dim] > 1 && ( co[dim] == 0 || co[dim] == (dims[dim]-1) ) )
	return true;
    }
    return false;
  }

  // Downsample
  template<class REAL, unsigned int D> 
  boost::shared_ptr< hoNDArray<REAL> > downsample( hoNDArray<REAL> *_in )
  {
    // A few sanity checks 

    if( _in == 0x0 ){
      throw std::runtime_error( "downsample(): illegal input provided.");
    }
    
    if( _in->get_number_of_dimensions() < D ){
      throw std::runtime_error( "downsample(): the number of array dimensions should be at least D");
    }
    
    for( size_t d=0; d<D; d++ ){
      if( (_in->get_size(d)%2) == 1 && _in->get_size(d) != 1 ){
	throw std::runtime_error( "downsample(): uneven array dimensions larger than one not accepted");
      }
    }
    
    typename uint64d<D>::Type matrix_size_in = from_std_vector<size_t,D>( *_in->get_dimensions() );
    typename uint64d<D>::Type matrix_size_out = matrix_size_in >> 1;

    for( size_t d=0; d<D; d++ ){
      if( matrix_size_out[d] == 0 ) 
	matrix_size_out[d] = 1;
    }
  
    size_t num_elements = prod(matrix_size_out);
    size_t num_batches = 1;

    for( size_t d=D; d<_in->get_number_of_dimensions(); d++ ){
      num_batches *= _in->get_size(d);
    }
  
    std::vector<size_t> dims = to_std_vector(matrix_size_out);
    for( size_t d=D; d<_in->get_number_of_dimensions(); d++ ){
      dims.push_back(_in->get_size(d));
    }
  
    REAL *in = _in->get_data_ptr();

    boost::shared_ptr< hoNDArray<REAL> > _out( new hoNDArray<REAL>(&dims) );
    REAL *out = _out->get_data_ptr();
    
    typedef vector_td<size_t,D> uint64d;

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for( long long idx=0; idx < num_elements*num_batches; idx++ ){

      const size_t frame_offset = idx/num_elements;
      const uint64d co_out = idx_to_co<D>( idx-frame_offset*num_elements, matrix_size_out );
      const uint64d co_in = co_out << 1;
      const uint64d twos(2);
      const size_t num_adds = 1 << D;

      size_t actual_adds = 0;
      REAL res = REAL(0);

      for( size_t i=0; i<num_adds; i++ ){
	const uint64d local_co = idx_to_co<D>( i, twos );
	if( weak_greater_equal( local_co, matrix_size_out ) ) continue; // To allow array dimensions of size 1
	const size_t in_idx = co_to_idx<D>(co_in+local_co, matrix_size_in)+frame_offset*prod(matrix_size_in);
	actual_adds++;
	res += in[in_idx];
      }    
      out[idx] = res/REAL(actual_adds);
    }

    return _out;
  }

  // Linear interpolation upsampling
  template<class REAL, unsigned int D> boost::shared_ptr< hoNDArray<REAL> >
  upsample( hoNDArray<REAL> *_in )
  {
    // A few sanity checks 

    if( _in == 0x0 ){
      throw std::runtime_error("upsample(): illegal input provided.");
    }

    if( _in->get_number_of_dimensions() < D ){
      throw std::runtime_error( "upsample(): the number of array dimensions should be at least D");
    }
    
    typename uint64d<D>::Type matrix_size_in = from_std_vector<size_t,D>( *_in->get_dimensions() );
    typename uint64d<D>::Type matrix_size_out = matrix_size_in << 1;

    for( size_t d=0; d<D; d++ ){
      if( matrix_size_in[d] == 1 )
	matrix_size_out[d] = 1;
    }
  
    size_t num_elements = prod(matrix_size_out);
    size_t num_batches = 1;

    for( size_t d=D; d<_in->get_number_of_dimensions(); d++ ){
      num_batches *= _in->get_size(d);
    }
  
    std::vector<size_t> dims = to_std_vector(matrix_size_out);
    for( size_t d=D; d<_in->get_number_of_dimensions(); d++ ){
      dims.push_back(_in->get_size(d));
    }

    REAL *in = _in->get_data_ptr();

    boost::shared_ptr< hoNDArray<REAL> > _out( new hoNDArray<REAL>(&dims) );
    REAL *out = _out->get_data_ptr();
    
    typedef vector_td<size_t,D> uint64d;

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for( long long idx=0; idx < num_elements*num_batches; idx++ ){
      
      REAL res = REAL(0);

      const size_t num_neighbors = 1 << D;
      const size_t frame_idx = idx/num_elements;
      const uint64d co_out = idx_to_co<D>( idx-frame_idx*num_elements, matrix_size_out );

      // We will only proceed if all neighbours exist (this adds a zero-boundary to the upsampled image/vector field)
      //
    
      if( !is_border_pixel<REAL,D>(co_out, matrix_size_out) ){
      
	for( size_t i=0; i<num_neighbors; i++ ){
	
	  // Determine coordinate of neighbor in input
	  //

	  const uint64d twos(2);
	  const uint64d stride = idx_to_co<D>( i, twos );

	  if( weak_greater_equal( stride, matrix_size_out ) ) continue; // To allow array dimensions of 1

	  // Be careful about dimensions of size 1
	  uint64d ones(1);
	  for( size_t d=0; d<D; d++ ){
	    if( matrix_size_out[d] == 1 )
	      ones[d] = 0;
	  }
	  uint64d co_in = ((co_out-ones)>>1)+stride;
	
	  // Read corresponding pixel value
	  //
	
	  const size_t in_idx = co_to_idx<D>(co_in, matrix_size_in)+frame_idx*prod(matrix_size_in);
	  REAL value = in[in_idx];
	
	  // Determine weight
	  //
	
	  REAL weight = REAL(1);
	
	  for( size_t dim=0; dim<D; dim++ ){	  
	    if( matrix_size_in[dim] > 1 ){
	      if( stride.vec[dim] == (co_out.vec[dim]%2) ) {
		weight *= REAL(0.25);
	      }
	      else{
		weight *= REAL(0.75);
	      }
	    }
	  }
	
	  // Accumulate result
	  //
	
	  res += weight*value;
	}
      }
      out[idx] = res;
    }
    
    return _out;
  }

}
