#pragma once

#include <boost/make_shared.hpp>
#include <boost/range/combine.hpp>
#include <numeric>
#include "hoNDArray.h"
#include "hoNDArray_iterators.h"
#include "vector_td_utilities.h"

#include <boost/version.hpp>

#if (BOOST_VERSION < 107200)
#include <boost/math/interpolators/cubic_b_spline.hpp>
namespace boost::math::interpolators {
    auto cardinal_cubic_b_spline = [](auto ... args){return boost::math::cubic_b_spline(args...);};
}
#else
#include <boost/math/interpolators/cardinal_cubic_b_spline.hpp>
#endif
#include <boost/math/special_functions/trunc.hpp>
#include <boost/range/adaptor/strided.hpp>
#include <range/v3/numeric.hpp>
#include <range/v3/view.hpp>
#include <range/v3/action.hpp>

#ifdef USE_OMP
#include <omp.h>
#endif
#ifdef max
#undef max
#endif

#ifdef min
#undef min
#endif

namespace Gadgetron {
  class ArrayIterator
  {
  public:

    ArrayIterator(std::vector<size_t> *dimensions, std::vector<size_t> *order)
    {
      block_sizes_.push_back(1);
      for (size_t i = 0; i < order->size(); i++) {
        dimensions_.push_back((*dimensions)[i]);
        order_.push_back((*order)[i]);
        current_.push_back(0);
        if (i > 0) {
          block_sizes_.push_back(block_sizes_[i-1]*dimensions_[i-1]);
        }
      }
      current_idx_ = 0;
    }

    inline size_t advance()
    {
      size_t order_index = 0;
      current_[order_[order_index]]++;
      while (current_[order_[order_index]] >= dimensions_[order_[order_index]]) {
        current_[order_[order_index]] = 0;
        order_index = (order_index+1)%dimensions_.size();
        current_[order_[order_index]]++;
      }

      current_idx_ = 0;
      for (size_t i = 0; i < dimensions_.size(); i++) {
        current_idx_ += current_[i]*block_sizes_[i];
      }
      return current_idx_;
    }

    inline size_t get_current_idx() const {
      return current_idx_;
    }

    std::vector<size_t> get_current_sub() {
      return current_;
    }

  protected:
    std::vector<size_t> dimensions_;
    std::vector<size_t> order_;
    std::vector<size_t> current_;
    std::vector<size_t> block_sizes_;
    size_t current_idx_;
  };

  template<class T> hoNDArray<T> shift_dim( const hoNDArray<T>& in, int shift )
  {
    std::vector<size_t> order;
    for (size_t i = 0; i < in.get_number_of_dimensions(); i++) {
      order.push_back(static_cast<size_t>((i+shift)%in.get_number_of_dimensions()));
    }
    return permute(in,order);
  }

  template<class T> void shift_dim(const hoNDArray<T>& in, hoNDArray<T>& out, int shift )
  {
    std::vector<size_t> order;
    for (size_t i = 0; i < in.get_number_of_dimensions(); i++) {
      order.push_back(static_cast<size_t>((i+shift)%in.get_number_of_dimensions()));
    }
    permute(in,out,order);
  }

  template<class T>  hoNDArray<T>
  permute( const hoNDArray<T>& in, const std::vector<size_t>& dim_order)
  {

    std::vector<size_t> dims;
    for (size_t i = 0; i < dim_order.size(); i++)
      dims.push_back(in.get_dimensions()->at(dim_order[i]));
    hoNDArray<T> out(dims);
    permute( in, out, dim_order);
    return out;
  }

  template<class T> void
  permute(const  hoNDArray<T>& in, hoNDArray<T>& out, const std::vector<size_t>& dim_order)
  {

    // Check ordering array
    if (dim_order.size() > in.get_number_of_dimensions()) {
      throw std::runtime_error("hoNDArray::permute - Invalid length of dimension ordering array");;
    }

    std::vector<size_t> dim_count(in.get_number_of_dimensions(),0);
    for (size_t i = 0; i < dim_order.size(); i++) {
      if (dim_order[i] >= in.get_number_of_dimensions()) {
        throw std::runtime_error("hoNDArray::permute - Invalid dimension order array");;
      }
      dim_count[dim_order[i]]++;
    }

    // Create an internal array to store the dimensions
    std::vector<size_t> dim_order_int;

    // Check that there are no duplicate dimensions
    for (size_t i = 0; i < dim_order.size(); i++) {
      if (dim_count[dim_order[i]] != 1) {
        throw std::runtime_error("hoNDArray::permute - Invalid dimension order array (duplicates)");;

      }
      dim_order_int.push_back(dim_order[i]);
    }

    for (size_t i = 0; i < dim_order_int.size(); i++) {
      if ((*in.get_dimensions())[dim_order_int[i]] != out.get_size(i)) {
        throw std::runtime_error("permute(): dimensions of output array do not match the input array");;
      }
    }

    // Pad dimension order array with dimension not mentioned in order array
    if (dim_order_int.size() < in.get_number_of_dimensions()) {
      for (size_t i = 0; i < dim_count.size(); i++) {
        if (dim_count[i] == 0) {
          dim_order_int.push_back(i);
        }
      }
    }

    T* o = out.get_data_ptr();

    // if memcpy can be used during permute
    size_t stride = 1;
    size_t num_dim_memcpy = 0;
    for (size_t i = 0; i < dim_order_int.size(); i++) {
        if (dim_order_int[i]==i){
            stride *= in.get_size(i);
            num_dim_memcpy = i;
        }
        else{
            break;
        }
    }

    if (stride == 1) {
        // point by point assignment is needed
        ArrayIterator it(in.get_dimensions().get(), &dim_order_int);
        for (size_t i = 0; i < in.get_number_of_elements(); i++) {
            o[i] = in.get_data_ptr()[it.get_current_idx()];
            it.advance();
        }
    }
    else {
        // memcpy can be used

        size_t nDim = in.get_number_of_dimensions();
        size_t num_memcpy = in.get_number_of_elements() / stride;

        if (num_dim_memcpy == nDim - 1){
            memcpy(out.begin(), in.begin(), in.get_number_of_bytes());
            return;
        }

        // for the array index calculation
        std::vector<size_t> dim_permute(nDim-num_dim_memcpy-1);
        for (size_t i = num_dim_memcpy+1; i < dim_order_int.size(); i++) {
            dim_permute[i - num_dim_memcpy - 1] = in.get_size(i);
        }

        size_t n;

        const hoNDArray<T> permuteArray(dim_permute, const_cast<T*>(in.get_data_ptr()), false);

        // starting index for in and out array for every permute memcpy operation
        std::vector<size_t> ind_permute_in(dim_permute.size(), 0), ind_in(nDim, 0), ind_out(nDim, 0);

        for (n = 0; n < num_memcpy; n++) {
            permuteArray.calculate_index(n, ind_permute_in);
            memcpy(&ind_in[0] + num_dim_memcpy + 1, &ind_permute_in[0], sizeof(size_t)*ind_permute_in.size());

            // permute the indexes
            for (size_t i = 0; i < nDim; i++) {
                ind_out[i] = ind_in[dim_order_int[i]];
            }

            size_t offset_in = in.calculate_offset(ind_in);
            size_t offset_out = out.calculate_offset(ind_out);

            memcpy(o + offset_out, in.begin() + offset_in, sizeof(T)*stride);
        }
    }
  }

  // Expand array to new dimension
  template<class T> hoNDArray<T>
  expand(const hoNDArray<T>& in, size_t new_dim_size )
  {

    const size_t number_of_elements_in = in.get_number_of_elements();

    std::vector<size_t> dims = in.dimensions();
    dims.push_back(new_dim_size);

    auto out = hoNDArray<T>(dims);

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for( long long int idx=0; idx<number_of_elements_in*new_dim_size; idx++ ){
      out[idx] = in[idx%number_of_elements_in];
    }
    return out;
  }

  namespace {
      template<class T, class ACCUMULATOR> hoNDArray<T>
      accumulate(const hoNDArray<T>& in, size_t dim, ACCUMULATOR acc )
      {
          if( !(in.get_number_of_dimensions()>1) ){
              throw std::runtime_error("sum(): underdimensioned.");;
          }

          if( dim > in.get_number_of_dimensions()-1 ){
              throw std::runtime_error( "sum(): dimension out of range.");;
          }

          size_t number_of_batches = in.get_size(dim);
          size_t number_of_elements = in.get_number_of_elements()/number_of_batches;
          std::vector<size_t> dims;
          for (auto i = 0; i < in.get_number_of_dimensions(); i++){
              if (i != dim) dims.push_back(in.get_size(i));
          }

          auto  out = hoNDArray<T>(dims);
          auto orig_dims = *in.get_dimensions();
          auto stride = std::accumulate(orig_dims.begin(),orig_dims.begin()+dim,1,std::multiplies<size_t>());



          size_t inner_elements = stride;
          size_t outer_elements = out.get_number_of_elements()/inner_elements;
//#ifdef USE_OMP
//#pragma omp parallel for schedule(dynamic,1) collapse(2)
//#endif
          for (size_t outer_idx = 0; outer_idx < outer_elements; outer_idx++) {

              for (size_t idx = 0; idx < inner_elements; idx++) {
                  size_t offset = outer_idx*inner_elements;
                  size_t old_offset = offset*number_of_batches;
                  T val = in.at(idx+old_offset);
                  for (size_t j = 1; j < number_of_batches; j++) {
                      size_t in_idx = j * stride + idx+ old_offset;
                      val = acc(val,in.at(in_idx));
                  }
                  out.at(idx + offset) = val;
              }
          }
          return out;
      }
  }

  // Sum over dimension
  template<class T> hoNDArray<T>
  sum(const hoNDArray<T>& in, size_t dim )
  {
      return accumulate(in, dim, std::plus<T>());
  }
    template<class T> boost::shared_ptr<hoNDArray<T>>
  sum(const hoNDArray<T>* in, size_t dim )
  {
      return boost::make_shared<hoNDArray<T>>(accumulate(*in, dim, std::plus<T>()));
  }

    template<class T> hoNDArray<T>
    max(const hoNDArray<T>& in, size_t dim )
    {
        return accumulate(in, dim, [](auto v1, auto v2){ return std::max(v1,v2);});
    }

    template<class T> hoNDArray<T>
    min(const hoNDArray<T>& in, size_t dim )
    {
        return accumulate(in, dim,  [](auto v1, auto v2){ return std::min(v1,v2);});
    }

    /**
  * @param[in] crop_offset starting position to crop
  * @param[in] crop_size Size of cropped array
  * @param[in] in input array
  * @param[out] out Output array after cropping
  */
  template<class T, unsigned int D> void
  crop(const vector_td<size_t, D>& crop_offset, const vector_td<size_t, D>& crop_size, const hoNDArray<T>& in, hoNDArray<T>& out)
  {

      if (in.get_number_of_dimensions() < D){
          std::stringstream ss;
          ss << "crop: number of image dimensions should be at least " << D;
          throw std::runtime_error(ss.str());;
      }

      std::vector<size_t> dims = to_std_vector(crop_size);
      for (unsigned int d = D; d<in.get_number_of_dimensions(); d++){
          dims.push_back(in.get_size(d));
      }

      if (!out.dimensions_equal(&dims)){
          out.create(dims);
      }

      typename uint64d<D>::Type matrix_size_in = from_std_vector<size_t, D>(*in.get_dimensions());
      typename uint64d<D>::Type matrix_size_out = from_std_vector<size_t, D>(*out.get_dimensions());

      if (weak_greater(crop_offset + matrix_size_out, matrix_size_in)){
          throw std::runtime_error("crop: cropping size mismatch");;
      }

      size_t len = out.get_size(0);
      size_t num = out.get_number_of_elements() / len;

      long long k;

      const T *in_ptr = in.get_data_ptr();
      T *out_ptr = out.get_data_ptr();

      #pragma omp parallel default(none) private(k) shared(in_ptr, out_ptr, num, len, in, out, crop_offset)
      {
          std::vector<size_t> ind;

      #pragma omp for
          for (k = 0; k < (long long)num; k++){
              ind = out.calculate_index(k*len);
              for (unsigned int d = 0; d < D; d++){
                  ind[d] += crop_offset[d];
              }

              const T* in_ptr_curr = in_ptr + in.calculate_offset(ind);
              memcpy(out_ptr + k*len, in_ptr_curr, sizeof(T)*len);
          }
      }
  }

  /**
  * @param[in] crop_size
  * @param[in] in input array
  * Crop the input array around its center N/2; that is, the center pixel of in array is the center pixel of out array
  */
  template<class T, unsigned int D> hoNDArray<T>
  crop(const vector_td<size_t, D>& crop_size, const hoNDArray<T>& in)
  {
    // compute crop offset, perserving the center
    hoNDArray<T> out;
    auto crop_offset = (from_std_vector<size_t,D>(*in.get_dimensions())-crop_size)/size_t(2);
    crop(crop_offset, crop_size, in, out);
    return out;
  }

  template<class T> void
  crop(size_t x, const hoNDArray<T>& in, hoNDArray<T>& out)
  {
      vector_td<size_t, 1> crop_size(x);

      auto crop_offset = (from_std_vector<size_t,1>(*in.get_dimensions())-crop_size)/size_t(2);
      crop(crop_offset, crop_size, in, out);
  }

  template<class T> void
  crop(size_t x, size_t y, const hoNDArray<T>& in, hoNDArray<T>& out)
  {
      vector_td<size_t, 2> crop_size(x, y);

      auto crop_offset = (from_std_vector<size_t,2>(*in.get_dimensions())-crop_size)/size_t(2);
      crop(crop_offset,crop_size, in, out);
  }

  template<class T> void
  crop(size_t x, size_t y, size_t z, const hoNDArray<T>& in, hoNDArray<T>& out)
  {
      vector_td<size_t, 3> crop_size(x, y, z);

      auto crop_offset = (from_std_vector<size_t,3>(*in.get_dimensions())-crop_size)/size_t(2);
      crop(crop_offset, crop_size, in, out);
  }

  template<class T, unsigned int D> hoNDArray<T>
  crop( const vector_td<size_t, D>& crop_offset, const vector_td<size_t, D>& crop_size, const hoNDArray<T>& in )
  {
    auto out =  hoNDArray<T>();
    crop(crop_offset, crop_size, in, out);
    return out;
  }

  /**
  * @param[in] offset_src starting position in src array
  * @param[in] size Size of subarray to be replaced
  * @param[in] src Src array to read in replaced content
  * @param[in] offset_dst starting position in dst array
  * @param[out] dst array to be replaced; other part outside the offset+size region will be unchanged
  */
  template<class T, unsigned int D> void
  fill(const vector_td<size_t, D>& offset_src, const vector_td<size_t, D>& size, hoNDArray<T> *src, const vector_td<size_t, D>& offset_dst, hoNDArray<T> *dst)
  {
      if (src == 0x0) {
          throw std::runtime_error("replace: 0x0 src array provided");;
      }

      if (src->get_number_of_dimensions() < D)
      {
          std::stringstream ss;
          ss << "fill: number of src image dimensions should be at least " << D;
          throw std::runtime_error(ss.str());;
      }

      if (dst == 0x0)
      {
          throw std::runtime_error("replace: 0x0 dst array provided");;
      }

      if (dst->get_number_of_dimensions() < D)
      {
          std::stringstream ss;
          ss << "fill: number of dst image dimensions should be at least " << D;
          throw std::runtime_error(ss.str());;
      }

      if (src->get_number_of_dimensions() != dst->get_number_of_dimensions())
      {
          std::stringstream ss;
          ss << "fill: src and dst array have different number of dimensions " << D;
          throw std::runtime_error(ss.str());;
      }

      std::vector<size_t> src_dim;
      src->get_dimensions(src_dim);

      std::vector<size_t> dst_dim;
      dst->get_dimensions(dst_dim);

      size_t d;
      for (d = 0; d < D; d++)
      {
          if (src_dim[d] < offset_src[d]+size[d]-1)
          {
              throw std::runtime_error("fill: src array is too small for provided offset and size");;
          }

          if (dst_dim[d] < offset_dst[d] + size[d] - 1)
          {
              throw std::runtime_error("fill: dst array is too small for provided offset and size");;
          }
      }

      size_t len = size[0];
      size_t num = 1;

      for (d = 1; d < D; d++) num *= size[d];

      long long k;

      T *src_ptr = src->get_data_ptr();
      T *dst_ptr = dst->get_data_ptr();

      std::vector<size_t> size_dim = to_std_vector(size);
      hoNDArray<T> array_size;
      array_size.create(size_dim, src->begin());

      {
          std::vector<size_t> ind_src = src->calculate_index(0);
          std::vector<size_t> ind_dst = dst->calculate_index(0);

          std::vector<size_t> ind_size(D, 0);

          for (k = 0; k < (long long)num; k++)
          {
              ind_size = array_size.calculate_index(k*len);

              for (unsigned int d = 0; d < D; d++)
              {
                  ind_src[d] = offset_src[d] + ind_size[d];
                  ind_dst[d] = offset_dst[d] + ind_size[d];
              }

              T* src_ptr_curr = src_ptr + src->calculate_offset(ind_src);
              T* dst_ptr_curr = dst_ptr + dst->calculate_offset(ind_dst);

              memcpy(dst_ptr_curr, src_ptr_curr, sizeof(T)*len);
          }
      }
  }

  template<class T, unsigned int D> void
  fill(const vector_td<size_t, D>& offset_src, hoNDArray<T>& src, const vector_td<size_t, D>& offset_dst, hoNDArray<T>& dst)
  {
      std::vector<size_t> dim;
      src.get_dimensions(dim);

      vector_td<size_t, D> size;

      if (dim.size() < D)
      {
          std::stringstream ss;
          ss << "fill: number of src image dimensions should be at least " << D;
          throw std::runtime_error(ss.str());;
      }

      size_t d;
      for (d = 0; d < D; d++) size[d] = dim[d];

      Gadgetron::fill(offset_src, size, &src, offset_dst, &dst);
  }

  template<class T, unsigned int D> void
  fill(hoNDArray<T>& src, const vector_td<size_t, D>& offset_dst, hoNDArray<T>& dst)
  {
      std::vector<size_t> dim;
      src.get_dimensions(dim);

      vector_td<size_t, D> offset_src, size;

      if (dim.size() < D)
      {
          std::stringstream ss;
          ss << "fill: number of src image dimensions should be at least " << D;
          throw std::runtime_error(ss.str());;
      }

      size_t d;
      for (d = 0; d < D; d++)
      {
          offset_src[d] = 0;
          size[d] = dim[d];
      }

      Gadgetron::fill(offset_src, size, &src, offset_dst, &dst);
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
  pad(const typename uint64d<D>::Type& size, const hoNDArray<T>& in, hoNDArray<T>& out, bool preset_out_with_val = true, T val = T(0))
  {

      if (in.get_number_of_dimensions() < D){
          std::stringstream ss;
          ss << "pad: number of image dimensions should be at least " << D;
          throw std::runtime_error(ss.str());;
      }

      unsigned int d;

      std::vector<size_t> dims = to_std_vector(size);
      for (d = D; d<in.get_number_of_dimensions(); d++){
          dims.push_back(in.get_size(d));
      }

      if (!out.dimensions_equal(&dims)){
          out.create(dims);
      }

      if (in.dimensions_equal(&dims)){
          memcpy(out.begin(), in.begin(), in.get_number_of_bytes());
          return;
      }

      const T *in_ptr = in.get_data_ptr();
      T *out_ptr = out.get_data_ptr();

      if (preset_out_with_val){
          if (val == T(0)){
              memset(out_ptr, 0, out.get_number_of_bytes());
          }
          else{
                size_t N = out.get_number_of_elements();
                long long n;
                #pragma omp parallel for default(none) private(n) shared(N, out_ptr, val)
                for (n = 0; n<(long long)N; n++)
                {
                    out_ptr[n] = val;
                }
          }
      }

      typename uint64d<D>::Type matrix_size_in = from_std_vector<size_t, D>(*in.get_dimensions());
      typename uint64d<D>::Type matrix_size_out = from_std_vector<size_t, D>(*out.get_dimensions());

      if (weak_greater(matrix_size_in, matrix_size_out)){
          throw std::runtime_error("pad: size mismatch, cannot expand");
      }

      typename uint64d<D>::Type offset(D);
      for (d = 0; d<D; d++){
          offset[d] = matrix_size_out[d]/2 - matrix_size_in[d]/2;
      }

      size_t len = in.get_size(0);
      size_t num = in.get_number_of_elements() / len;

      long long k;

#pragma omp parallel default(none) private(k, d) shared(in_ptr, out_ptr, num, len, in, out, offset)
      {
          std::vector<size_t> ind;

#pragma omp for
          for (k = 0; k < (long long)num; k++){
              ind = in.calculate_index(k*len);
              for (d = 0; d < D; d++){
                  ind[d] += offset[d];
              }

              T* out_ptr_curr = out_ptr + out.calculate_offset(ind);
              memcpy(out_ptr_curr, in_ptr + k*len, sizeof(T)*len);
          }
      }
  }

  template<class T, unsigned int D> void pad(const hoNDArray<T>& in, hoNDArray<T>& out, T val = T(0)){
        vector_td<size_t,D> dims = from_std_vector<size_t,D>(*out.get_dimensions());
        pad<T,D>(dims,in,out,true, val);
  }

  template<class T> void
  pad(size_t x, const hoNDArray<T>& in, hoNDArray<T>& out, bool preset_out_with_val = true, T val = T(0))
  {
      typename uint64d<1>::Type padSize(x);
      pad<T, 1>(padSize, in, out, preset_out_with_val, val);
  }

  template<class T> void
  pad(size_t x, size_t y, const hoNDArray<T>& in, hoNDArray<T>& out, bool preset_out_with_val = true, T val = T(0))
  {
      typename uint64d<2>::Type padSize(x, y);
      pad<T, 2>(padSize, in, out, preset_out_with_val, val);
  }

  template<class T> void
  pad(size_t x, size_t y, size_t z, const hoNDArray<T> &in, hoNDArray<T>& out, bool preset_out_with_val = true, T val = T(0))
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
  template<class T, unsigned int D>  hoNDArray<T>
  pad(const typename uint64d<D>::Type& size, const hoNDArray<T> & in, T val = T(0))
  {
    auto out = hoNDArray<T>();
    pad<T,D>(size, in, out, true, val);
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
   hoNDArray<REAL>  downsample(const  hoNDArray<REAL>& _in )
  {
    // A few sanity checks
    if( _in.get_number_of_dimensions() < D ){
      throw std::runtime_error( "downsample(): the number of array dimensions should be at least D");
    }

    for( size_t d=0; d<D; d++ ){
      if( (_in.get_size(d)%2) == 1 && _in.get_size(d) != 1 ){
	throw std::runtime_error( "downsample(): uneven array dimensions larger than one not accepted");
      }
    }

    typename uint64d<D>::Type matrix_size_in = from_std_vector<size_t,D>( *_in.get_dimensions() );
    typename uint64d<D>::Type matrix_size_out = matrix_size_in >> 1;

    for( size_t d=0; d<D; d++ ){
      if( matrix_size_out[d] == 0 )
	matrix_size_out[d] = 1;
    }

    size_t num_elements = prod(matrix_size_out);
    size_t num_batches = 1;

    for( size_t d=D; d<_in.get_number_of_dimensions(); d++ ){
      num_batches *= _in.get_size(d);
    }

    std::vector<size_t> dims = to_std_vector(matrix_size_out);
    for( size_t d=D; d<_in.get_number_of_dimensions(); d++ ){
      dims.push_back(_in.get_size(d));
    }

    const REAL *in = _in.get_data_ptr();

     hoNDArray<REAL>  _out( dims );
    REAL *out = _out.get_data_ptr();

    typedef vector_td<size_t,D> uint64d;

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for( int64_t idx=0; idx < num_elements*num_batches; idx++ ){

      const size_t frame_offset = idx/num_elements;
      const uint64d co_out = idx_to_co<uint64_t,D>( idx-frame_offset*num_elements, matrix_size_out );
      const uint64d co_in = co_out << 1;
      const uint64d twos(2);
      const size_t num_adds = 1 << D;

      size_t actual_adds = 0;
      REAL res = REAL(0);

      for( size_t i=0; i<num_adds; i++ ){
	const uint64d local_co = idx_to_co( i, twos );
	if( weak_greater_equal( local_co, matrix_size_out ) ) continue; // To allow array dimensions of size 1
	const size_t in_idx = co_to_idx(co_in+local_co, matrix_size_in)+frame_offset*prod(matrix_size_in);
	actual_adds++;
	res += in[in_idx];
      }
      out[idx] = res/REAL(actual_adds);
    }

    return _out;
  }

  namespace {
          template<class T> hoNDArray<T> upsample_along_dimension(const hoNDArray<T>& array,int dim){
              auto new_dims = *array.get_dimensions();
              auto old_dim = new_dims[dim];
              new_dims[dim] *= 2;
              hoNDArray<T> result(new_dims);
              size_t stride = std::accumulate(new_dims.begin(),new_dims.begin()+dim,size_t(1),std::multiplies<size_t>());

              size_t nbatches = result.get_number_of_elements()/stride/new_dims[dim];
              size_t batch_size = stride*new_dims[dim];
              size_t old_batch_size = batch_size/2;

#pragma omp parallel for
              for (int batch = 0; batch < nbatches; batch++){
                  T* result_ptr = result.get_data_ptr()+batch_size*batch;
                  const T* input_ptr = array.get_data_ptr()+batch*old_batch_size;
                  for (size_t i = 0; i < old_dim-1; i++){
                      for (size_t k = 0; k < stride; k++){
                        result_ptr[2*i*stride+k] = input_ptr[i*stride+k];
                        result_ptr[(2*i+1)*stride+k] = (input_ptr[i*stride+k]+input_ptr[i*stride+k])/2;
                      }
                  }

                  size_t i = old_dim-1;
                  for (size_t k = 0; k < stride; k++){
                    result_ptr[2*i*stride+k] = input_ptr[i*stride+k];
                    result_ptr[(2*i+1)*stride+k] = input_ptr[i*stride+k];
                  }

              }
              return result;



          }


          template<class T> hoNDArray<T> upsample_spline_along_dimension(const hoNDArray<T>& array,int dim,int scale){
              namespace ba = boost::adaptors;
              namespace bm = boost::math;
              auto new_dims = *array.get_dimensions();
              auto old_dim = new_dims[dim];
              new_dims[dim] *= 2;
              hoNDArray<T> result(new_dims);
              size_t stride = std::accumulate(new_dims.begin(),new_dims.begin()+dim,size_t(1),std::multiplies<size_t>());
              size_t nbatches = result.get_number_of_elements()/stride/new_dims[dim];
              size_t batch_size = stride*new_dims[dim];
              size_t old_batch_size = batch_size/2;

#pragma omp parallel for
              for (int batch = 0; batch < (int)nbatches; batch++){
                  T* result_ptr = result.get_data_ptr()+batch_size*batch;
                  const T* input_ptr = array.get_data_ptr()+batch*old_batch_size;

                  for (size_t k = 0; k < stride; k++){
                      auto strided_iterator = std::make_pair(input_ptr+k,input_ptr+k+old_batch_size) | ba::strided(stride);
                      auto spline = bm::interpolators::cardinal_cubic_b_spline(
                          boost::begin(strided_iterator),
                          boost::end(strided_iterator),
                          T(0.25)*scale, T(scale), T(0), T(0)
                      );
                      for (int i = 0; i < new_dims[dim]; i++){
                          result_ptr[k+i*stride] = spline(i);
                      }

                  }

              }
              return result;



          }
      }
  // Linear interpolation upsampling
  template<class T, unsigned int D> hoNDArray<T>
  upsample( const  hoNDArray<T>& in )
  {


    if( in.get_number_of_dimensions() < D ){
      throw std::runtime_error( "upsample(): the number of array dimensions should be at least D");
    }

    hoNDArray<T> result = in;
    for (int i = D-1; i >= 0; i--){
        result = upsample_along_dimension<T>(result,i);
    }
    return result;

  }

    template<class T, unsigned int D> hoNDArray<T>
  upsample_spline( const  hoNDArray<T>& in, int scale = 2 )
  {

    if( in.get_number_of_dimensions() < D ){
      throw std::runtime_error( "upsample(): the number of array dimensions should be at least D");
    }

    hoNDArray<T> result = in;
    for (int i = D-1; i >= 0; i--){
        result = upsample_spline_along_dimension<T>(result,i,scale);
    }
    return result;

  }


  // Linear interpolation upsampling
  template<class T, unsigned int D> hoNDArray<T>
  upsample_nearest( const  hoNDArray<T>& in )
  {
    // A few sanity checks

    if( in.get_number_of_dimensions() < D ){
      throw std::runtime_error( "upsample(): the number of array dimensions should be at least D");
    }

    typename uint64d<D>::Type matrix_size_in = from_std_vector<size_t,D>( *in.get_dimensions() );
    typename uint64d<D>::Type matrix_size_out = matrix_size_in << 1;

    for( size_t d=0; d<D; d++ ){
      if( matrix_size_in[d] == 1 )
	matrix_size_out[d] = 1;
    }

    size_t num_elements = prod(matrix_size_out);
    size_t num_batches = 1;

    for( size_t d=D; d<in.get_number_of_dimensions(); d++ ){
      num_batches *= in.get_size(d);
    }

    std::vector<size_t> dims = to_std_vector(matrix_size_out);
    for( size_t d=D; d<in.get_number_of_dimensions(); d++ ){
      dims.push_back(in.get_size(d));
    }

    const T *in_ptr = in.get_data_ptr();

    hoNDArray<T> out(&dims);
    T *out_ptr = out.get_data_ptr();

    typedef vector_td<size_t,D> uint64d;

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for( long long idx=0; idx < num_elements*num_batches; idx++ ){

      const size_t frame_idx = idx/num_elements;
      const uint64d co_out = idx_to_co<D>( idx-frame_idx*num_elements, matrix_size_out );
      uint64d co_in = co_out/uint64_t(2);

      const size_t in_idx = co_to_idx<D>(co_in, matrix_size_in)+frame_idx*prod(matrix_size_in);
      out_ptr[idx] = in_ptr[in_idx];
	}

    return out;
  }

  template<class T> hoNDArray<T> repeat(const hoNDArray<T>& array,unsigned int repeats){
      auto dims = array.dimensions();
      dims.push_back(repeats);

      hoNDArray<T> output(dims);

      for (auto span : spans(output, array.get_number_of_dimensions())) {
          span = array;
      }

      return output;
  }

  /**
   * This functions takes a collection of hoNDArrays and concatenates them along the specified dimension
   * @tparam COLL Collection of hoNDArray such as std::vector<hoNDArray<float>>.
   * @param arrays The hoNDArrays. Must be of equal size, except along the concat dimension
   * @param dimension Dimension along which to concatenate
   * @return The concatenated arrays.
   */

  template <class COLL> auto concat_along_dimension(const COLL& arrays, size_t dimension) {
      using namespace ranges;
      using T = std::decay_t<decltype(*std::begin(*std::begin(arrays)))>;
      if (arrays.empty())
          return hoNDArray<T>();

      const hoNDArray<T>& first = *std::begin(arrays);
      std::vector dims = first.dimensions();

      size_t count = ranges::accumulate(arrays | views::transform([dimension](const auto& array) {
                                            return array.dimensions().at(dimension);
                                        }),
                                        size_t(0));
      dims[dimension] = count;

      auto dimensions_valid = [&dims,dimension](const auto& array) {
          bool result = true;
          const auto& d = array.dimensions();
          for (size_t i = 0; i < d.size(); i++) {
              if (i == dimension)
                  continue;
              result &= d[i] == dims[i];
          }
          return result && (d.size() == dims.size());
      };

      bool all_dimensions_valid = ranges::accumulate(arrays | views::transform(dimensions_valid), true, std::logical_and() );
      if (!all_dimensions_valid) throw std::runtime_error("The dimensions of all provided arrays must be equal except along the concatenate dimension");

      auto result = hoNDArray<T>(dims);

      const size_t inner_stride = ranges::accumulate(dims | views::slice(size_t(0), dimension),
                                                     size_t(1), std::multiplies());
      const size_t outer_stride = inner_stride * count;
      size_t current_slice = 0;

      for (const auto& array : arrays) {
          size_t slice_count = array.dimensions()[dimension];
          auto array_inner_stride = slice_count * inner_stride;
          auto repetitions = array.size() / array_inner_stride;

          for (int i = 0; i < repetitions; i++) {
              std::copy_n(array.begin() + i * array_inner_stride, array_inner_stride,
                          result.begin() + current_slice * inner_stride + outer_stride * i);
          }
          current_slice += slice_count;
      }
      return result;
  }
  template<class COLL>
  auto concat(const COLL &arrays) {

      using T = std::decay_t<decltype(*std::begin(*std::begin(arrays)))>;
      if (arrays.empty()) return hoNDArray<T>();

      const hoNDArray<T> &first = *std::begin(arrays);

      auto dims = first.dimensions();
      auto size = first.size();
      using std::begin;
      using std::end;

      if (!std::all_of(begin(arrays), end(arrays), [&](const auto &array) { return dims == array.dimensions(); }) ||
          !std::all_of(begin(arrays), end(arrays), [&](const auto &array) { return size == array.size(); })) {
          throw std::runtime_error("Array size or dimensions do not match.");
      }

      dims.push_back(arrays.size());
      hoNDArray<T> output(dims);

      auto output_iterator = spans(output, first.get_number_of_dimensions()).begin();

      for (const auto& array : arrays) {
          *output_iterator = array;
          ++output_iterator;
      }

      return output;
  }

  template<class T, class...  ARRAYS>
  hoNDArray<T> concat(const hoNDArray<T>& first_array, const ARRAYS& ... arrays){

      static_assert((std::is_same_v<hoNDArray<T>,std::decay_t<ARRAYS>> && ...));
      using namespace ranges;
      return concat(views::concat(views::single(first_array),views::single(arrays)...));
  }
}

