#include "hoNDArray_utils.h"
#include "hoNDArray_elemwise.h"
#include "complext.h"
#include "GadgetronTimer.h"

#include <gtest/gtest.h>
#include <complex>
#include <vector>

using namespace Gadgetron;
using testing::Types;

template <typename T> class hoNDArray_utils_TestReal : public ::testing::Test {
protected:
  virtual void SetUp() {
    size_t vdims[] = {37, 49, 23, 19}; //Using prime numbers for setup because they are messy
    dims = std::vector<size_t>(vdims,vdims+sizeof(vdims)/sizeof(size_t));
    Array = hoNDArray<T>(&dims);
  }
  std::vector<size_t> dims;
  hoNDArray<T> Array;
};

template <typename T> class hoNDArray_utils_TestCplx : public ::testing::Test {
protected:
  virtual void SetUp() {
    size_t vdims[] = {37, 49, 23, 19}; //Using prime numbers for setup because they are messy
    dims = std::vector<size_t>(vdims,vdims+sizeof(vdims)/sizeof(size_t));
    Array = hoNDArray<T>(&dims);
  }
  std::vector<size_t> dims;
  hoNDArray<T> Array;
};

template <typename T> class hoNDArray_utils_SpeedTestCplx : public ::testing::Test {
protected:
    virtual void SetUp() {
        size_t vdims[] = { 256, 144, 24, 192 };
        dims = std::vector<size_t>(vdims, vdims + sizeof(vdims) / sizeof(size_t));
        Array = hoNDArray<T>(&dims);
        Gadgetron::clear(Array);
    }
    std::vector<size_t> dims;
    hoNDArray<T> Array;
    hoNDArray<T> Array_out;
};

typedef Types<float, double> realImplementations;
typedef Types</*std::complex<float>, std::complex<double>,*/ float_complext, double_complext> cplxImplementations;

TYPED_TEST_CASE(hoNDArray_utils_TestReal, realImplementations);

TYPED_TEST(hoNDArray_utils_TestReal,permuteTest){

  fill(&this->Array,TypeParam(1));

  std::vector<size_t> order;
  order.push_back(0); order.push_back(1); order.push_back(2); order.push_back(3);

  this->Array.get_data_ptr()[37] = TypeParam(2);

  EXPECT_FLOAT_EQ(1, permute(&this->Array,&order)->at(0));
  EXPECT_FLOAT_EQ(2, permute(&this->Array,&order)->at(37));

  order.clear();
  order.push_back(1); order.push_back(0); order.push_back(2); order.push_back(3);

  EXPECT_FLOAT_EQ(2, permute(&this->Array,&order)->at(1));

  order.clear();
  order.push_back(3); order.push_back(1); order.push_back(2); order.push_back(0);

  EXPECT_FLOAT_EQ(2, permute(&this->Array,&order)->at(19));

  order.clear();
  order.push_back(2); order.push_back(0); order.push_back(1); order.push_back(3);

  EXPECT_FLOAT_EQ(2, permute(&this->Array,&order)->at(851));
}

TYPED_TEST(hoNDArray_utils_TestReal,shiftDimTest){

  fill(&this->Array,TypeParam(1));
  this->Array.get_data_ptr()[37] = 2;

  EXPECT_FLOAT_EQ(1, shift_dim(&this->Array,0)->at(0));
  EXPECT_FLOAT_EQ(2, shift_dim(&this->Array,0)->at(37));
  EXPECT_FLOAT_EQ(2, shift_dim(&this->Array,1)->at(1));
  EXPECT_FLOAT_EQ(2, shift_dim(&this->Array,-1)->at(37*19));
  EXPECT_FLOAT_EQ(2, shift_dim(&this->Array,2)->at(23*37*19));
  EXPECT_FLOAT_EQ(2, shift_dim(&this->Array,3)->at(37*19));
  EXPECT_FLOAT_EQ(2, shift_dim(&this->Array,4)->at(37));
}

TYPED_TEST(hoNDArray_utils_TestReal,sumTest){
  TypeParam v1 = TypeParam(12.34);
  unsigned int idx = 0;

  fill(&this->Array,v1);
  EXPECT_FLOAT_EQ(49*v1,sum(&this->Array,1)->get_data_ptr()[idx]);

  fill(&this->Array,v1);
  EXPECT_FLOAT_EQ(23*v1,sum(&this->Array,2)->get_data_ptr()[idx]);

  fill(&this->Array,v1);
  EXPECT_FLOAT_EQ(19*v1,sum(&this->Array,3)->get_data_ptr()[idx]);
}

TYPED_TEST_CASE(hoNDArray_utils_TestCplx, cplxImplementations);

TYPED_TEST(hoNDArray_utils_TestCplx,permuteTest){

  fill(&this->Array,TypeParam(1,1));

  std::vector<size_t> order;
  order.push_back(0); order.push_back(1); order.push_back(2); order.push_back(3);
  
  this->Array.get_data_ptr()[37] = TypeParam(2,3);

  EXPECT_FLOAT_EQ(1, real(permute(&this->Array,&order)->at(0)));
  EXPECT_FLOAT_EQ(1, imag(permute(&this->Array,&order)->at(0)));

  EXPECT_FLOAT_EQ(2, real(permute(&this->Array,&order)->at(37)));
  EXPECT_FLOAT_EQ(3, imag(permute(&this->Array,&order)->at(37)));

  order.clear();
  order.push_back(1); order.push_back(0); order.push_back(2); order.push_back(3);

  EXPECT_FLOAT_EQ(2, real(permute(&this->Array,&order)->at(1)));
  EXPECT_FLOAT_EQ(3, imag(permute(&this->Array,&order)->at(1)));

  order.clear();
  order.push_back(3); order.push_back(1); order.push_back(2); order.push_back(0);

  EXPECT_FLOAT_EQ(2, real(permute(&this->Array,&order)->at(19)));
  EXPECT_FLOAT_EQ(3, imag(permute(&this->Array,&order)->at(19)));

  order.clear();
  order.push_back(2); order.push_back(0); order.push_back(1); order.push_back(3);

  EXPECT_FLOAT_EQ(2, real(permute(&this->Array,&order)->at(851)));
  EXPECT_FLOAT_EQ(3, imag(permute(&this->Array,&order)->at(851)));

  order.clear();
  order.push_back(0); order.push_back(1); order.push_back(3); order.push_back(2);

  EXPECT_FLOAT_EQ(2, real(permute(&this->Array, &order)->at(37)));
  EXPECT_FLOAT_EQ(3, imag(permute(&this->Array, &order)->at(37)));

  order.clear();
  order.push_back(0); order.push_back(2); order.push_back(3); order.push_back(1);

  EXPECT_FLOAT_EQ(2, real(permute(&this->Array, &order)->at(37*23*19)));
  EXPECT_FLOAT_EQ(3, imag(permute(&this->Array, &order)->at(37*23*19)));
}

TYPED_TEST(hoNDArray_utils_TestCplx,shiftDimTest){

  fill(&this->Array,TypeParam(1,1));
  this->Array.get_data_ptr()[37]=TypeParam(2,3);

  EXPECT_FLOAT_EQ(1, real(shift_dim(&this->Array,0)->at(0)));
  EXPECT_FLOAT_EQ(1, imag(shift_dim(&this->Array,0)->at(0)));

  EXPECT_FLOAT_EQ(2, real(shift_dim(&this->Array,0)->at(37)));
  EXPECT_FLOAT_EQ(3, imag(shift_dim(&this->Array,0)->at(37)));

  EXPECT_FLOAT_EQ(2, real(shift_dim(&this->Array,1)->at(1)));
  EXPECT_FLOAT_EQ(3, imag(shift_dim(&this->Array,1)->at(1)));

  EXPECT_FLOAT_EQ(2, real(shift_dim(&this->Array,-1)->at(37*19)));
  EXPECT_FLOAT_EQ(3, imag(shift_dim(&this->Array,-1)->at(37*19)));

  EXPECT_FLOAT_EQ(2, real(shift_dim(&this->Array,2)->at(23*37*19)));
  EXPECT_FLOAT_EQ(3, imag(shift_dim(&this->Array,2)->at(23*37*19)));

  EXPECT_FLOAT_EQ(2, real(shift_dim(&this->Array,3)->at(37*19)));
  EXPECT_FLOAT_EQ(3, imag(shift_dim(&this->Array,3)->at(37*19)));

  EXPECT_FLOAT_EQ(2, real(shift_dim(&this->Array,4)->at(37)));
  EXPECT_FLOAT_EQ(3, imag(shift_dim(&this->Array,4)->at(37)));
}

TYPED_TEST(hoNDArray_utils_TestCplx,sumTest){
  TypeParam v1 = TypeParam(12.34, 56.78);
  unsigned int idx = 0;

  fill(&this->Array,v1);
  EXPECT_FLOAT_EQ(real(TypeParam(49)*v1),real(sum(&this->Array,1)->get_data_ptr()[idx]));
  EXPECT_FLOAT_EQ(imag(TypeParam(49)*v1),imag(sum(&this->Array,1)->get_data_ptr()[idx]));

  fill(&this->Array,v1);
  EXPECT_FLOAT_EQ(real(TypeParam(23)*v1),real(sum(&this->Array,2)->get_data_ptr()[idx]));
  EXPECT_FLOAT_EQ(imag(TypeParam(23)*v1),imag(sum(&this->Array,2)->get_data_ptr()[idx]));

  fill(&this->Array,v1);
  EXPECT_FLOAT_EQ(real(TypeParam(19)*v1),real(sum(&this->Array,3)->get_data_ptr()[idx]));
  EXPECT_FLOAT_EQ(imag(TypeParam(19)*v1),imag(sum(&this->Array,3)->get_data_ptr()[idx]));
}

TYPED_TEST_CASE(hoNDArray_utils_SpeedTestCplx, cplxImplementations);

template<class T> void
permuteTest(hoNDArray<T> *in, hoNDArray<T> *out, std::vector<size_t> *dim_order, int shift_mode = 0)
{
    if (in == 0x0 || out == 0x0 || dim_order == 0x0) {
        throw std::runtime_error("permute(): invalid pointer provided");;
    }

    if (in == out){
        throw std::runtime_error("permute(): in-place permutation not supported");;
    }

    // Check ordering array
    if (dim_order->size() > in->get_number_of_dimensions()) {
        throw std::runtime_error("hoNDArray::permute - Invalid length of dimension ordering array");;
    }

    std::vector<size_t> dim_count(in->get_number_of_dimensions(), 0);
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
        if (dim_order_int[i] == i){
            stride *= in->get_size(i);
            num_dim_memcpy = i;
        }
        else{
            break;
        }
    }

    if (stride == 1) {

        // check whether permute can be done with multi-threading
        size_t num_permute_task = 1;
        size_t num_dim_task = 0;
        for (long long i = dim_order_int.size() - 1; i >= 0; i--) {
            if (dim_order_int[i] == i){
                num_permute_task *= in->get_size(i);
                num_dim_task = i;
            }
            else{
                break;
            }
        }

        // if ((num_permute_task==1) || (num_dim_task==0)) {
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
        std::vector<size_t> dim_permute(nDim - num_dim_memcpy - 1);
        for (size_t i = num_dim_memcpy + 1; i < dim_order_int.size(); i++) {
            dim_permute[i - num_dim_memcpy - 1] = in->get_size(i);
        }

        long long n;

#ifdef USE_OMP
#pragma omp parallel default(none) shared(stride, num_dim_memcpy, num_memcpy, dim_permute, dim_order_int, nDim, in, out, o) private(n)
#endif // USE_OMP
        {
            hoNDArray<T> permuteArray(dim_permute, in->begin(), false);

            // starting index for in and out array for every permute memcpy operation
            std::vector<size_t> ind_permute_in(dim_permute.size(), 0), ind_in(nDim, 0), ind_out(nDim, 0);

#ifdef USE_OMP
#pragma omp for
#endif // USE_OMP
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
}

TYPED_TEST(hoNDArray_utils_SpeedTestCplx, permuteSpeedTest){

    fill(&this->Array, TypeParam(1, 1));

    std::vector<size_t> order;
    order.push_back(0); order.push_back(1); order.push_back(2); order.push_back(3);

    order.clear();
    order.push_back(1); order.push_back(0); order.push_back(2); order.push_back(3);

    GadgetronTimer timer(false);

    this->Array_out.create(this->Array.get_size(1), this->Array.get_size(0), this->Array.get_size(2), this->Array.get_size(3));

    timer.start("permute using openmp ... ");
    permute(&this->Array, &this->Array_out, &order);
    timer.stop();

    timer.start("permute not using openmp ... ");
    permuteTest(&this->Array, &this->Array_out, &order);
    timer.stop();
}

