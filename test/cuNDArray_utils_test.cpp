#include "cuNDArray_utils.h"
#include "cuNDArray_reductions.h"
#include "cuNDArray_blas.h"
#include "cuNDArray_elemwise.h"
#include "complext.h"

#include <gtest/gtest.h>
#include <complex>
#include <vector>
#include "vector_td_utilities.h"

using namespace Gadgetron;
using testing::Types;

template <typename T> class cuNDArray_utils_TestReal : public ::testing::Test {
protected:
  virtual void SetUp() {
    size_t vdims[] = {37, 49, 23, 19}; //Using prime numbers for setup because they are messy
    dims = std::vector<size_t>(vdims,vdims+sizeof(vdims)/sizeof(size_t));
    Array = cuNDArray<T>(&dims);
    Array2 = cuNDArray<T>(&dims);
  }
  std::vector<size_t> dims;
  cuNDArray<T> Array;
  cuNDArray<T> Array2;
};

template <typename T> class cuNDArray_utils_TestCplx : public ::testing::Test {
protected:
  virtual void SetUp() {
    size_t vdims[] = {37, 49, 23, 19}; //Using prime numbers for setup because they are messy
    dims = std::vector<size_t>(vdims,vdims+sizeof(vdims)/sizeof(size_t));
    Array = cuNDArray<T>(&dims);
    Array2 = cuNDArray<T>(&dims);
  }
  std::vector<size_t> dims;
  cuNDArray<T> Array;
  cuNDArray<T> Array2;
};

typedef Types<float, double> realImplementations;
typedef Types</*std::complex<float>, std::complex<double>,*/ float_complext, double_complext> cplxImplementations;

TYPED_TEST_SUITE(cuNDArray_utils_TestReal, realImplementations);

TYPED_TEST(cuNDArray_utils_TestReal,permuteTest){

  fill(&this->Array,TypeParam(1));

  std::vector<size_t> order;
  order.push_back(0); order.push_back(1); order.push_back(2); order.push_back(3);
  
  TypeParam tmp(2);
  CUDA_CALL(cudaMemcpy(&this->Array.get_data_ptr()[37], &tmp, sizeof(TypeParam), cudaMemcpyHostToDevice));

  EXPECT_FLOAT_EQ(1, permute(this->Array,order).at(0));
  EXPECT_FLOAT_EQ(2, permute(this->Array,order).at(37));

  order.clear();
  order.push_back(1); order.push_back(0); order.push_back(2); order.push_back(3);

  EXPECT_FLOAT_EQ(2, permute(this->Array,order).at(1));

  order.clear();
  order.push_back(3); order.push_back(1); order.push_back(2); order.push_back(0);

  EXPECT_FLOAT_EQ(2, permute(this->Array,order).at(19));

  order.clear();
  order.push_back(2); order.push_back(0); order.push_back(1); order.push_back(3);

  EXPECT_FLOAT_EQ(2, permute(this->Array,order).at(851));
}

TYPED_TEST(cuNDArray_utils_TestReal,shiftDimTest){

  fill(&this->Array,TypeParam(1));

  TypeParam tmp(2);
  CUDA_CALL(cudaMemcpy(&this->Array.get_data_ptr()[37], &tmp, sizeof(TypeParam), cudaMemcpyHostToDevice));

  EXPECT_FLOAT_EQ(1, shift_dim(this->Array,0).at(0));
  EXPECT_FLOAT_EQ(2, shift_dim(this->Array,0).at(37));
  EXPECT_FLOAT_EQ(2, shift_dim(this->Array,1).at(1));
  EXPECT_FLOAT_EQ(2, shift_dim(this->Array,-1).at(37*19));
  EXPECT_FLOAT_EQ(2, shift_dim(this->Array,2).at(23*37*19));
  EXPECT_FLOAT_EQ(2, shift_dim(this->Array,3).at(37*19));
  EXPECT_FLOAT_EQ(2, shift_dim(this->Array,4).at(37));
}

TYPED_TEST(cuNDArray_utils_TestReal,sumTest){
  TypeParam v1 = TypeParam(12.34);
  unsigned int idx = 0;

  fill(&this->Array,v1);
  EXPECT_FLOAT_EQ(49*v1,sum(&this->Array,1)->at(idx));

  fill(&this->Array,v1);
  EXPECT_FLOAT_EQ(23*v1,sum(&this->Array,2)->at(idx));

  fill(&this->Array,v1);
  EXPECT_FLOAT_EQ(19*v1,sum(&this->Array,3)->at(idx));
}


TYPED_TEST(cuNDArray_utils_TestReal,meanTest){
  TypeParam v1 = TypeParam(12.34);
  unsigned int idx = 0;

  fill(&this->Array,v1);
  EXPECT_NEAR(v1,mean(&this->Array), 0.001);

}
TYPED_TEST_SUITE(cuNDArray_utils_TestCplx, cplxImplementations);



TYPED_TEST(cuNDArray_utils_TestCplx,meanTest){
  TypeParam v1 = TypeParam(12.34);
  unsigned int idx = 0;

  fill(&this->Array,v1);
  EXPECT_FLOAT_EQ(real(v1),real(mean(&this->Array)));

}

TYPED_TEST(cuNDArray_utils_TestCplx,permuteTest){
  
  fill(&this->Array,TypeParam(1,1));

  std::vector<size_t> order;
  order.push_back(0); order.push_back(1); order.push_back(2); order.push_back(3);
  
  TypeParam tmp(2,3);
  CUDA_CALL(cudaMemcpy(&this->Array.get_data_ptr()[37], &tmp, sizeof(TypeParam), cudaMemcpyHostToDevice));

  EXPECT_FLOAT_EQ(1, real(permute(this->Array,order).at(0)));
  EXPECT_FLOAT_EQ(1, imag(permute(this->Array,order).at(0)));

  EXPECT_FLOAT_EQ(2, real(permute(this->Array,order).at(37)));
  EXPECT_FLOAT_EQ(3, imag(permute(this->Array,order).at(37)));

  order.clear();
  order.push_back(1); order.push_back(0); order.push_back(2); order.push_back(3);

  EXPECT_FLOAT_EQ(2, real(permute(this->Array,order).at(1)));
  EXPECT_FLOAT_EQ(3, imag(permute(this->Array,order).at(1)));

  order.clear();
  order.push_back(3); order.push_back(1); order.push_back(2); order.push_back(0);

  EXPECT_FLOAT_EQ(2, real(permute(this->Array,order).at(19)));
  EXPECT_FLOAT_EQ(3, imag(permute(this->Array,order).at(19)));

  order.clear();
  order.push_back(2); order.push_back(0); order.push_back(1); order.push_back(3);

  EXPECT_FLOAT_EQ(2, real(permute(this->Array,order).at(851)));
  EXPECT_FLOAT_EQ(3, imag(permute(this->Array,order).at(851)));
}

TYPED_TEST(cuNDArray_utils_TestCplx,shiftDimTest){

  fill(&this->Array,TypeParam(1,1));

  TypeParam tmp(2,3);
  CUDA_CALL(cudaMemcpy(&this->Array.get_data_ptr()[37], &tmp, sizeof(TypeParam), cudaMemcpyHostToDevice));

  EXPECT_FLOAT_EQ(1, real(shift_dim(this->Array,0).at(0)));
  EXPECT_FLOAT_EQ(1, imag(shift_dim(this->Array,0).at(0)));

  EXPECT_FLOAT_EQ(2, real(shift_dim(this->Array,0).at(37)));
  EXPECT_FLOAT_EQ(3, imag(shift_dim(this->Array,0).at(37)));

  EXPECT_FLOAT_EQ(2, real(shift_dim(this->Array,1).at(1)));
  EXPECT_FLOAT_EQ(3, imag(shift_dim(this->Array,1).at(1)));

  EXPECT_FLOAT_EQ(2, real(shift_dim(this->Array,-1).at(37*19)));
  EXPECT_FLOAT_EQ(3, imag(shift_dim(this->Array,-1).at(37*19)));

  EXPECT_FLOAT_EQ(2, real(shift_dim(this->Array,2).at(23*37*19)));
  EXPECT_FLOAT_EQ(3, imag(shift_dim(this->Array,2).at(23*37*19)));

  EXPECT_FLOAT_EQ(2, real(shift_dim(this->Array,3).at(37*19)));
  EXPECT_FLOAT_EQ(3, imag(shift_dim(this->Array,3).at(37*19)));

  EXPECT_FLOAT_EQ(2, real(shift_dim(this->Array,4).at(37)));
  EXPECT_FLOAT_EQ(3, imag(shift_dim(this->Array,4).at(37)));
}

TYPED_TEST(cuNDArray_utils_TestCplx,sumTest){
  TypeParam v1 = TypeParam(12.34, 56.78);
  unsigned int idx = 0;

  fill(&this->Array,v1);
  EXPECT_FLOAT_EQ(real(TypeParam(49)*v1),real(sum(&this->Array,1)->at(idx)));
  EXPECT_FLOAT_EQ(imag(TypeParam(49)*v1),imag(sum(&this->Array,1)->at(idx)));

  fill(&this->Array,v1);
  EXPECT_FLOAT_EQ(real(TypeParam(23)*v1),real(sum(&this->Array,2)->at(idx)));
  EXPECT_FLOAT_EQ(imag(TypeParam(23)*v1),imag(sum(&this->Array,2)->at(idx)));

  fill(&this->Array,v1);
  EXPECT_FLOAT_EQ(real(TypeParam(19)*v1),real(sum(&this->Array,3)->at(idx)));
  EXPECT_FLOAT_EQ(imag(TypeParam(19)*v1),imag(sum(&this->Array,3)->at(idx)));
}

TYPED_TEST(cuNDArray_utils_TestCplx,padTest){
  TypeParam v1 = TypeParam(12.34, 56.78);
  unsigned int idx = 0;

  fill(&this->Array,v1);

  vector_td<size_t,4> size = from_std_vector<size_t,4>(this->dims);
  size *= 2;

  auto out = pad<TypeParam,4>(size,this->Array);

  double scale = std::pow(2.0,4);
  EXPECT_EQ(out.get_number_of_elements(),this->Array.get_number_of_elements()*scale);
  EXPECT_NEAR(real(mean(&out))*scale,real(mean(&this->Array)), 0.001);
  EXPECT_NEAR(imag(mean(&out))*scale,imag(mean(&this->Array)), 0.001);
}


TEST(padTest,largeSize){
// So, this test is mainly here because pad apparently fails for large sized arrays.
	size_t vdims[] = {192,192,50};
	std::vector<size_t> dims(vdims,vdims+sizeof(vdims)/sizeof(size_t));
	size_t vdims2[] = {256,256,256};
	std::vector<size_t> dims2(vdims2,vdims2+sizeof(vdims2)/sizeof(size_t));

	cuNDArray<float_complext> in(&dims);
	fill(&in,float_complext(1));
	cuNDArray<float_complext> out(&dims2);

	pad<float_complext,3>(in,out);

	EXPECT_FLOAT_EQ(nrm2(&in),nrm2(&out));

}
