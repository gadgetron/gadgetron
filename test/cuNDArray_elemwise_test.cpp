#include "cuNDArray_elemwise.h"
#include "complext.h"

#include <gtest/gtest.h>
#include <complex>
#include <vector>

using namespace Gadgetron;
using testing::Types;

template <typename T> class cuNDArray_elemwise_TestReal : public ::testing::Test {
protected:
  virtual void SetUp() {
    size_t vdims[] = {37, 49, 23, 19}; //Using prime numbers for setup because they are messy
    dims = std::vector<size_t>(vdims,vdims+sizeof(vdims)/sizeof(size_t));
    Array = cuNDArray<T>(dims);
    Array2 = cuNDArray<T>(dims);
  }
  std::vector<size_t> dims;
  cuNDArray<T> Array;
  cuNDArray<T> Array2;
};

template <typename T> class cuNDArray_elemwise_TestCplx : public ::testing::Test {
protected:
  virtual void SetUp() {
    size_t vdims[] = {37, 49, 23, 19}; //Using prime numbers for setup because they are messy
    dims = std::vector<size_t>(vdims,vdims+sizeof(vdims)/sizeof(size_t));
    Array = cuNDArray<T>(dims);
    Array2 = cuNDArray<T>(dims);
  }
  std::vector<size_t> dims;
  cuNDArray<T> Array;
  cuNDArray<T> Array2;
};

template <typename T> class cuNDArray_elemwise_TestCplx2 : public ::testing::Test {
protected:
  virtual void SetUp() {
    size_t vdims[] = {37, 49, 23, 19}; //Using prime numbers for setup because they are messy
    dims = std::vector<size_t>(vdims,vdims+sizeof(vdims)/sizeof(size_t));
    Array = cuNDArray<T>(dims);
    Array2 = cuNDArray<T>(dims);
  }
  std::vector<size_t> dims;
  cuNDArray<T> Array;
  cuNDArray<T> Array2;
};

template <typename T> class cuNDArray_elemwise_TestCplx3 : public ::testing::Test {
protected:
  virtual void SetUp() {
    size_t vdims[] = {37, 49, 23, 19}; //Using prime numbers for setup because they are messy
    dims = std::vector<size_t>(vdims,vdims+sizeof(vdims)/sizeof(size_t));
    Array = cuNDArray<T>(dims);
    Array2 = cuNDArray<T>(dims);
  }
  std::vector<size_t> dims;
  cuNDArray<T> Array;
  cuNDArray<T> Array2;
};

template <typename T> class cuNDArray_elemwise_TestCplx4 : public ::testing::Test {
protected:
  virtual void SetUp() {
    size_t vdims[] = {37, 49, 23, 19}; //Using prime numbers for setup because they are messy
    dims = std::vector<size_t>(vdims,vdims+sizeof(vdims)/sizeof(size_t));
    Array = cuNDArray<T>(dims);
    Array2 = cuNDArray<typename realType<T>::Type>(dims);
  }
  std::vector<size_t> dims;
  cuNDArray<T> Array;
  cuNDArray<typename realType<T>::Type> Array2;
};

typedef Types<float, double> realImplementations;
typedef Types</*std::complex<float>, std::complex<double>,*/ float_complext> cplxImplementations;

TYPED_TEST_SUITE(cuNDArray_elemwise_TestReal, realImplementations);

TYPED_TEST(cuNDArray_elemwise_TestReal,fillTest){
  fill(&this->Array,TypeParam(1.1));
  EXPECT_FLOAT_EQ(1.1,TypeParam(this->Array[5]));
  fill(&this->Array,TypeParam(27.45));
  EXPECT_FLOAT_EQ(27.45,TypeParam(this->Array[3242]));
}

TYPED_TEST(cuNDArray_elemwise_TestReal,clearTest){
  fill(&this->Array,TypeParam(1));
  EXPECT_FLOAT_EQ(1,TypeParam(this->Array[5324]));
  clear(&this->Array);
  EXPECT_FLOAT_EQ(0,TypeParam(this->Array[5324]));
}

TYPED_TEST(cuNDArray_elemwise_TestReal,absTest){
  fill(&this->Array,TypeParam(-5.5));
  EXPECT_FLOAT_EQ(TypeParam(-5.5),TypeParam(this->Array[13]));
  EXPECT_FLOAT_EQ(TypeParam(5.5),TypeParam(abs(&this->Array)->at(13)));
  fill(&this->Array,TypeParam(-1.3));
  EXPECT_FLOAT_EQ(TypeParam(-1.3),TypeParam(this->Array[2454]));
  abs_inplace(&this->Array);
  EXPECT_FLOAT_EQ(TypeParam(1.3),TypeParam(this->Array[2454]));
}

TYPED_TEST(cuNDArray_elemwise_TestReal,absSquareTest){
  fill(&this->Array,TypeParam(-5.5));
  EXPECT_FLOAT_EQ(TypeParam(-5.5),TypeParam(this->Array[13]));
  EXPECT_FLOAT_EQ(TypeParam(-5.5*-5.5),TypeParam(abs_square(&this->Array)->at(13)));
}

TYPED_TEST(cuNDArray_elemwise_TestReal,sqrtTest){
  fill(&this->Array,TypeParam(17.9));
  EXPECT_FLOAT_EQ(std::sqrt(TypeParam(17.9)),TypeParam(sqrt(&this->Array)->at(23433)));
  fill(&this->Array,TypeParam(3.14));
  sqrt_inplace(&this->Array);
  EXPECT_FLOAT_EQ(std::sqrt(TypeParam(3.14)),TypeParam(this->Array[32343]));
}

TYPED_TEST(cuNDArray_elemwise_TestReal,squareTest){
  fill(&this->Array,TypeParam(1.7));
  EXPECT_FLOAT_EQ(TypeParam(1.7)*TypeParam(1.7),TypeParam(square(&this->Array)->at(22542)));
  fill(&this->Array,TypeParam(31.4));
  square_inplace(&this->Array);
  EXPECT_FLOAT_EQ(TypeParam(31.4)*TypeParam(31.4),TypeParam(this->Array[652252]));
}

TYPED_TEST(cuNDArray_elemwise_TestReal,reciprocalTest){
  fill(&this->Array,TypeParam(11.7));
  EXPECT_FLOAT_EQ(TypeParam(1)/TypeParam(11.7),TypeParam(reciprocal(&this->Array)->at(45452)));
  fill(&this->Array,TypeParam(314.114));
  reciprocal_inplace(&this->Array);
  EXPECT_FLOAT_EQ(TypeParam(1)/TypeParam(314.114),TypeParam(this->Array[43432]));
}

TYPED_TEST(cuNDArray_elemwise_TestReal,reciprocal_sqrtTest){
  fill(&this->Array,TypeParam(1.9));
  EXPECT_FLOAT_EQ(TypeParam(1)/std::sqrt(TypeParam(1.9)),TypeParam(reciprocal_sqrt(&this->Array)->at(12345)));
  fill(&this->Array,TypeParam(1.14));
  reciprocal_sqrt_inplace(&this->Array);
  EXPECT_FLOAT_EQ(TypeParam(1)/std::sqrt(TypeParam(1.14)),TypeParam(this->Array[0]));
}

TYPED_TEST(cuNDArray_elemwise_TestReal,sgnTest){
  fill(&this->Array,TypeParam(-5.7));
  TypeParam tmp(101.1);
  TypeParam tmp2(0.0);
  CUDA_CALL(cudaMemcpy(&this->Array.get_data_ptr()[91], &tmp, sizeof(TypeParam), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(&this->Array.get_data_ptr()[19100], &tmp2, sizeof(TypeParam), cudaMemcpyHostToDevice));
  EXPECT_FLOAT_EQ(TypeParam(-1),TypeParam(sgn(&this->Array)->at(28)));
  EXPECT_FLOAT_EQ(TypeParam(1),TypeParam(sgn(&this->Array)->at(91)));
  EXPECT_FLOAT_EQ(TypeParam(0),TypeParam(sgn(&this->Array)->at(19100)));
  fill(&this->Array,TypeParam(-5.7));
  CUDA_CALL(cudaMemcpy(&this->Array.get_data_ptr()[9100], &tmp, sizeof(TypeParam), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(&this->Array.get_data_ptr()[19100], &tmp2, sizeof(TypeParam), cudaMemcpyHostToDevice));
  sgn_inplace(&this->Array);
  EXPECT_FLOAT_EQ(TypeParam(-1),TypeParam(this->Array[2800]));
  EXPECT_FLOAT_EQ(TypeParam(1),TypeParam(this->Array[9100]));
  EXPECT_FLOAT_EQ(TypeParam(0),TypeParam(this->Array[19100]));
}

TYPED_TEST(cuNDArray_elemwise_TestReal,clampTest){
  fill(&this->Array,TypeParam(-5.7));
  TypeParam tmp(101.3);
  CUDA_CALL(cudaMemcpy(&this->Array.get_data_ptr()[354222], &tmp, sizeof(TypeParam), cudaMemcpyHostToDevice));
  clamp(&this->Array,TypeParam(4.9),TypeParam(100.0));
  EXPECT_FLOAT_EQ(TypeParam(4.9),this->Array[3435]);
  EXPECT_FLOAT_EQ(TypeParam(100.0),this->Array[354222]);
}

TYPED_TEST(cuNDArray_elemwise_TestReal,clamp_minTest){
  fill(&this->Array,TypeParam(-5.7));
  TypeParam tmp(-101.3);
  CUDA_CALL(cudaMemcpy(&this->Array.get_data_ptr()[91], &tmp, sizeof(TypeParam), cudaMemcpyHostToDevice));
  clamp_min(&this->Array,TypeParam(-10.6));
  EXPECT_FLOAT_EQ(TypeParam(-5.7),this->Array[28]);
  EXPECT_FLOAT_EQ(TypeParam(-10.6),this->Array[91]);
}

TYPED_TEST(cuNDArray_elemwise_TestReal,clamp_maxTest){
  fill(&this->Array,TypeParam(5.7));
  TypeParam tmp(101.3);
  CUDA_CALL(cudaMemcpy(&this->Array.get_data_ptr()[91], &tmp, sizeof(TypeParam), cudaMemcpyHostToDevice));
  clamp_max(&this->Array,TypeParam(10.6));
  EXPECT_FLOAT_EQ(TypeParam(5.7),this->Array[28]);
  EXPECT_FLOAT_EQ(TypeParam(10.6),this->Array[91]);
}

TYPED_TEST(cuNDArray_elemwise_TestReal,normalizeTest){
  fill(&this->Array,TypeParam(50));
  TypeParam tmp(-200);
  CUDA_CALL(cudaMemcpy(&this->Array.get_data_ptr()[23], &tmp, sizeof(TypeParam), cudaMemcpyHostToDevice));
  normalize(&this->Array,110);
  EXPECT_FLOAT_EQ(TypeParam(50)*TypeParam(110)/abs(TypeParam(-200)),this->Array[12345]);
}

TYPED_TEST(cuNDArray_elemwise_TestReal,shrink1Test){
  fill(&this->Array,TypeParam(1.2));
  shrink1(&this->Array,0.75);
  EXPECT_FLOAT_EQ(TypeParam(1.2)/abs(TypeParam(1.2))*std::max(abs(TypeParam(1.2))-0.75,0.0),this->Array[125]);
  fill(&this->Array,TypeParam(1));
  shrink1(&this->Array,2.0);
  EXPECT_FLOAT_EQ(0.0,this->Array[125]);
}

TYPED_TEST(cuNDArray_elemwise_TestReal,shrinkdTest){
  fill(&this->Array,TypeParam(1.2));
  fill(&this->Array2,TypeParam(4.0));
  shrinkd(&this->Array,&this->Array2,1.0);
  EXPECT_FLOAT_EQ(TypeParam(1.2)/TypeParam(4.0)*std::max(4.0-1.0,0.0),this->Array[125]);
  shrinkd(&this->Array,&this->Array2,8.0);
  EXPECT_FLOAT_EQ(0.0,this->Array[125]);
}

TYPED_TEST(cuNDArray_elemwise_TestReal,realTest){
  fill(&this->Array,TypeParam(1.2));
  EXPECT_FLOAT_EQ(TypeParam(1.2),real(&this->Array)->at(125));
}

TYPED_TEST(cuNDArray_elemwise_TestReal,imagTest){
  fill(&this->Array,TypeParam(1.2));
  EXPECT_FLOAT_EQ(TypeParam(0.0),imag(&this->Array)->at(125));
}

TYPED_TEST(cuNDArray_elemwise_TestReal,conjTest){
  fill(&this->Array,TypeParam(1.2));
  EXPECT_FLOAT_EQ(TypeParam(1.2),real(&this->Array)->at(125));
}

TYPED_TEST_SUITE(cuNDArray_elemwise_TestCplx, cplxImplementations);

TYPED_TEST(cuNDArray_elemwise_TestCplx,fillTest){
  fill(&this->Array,TypeParam(1.1,2.2));
  EXPECT_FLOAT_EQ(1.1,real(TypeParam(this->Array[52323])));
  EXPECT_FLOAT_EQ(2.2,imag(TypeParam(this->Array[52323])));
}

TYPED_TEST(cuNDArray_elemwise_TestCplx,clearTest){
  fill(&this->Array,TypeParam(1,1));
  clear(&this->Array);
  EXPECT_FLOAT_EQ(0,real(TypeParam(this->Array[325])));
  EXPECT_FLOAT_EQ(0,imag(TypeParam(this->Array[325])));
}

TYPED_TEST(cuNDArray_elemwise_TestCplx,absTest){
  fill(&this->Array,TypeParam(-5.5,7.7));
  EXPECT_FLOAT_EQ(std::sqrt(5.5*5.5+7.7*7.7),abs(&this->Array)->at(32113));
}

TYPED_TEST(cuNDArray_elemwise_TestCplx,absSquareTest){
  fill(&this->Array,TypeParam(-5.5,7.7));
  EXPECT_FLOAT_EQ(5.5*5.5+7.7*7.7,abs_square(&this->Array)->at(32113));
}

TYPED_TEST(cuNDArray_elemwise_TestCplx,sqrtTest){
  fill(&this->Array,TypeParam(17.9,3.5));
  EXPECT_NEAR(real(sqrt(TypeParam(17.9,3.5))),real(sqrt(&this->Array)->at(2131)),0.00001);
  EXPECT_NEAR(imag(sqrt(TypeParam(17.9,3.5))),imag(sqrt(&this->Array)->at(2131)),0.00001);
  fill(&this->Array,TypeParam(3.14,4.13));
  sqrt_inplace(&this->Array);
  EXPECT_NEAR(real(sqrt(TypeParam(3.14,4.13))),real(this->Array[120000]),0.00001);
  EXPECT_NEAR(imag(sqrt(TypeParam(3.14,4.13))),imag(this->Array[120000]),0.00001);
}

TYPED_TEST(cuNDArray_elemwise_TestCplx,squareTest){
  fill(&this->Array,TypeParam(1.7,7.1));
  EXPECT_FLOAT_EQ(real(TypeParam(1.7,7.1)*TypeParam(1.7,7.1)),real(square(&this->Array)->at(22123)));
  EXPECT_FLOAT_EQ(imag(TypeParam(1.7,7.1)*TypeParam(1.7,7.1)),imag(square(&this->Array)->at(22123)));
  fill(&this->Array,TypeParam(31.4,4.31));
  square_inplace(&this->Array);
  EXPECT_FLOAT_EQ(real(TypeParam(31.4,4.31)*TypeParam(31.4,4.31)),real(this->Array[51234]));
  EXPECT_FLOAT_EQ(imag(TypeParam(31.4,4.31)*TypeParam(31.4,4.31)),imag(this->Array[51234]));
}

TYPED_TEST(cuNDArray_elemwise_TestCplx,reciprocalTest){
  fill(&this->Array,TypeParam(1.9,2.7));
  EXPECT_FLOAT_EQ(real(TypeParam(1,0)/TypeParam(1.9,2.7)),real(reciprocal(&this->Array)->at(11232)));
  EXPECT_FLOAT_EQ(imag(TypeParam(1,0)/TypeParam(1.9,2.7)),imag(reciprocal(&this->Array)->at(11232)));
  fill(&this->Array,TypeParam(1.14,4.32));
  reciprocal_inplace(&this->Array);
  EXPECT_FLOAT_EQ(real(TypeParam(1,0)/TypeParam(1.14,4.32)),real(this->Array[10]));
  EXPECT_FLOAT_EQ(imag(TypeParam(1,0)/TypeParam(1.14,4.32)),imag(this->Array[10]));
}

TYPED_TEST(cuNDArray_elemwise_TestCplx,reciprocal_sqrtTest){
  fill(&this->Array,TypeParam(1.9,2.7));
  EXPECT_FLOAT_EQ(real(TypeParam(1,0)/sqrt(TypeParam(1.9,2.7))),real(reciprocal_sqrt(&this->Array)->at(12543)));
  EXPECT_FLOAT_EQ(imag(TypeParam(1,0)/sqrt(TypeParam(1.9,2.7))),imag(reciprocal_sqrt(&this->Array)->at(12543)));
  fill(&this->Array,TypeParam(1.14,4.32));
  reciprocal_sqrt_inplace(&this->Array);
  EXPECT_FLOAT_EQ(real(TypeParam(1,0)/sqrt(TypeParam(1.14,4.32))),real(this->Array[10000]));
  EXPECT_FLOAT_EQ(imag(TypeParam(1,0)/sqrt(TypeParam(1.14,4.32))),imag(this->Array[10000]));
}

TYPED_TEST(cuNDArray_elemwise_TestCplx,realImagTest){
  fill(&this->Array,TypeParam(3.4,4.2));
  EXPECT_FLOAT_EQ(3.4,real(&this->Array)->at(33425));
  EXPECT_FLOAT_EQ(4.2,imag(&this->Array)->at(45));
}

TYPED_TEST(cuNDArray_elemwise_TestCplx,conjTest){
  fill(&this->Array,TypeParam(3.4,4.2));
  EXPECT_FLOAT_EQ(3.4,real(conj(&this->Array)->at(33425)));
  EXPECT_FLOAT_EQ(-4.2,imag(conj(&this->Array)->at(45)));
}

TYPED_TEST(cuNDArray_elemwise_TestCplx,normalizeTest){
  fill(&this->Array,TypeParam(50,50));
  TypeParam tmp(-200,-200);
  CUDA_CALL(cudaMemcpy(&this->Array.get_data_ptr()[23], &tmp, sizeof(TypeParam), cudaMemcpyHostToDevice));
  normalize(&this->Array,110);
  EXPECT_FLOAT_EQ(real(TypeParam(50,50)*real(TypeParam(110,110))/abs(TypeParam(-200,-200))),real(&this->Array)->at(12345));
  EXPECT_FLOAT_EQ(imag(TypeParam(50,50)*real(TypeParam(110,110))/abs(TypeParam(-200,-200))),imag(&this->Array)->at(12345));
}

TYPED_TEST(cuNDArray_elemwise_TestCplx,clampTest){
  fill(&this->Array,TypeParam(-5.7, -4.6));
  TypeParam tmp(101.3,203.4);
  CUDA_CALL(cudaMemcpy(&this->Array.get_data_ptr()[354222], &tmp, sizeof(TypeParam), cudaMemcpyHostToDevice));  
  clamp(&this->Array,real(TypeParam(4.9,0)),real(TypeParam(100.0,0)));
  EXPECT_FLOAT_EQ(real(TypeParam(4.9,0)),real(&this->Array)->at(3435));
  EXPECT_FLOAT_EQ(real(TypeParam(100.0,0)),real(&this->Array)->at(354222));
  EXPECT_FLOAT_EQ(imag(TypeParam(4.9,0)),imag(&this->Array)->at(3435));
  EXPECT_FLOAT_EQ(imag(TypeParam(100.0,0)),imag(&this->Array)->at(354222));
}

TYPED_TEST(cuNDArray_elemwise_TestCplx,clamp_minTest){
  fill(&this->Array,TypeParam(-5.7, -4.6));
  TypeParam tmp(-101.3,-203.4);
  CUDA_CALL(cudaMemcpy(&this->Array.get_data_ptr()[91], &tmp, sizeof(TypeParam), cudaMemcpyHostToDevice));  
  clamp_min(&this->Array, real(TypeParam(-10.6,0)));
  EXPECT_FLOAT_EQ(real(TypeParam(-5.7,0)),real(&this->Array)->at(28));
  EXPECT_FLOAT_EQ(real(TypeParam(-10.6,0)),real(&this->Array)->at(91));
  EXPECT_FLOAT_EQ(imag(TypeParam(-5.7,0)),imag(&this->Array)->at(28));
  EXPECT_FLOAT_EQ(imag(TypeParam(-10.6,0)),imag(&this->Array)->at(91));
}

TYPED_TEST(cuNDArray_elemwise_TestCplx,clamp_maxTest){
  fill(&this->Array,TypeParam(5.7, 4.6));
  TypeParam tmp(101.3,203.4);
  CUDA_CALL(cudaMemcpy(&this->Array.get_data_ptr()[91], &tmp, sizeof(TypeParam), cudaMemcpyHostToDevice));  
  clamp_max(&this->Array,real(TypeParam(10.6,0)));
  EXPECT_FLOAT_EQ(real(TypeParam(5.7,0)),real(&this->Array)->at(28));
  EXPECT_FLOAT_EQ(real(TypeParam(10.6,0)),real(&this->Array)->at(91));
  EXPECT_FLOAT_EQ(imag(TypeParam(5.7,0)),imag(&this->Array)->at(28));
  EXPECT_FLOAT_EQ(imag(TypeParam(10.6,0)),imag(&this->Array)->at(91));
}

TYPED_TEST(cuNDArray_elemwise_TestCplx,shrink1Test){
  fill(&this->Array,TypeParam(1.2,1.4));
  shrink1(&this->Array,0.75);
  EXPECT_FLOAT_EQ(real(TypeParam(1.2,1.4)/abs(TypeParam(1.2,1.4)))*std::max(abs(TypeParam(1.2,1.4))-0.75,0.0),real(&this->Array)->at(125));
  EXPECT_FLOAT_EQ(imag(TypeParam(1.2,1.4)/abs(TypeParam(1.2,1.4)))*std::max(abs(TypeParam(1.2,1.4))-0.75,0.0),imag(&this->Array)->at(125));
  fill(&this->Array,TypeParam(1,1));
  shrink1(&this->Array,2.0);
  EXPECT_FLOAT_EQ(0.0,real(&this->Array)->at(125));
  EXPECT_FLOAT_EQ(0.0,imag(&this->Array)->at(23125));
}

TYPED_TEST_SUITE(cuNDArray_elemwise_TestCplx4, cplxImplementations);

TYPED_TEST(cuNDArray_elemwise_TestCplx4,shrinkdTest){
  fill(&this->Array,TypeParam(1.2,1.4));
  fill(&this->Array2,real(TypeParam(4.0,4.0)));
  shrinkd(&this->Array,&this->Array2,1.0);
  EXPECT_FLOAT_EQ(real(TypeParam(1.2,1.4)/real(TypeParam(4.0,4.0)))*std::max(4.0-1.0,0.0),real(&this->Array)->at(125));
  EXPECT_FLOAT_EQ(imag(TypeParam(1.2,1.4)/imag(TypeParam(4.0,4.0)))*std::max(4.0-1.0,0.0),imag(&this->Array)->at(125));
  shrinkd(&this->Array,&this->Array2,8.0);
  EXPECT_FLOAT_EQ(0.0,real(&this->Array)->at(125));
  EXPECT_FLOAT_EQ(0.0,imag(&this->Array)->at(23125));
}

TYPED_TEST_SUITE(cuNDArray_elemwise_TestCplx3, cplxImplementations);

TYPED_TEST(cuNDArray_elemwise_TestCplx3,realToCplxTest){
  fill(&this->Array,TypeParam(3.4,4.2));
  EXPECT_FLOAT_EQ(3.4,real(real_to_complex<TypeParam>(real(&this->Array).get())->at(33425)));
  EXPECT_FLOAT_EQ(0.0,imag(real_to_complex<TypeParam>(real(&this->Array).get())->at(33425)));
}
