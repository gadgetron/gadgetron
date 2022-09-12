#include "hoNDArray_elemwise.h"
#include "complext.h"

#include <gtest/gtest.h>
#include <complex>
#include <vector>

using namespace Gadgetron;
using testing::Types;

template <typename T> class hoNDArray_elemwise_TestReal : public ::testing::Test {
protected:
  virtual void SetUp() {
    size_t vdims[] = {37, 49, 23, 19}; //Using prime numbers for setup because they are messy
    dims = std::vector<size_t>(vdims,vdims+sizeof(vdims)/sizeof(size_t));
    Array = hoNDArray<T>(dims);
    Array2 = hoNDArray<T>(dims);
  }
  std::vector<size_t> dims;
  hoNDArray<T> Array;
  hoNDArray<T> Array2;
};

template <typename T> class hoNDArray_elemwise_TestCplx : public ::testing::Test {
protected:
  virtual void SetUp() {
    size_t vdims[] = {37, 49, 23, 19}; //Using prime numbers for setup because they are messy
    dims = std::vector<size_t>(vdims,vdims+sizeof(vdims)/sizeof(size_t));
    Array = hoNDArray<T>(dims);
    Array2 = hoNDArray<T>(dims);
  }
  std::vector<size_t> dims;
  hoNDArray<T> Array;
  hoNDArray<T> Array2;
};

template <typename T> class hoNDArray_elemwise_TestCplx2 : public ::testing::Test {
protected:
  virtual void SetUp() {
    size_t vdims[] = {37, 49, 23, 19}; //Using prime numbers for setup because they are messy
    dims = std::vector<size_t>(vdims,vdims+sizeof(vdims)/sizeof(size_t));
    Array = hoNDArray<T>(dims);
    Array2 = hoNDArray<T>(dims);
  }
  std::vector<size_t> dims;
  hoNDArray<T> Array;
  hoNDArray<T> Array2;
};

template <typename T> class hoNDArray_elemwise_TestCplx3 : public ::testing::Test {
protected:
  virtual void SetUp() {
    size_t vdims[] = {37, 49, 23, 19}; //Using prime numbers for setup because they are messy
    dims = std::vector<size_t>(vdims,vdims+sizeof(vdims)/sizeof(size_t));
    Array = hoNDArray<T>(dims);
    Array2 = hoNDArray<T>(dims);
  }
  std::vector<size_t> dims;
  hoNDArray<T> Array;
  hoNDArray<T> Array2;
};

template <typename T> class hoNDArray_elemwise_TestCplx4 : public ::testing::Test {
protected:
  virtual void SetUp() {
    size_t vdims[] = {37, 49, 23, 19}; //Using prime numbers for setup because they are messy
    dims = std::vector<size_t>(vdims,vdims+sizeof(vdims)/sizeof(size_t));
    Array = hoNDArray<T>(dims);
    Array2 = hoNDArray<typename realType<T>::Type>(dims);
  }
  std::vector<size_t> dims;
  hoNDArray<T> Array;
  hoNDArray<typename realType<T>::Type> Array2;
};

typedef Types<float, double> realImplementations;
typedef Types<std::complex<float>, std::complex<double>, float_complext, double_complext> cplxImplementations;
typedef Types<std::complex<float>, std::complex<double> > stdCplxImplementations;
typedef Types<float_complext, double_complext> cplxtImplementations;

TYPED_TEST_SUITE(hoNDArray_elemwise_TestReal, realImplementations);

TYPED_TEST(hoNDArray_elemwise_TestReal,fillTest){
  fill(&this->Array,TypeParam(1.1));
  EXPECT_FLOAT_EQ(1.1,TypeParam(this->Array.get_data_ptr()[5]));
  fill(&this->Array,TypeParam(27.45));
  EXPECT_FLOAT_EQ(27.45,TypeParam(this->Array.get_data_ptr()[3242]));
}

TYPED_TEST(hoNDArray_elemwise_TestReal,clearTest){
  fill(&this->Array,TypeParam(1));
  EXPECT_FLOAT_EQ(1,TypeParam(this->Array.get_data_ptr()[5324]));
  clear(&this->Array);
  EXPECT_FLOAT_EQ(0,TypeParam(this->Array.get_data_ptr()[5324]));
}

TYPED_TEST(hoNDArray_elemwise_TestReal,absTest){
  fill(&this->Array,TypeParam(-5.5));
  EXPECT_FLOAT_EQ(TypeParam(-5.5),TypeParam(this->Array.get_data_ptr()[13]));
  EXPECT_FLOAT_EQ(TypeParam(5.5),TypeParam(abs(&this->Array)->get_data_ptr()[13]));
  fill(&this->Array,TypeParam(-1.3));
  EXPECT_FLOAT_EQ(TypeParam(-1.3),TypeParam(this->Array.get_data_ptr()[2454]));
  abs_inplace(&this->Array);
  EXPECT_FLOAT_EQ(TypeParam(1.3),TypeParam(this->Array.get_data_ptr()[2454]));
}

TYPED_TEST(hoNDArray_elemwise_TestReal,absSquareTest){
  fill(&this->Array,TypeParam(-5.5));
  EXPECT_FLOAT_EQ(TypeParam(-5.5),TypeParam(this->Array.get_data_ptr()[13]));
  EXPECT_FLOAT_EQ(TypeParam(-5.5*-5.5),TypeParam(abs_square(&this->Array)->get_data_ptr()[13]));
}

TYPED_TEST(hoNDArray_elemwise_TestReal,sqrtTest){
  fill(&this->Array,TypeParam(17.9));
  EXPECT_FLOAT_EQ(std::sqrt(TypeParam(17.9)),TypeParam(sqrt(&this->Array)->get_data_ptr()[23433]));
  fill(&this->Array,TypeParam(3.14));
  sqrt_inplace(&this->Array);
  EXPECT_FLOAT_EQ(std::sqrt(TypeParam(3.14)),TypeParam(this->Array.get_data_ptr()[32343]));
}

TYPED_TEST(hoNDArray_elemwise_TestReal,squareTest){
  fill(&this->Array,TypeParam(1.7));
  EXPECT_FLOAT_EQ(TypeParam(1.7)*TypeParam(1.7),TypeParam(square(&this->Array)->get_data_ptr()[22542]));
  fill(&this->Array,TypeParam(31.4));
  square_inplace(&this->Array);
  EXPECT_FLOAT_EQ(TypeParam(31.4)*TypeParam(31.4),TypeParam(this->Array.get_data_ptr()[652252]));
}

TYPED_TEST(hoNDArray_elemwise_TestReal,reciprocalTest){
  fill(&this->Array,TypeParam(11.7));
  EXPECT_FLOAT_EQ(TypeParam(1)/TypeParam(11.7),TypeParam(reciprocal(&this->Array)->get_data_ptr()[45452]));
  fill(&this->Array,TypeParam(314.114));
  reciprocal_inplace(&this->Array);
  EXPECT_FLOAT_EQ(TypeParam(1)/TypeParam(314.114),TypeParam(this->Array.get_data_ptr()[43432]));
}

TYPED_TEST(hoNDArray_elemwise_TestReal,reciprocal_sqrtTest){
  fill(&this->Array,TypeParam(1.9));
  EXPECT_FLOAT_EQ(TypeParam(1)/std::sqrt(TypeParam(1.9)),TypeParam(reciprocal_sqrt(&this->Array)->get_data_ptr()[12345]));
  fill(&this->Array,TypeParam(1.14));
  reciprocal_sqrt_inplace(&this->Array);
  EXPECT_FLOAT_EQ(TypeParam(1)/std::sqrt(TypeParam(1.14)),TypeParam(this->Array.get_data_ptr()[0]));
}

TYPED_TEST(hoNDArray_elemwise_TestReal,sgnTest){
  fill(&this->Array,TypeParam(-5.7));
  this->Array.get_data_ptr()[91] = TypeParam(101.1);
  this->Array.get_data_ptr()[19100] = TypeParam(0);
  EXPECT_FLOAT_EQ(TypeParam(-1),TypeParam(sgn(&this->Array)->get_data_ptr()[28]));
  EXPECT_FLOAT_EQ(TypeParam(1),TypeParam(sgn(&this->Array)->get_data_ptr()[91]));
  EXPECT_FLOAT_EQ(TypeParam(0),TypeParam(sgn(&this->Array)->get_data_ptr()[19100]));
  fill(&this->Array,TypeParam(-5.7));
  this->Array.get_data_ptr()[9100] = TypeParam(101.1);
  this->Array.get_data_ptr()[19100] = TypeParam(0);
  sgn_inplace(&this->Array);
  EXPECT_FLOAT_EQ(TypeParam(-1),TypeParam(this->Array.get_data_ptr()[2800]));
  EXPECT_FLOAT_EQ(TypeParam(1),TypeParam(this->Array.get_data_ptr()[9100]));
  EXPECT_FLOAT_EQ(TypeParam(0),TypeParam(this->Array.get_data_ptr()[19100]));
}

TYPED_TEST(hoNDArray_elemwise_TestReal,clampTest){
  fill(&this->Array,TypeParam(-5.7));
  this->Array.get_data_ptr()[354222] = TypeParam(101.3);
  clamp(&this->Array,TypeParam(4.9),TypeParam(100.0));
  EXPECT_FLOAT_EQ(TypeParam(4.9),this->Array.get_data_ptr()[3435]);
  EXPECT_FLOAT_EQ(TypeParam(100.0),this->Array.get_data_ptr()[354222]);
}

TYPED_TEST(hoNDArray_elemwise_TestReal,clamp_minTest){
  fill(&this->Array,TypeParam(-5.7));
  this->Array.get_data_ptr()[91] = TypeParam(-101.3);
  clamp_min(&this->Array,TypeParam(-10.6));
  EXPECT_FLOAT_EQ(TypeParam(-5.7),this->Array.get_data_ptr()[28]);
  EXPECT_FLOAT_EQ(TypeParam(-10.6),this->Array.get_data_ptr()[91]);
}

TYPED_TEST(hoNDArray_elemwise_TestReal,clamp_maxTest){
  fill(&this->Array,TypeParam(5.7));
  this->Array.get_data_ptr()[91] = TypeParam(101.3);
  clamp_max(&this->Array,TypeParam(10.6));
  EXPECT_FLOAT_EQ(TypeParam(5.7),this->Array.get_data_ptr()[28]);
  EXPECT_FLOAT_EQ(TypeParam(10.6),this->Array.get_data_ptr()[91]);
}

TYPED_TEST(hoNDArray_elemwise_TestReal,normalizeTest){
  fill(&this->Array,TypeParam(50));
  this->Array.get_data_ptr()[23]=TypeParam(-200);
  normalize(&this->Array,110);
  EXPECT_FLOAT_EQ(TypeParam(50)*TypeParam(110)/abs(TypeParam(-200)),this->Array.get_data_ptr()[12345]);
}

TYPED_TEST(hoNDArray_elemwise_TestReal,shrink1Test){
  fill(&this->Array,TypeParam(1.2));
  shrink1(&this->Array,0.75);
  EXPECT_FLOAT_EQ(TypeParam(1.2)/abs(TypeParam(1.2))*std::max(abs(TypeParam(1.2))-0.75,0.0),this->Array.get_data_ptr()[125]);
  fill(&this->Array,TypeParam(1));
  shrink1(&this->Array,2.0);
  EXPECT_FLOAT_EQ(0.0,this->Array.get_data_ptr()[125]);
}

TYPED_TEST(hoNDArray_elemwise_TestReal,shrinkdTest){
  fill(&this->Array,TypeParam(1.2));
  fill(&this->Array2,TypeParam(4.0));
  shrinkd(&this->Array,&this->Array2,1.0);
  EXPECT_FLOAT_EQ(TypeParam(1.2)/TypeParam(4.0)*std::max(4.0-1.0,0.0),this->Array.get_data_ptr()[125]);
  shrinkd(&this->Array,&this->Array2,8.0);
  EXPECT_FLOAT_EQ(0.0,this->Array.get_data_ptr()[125]);
}

TYPED_TEST(hoNDArray_elemwise_TestReal,realTest){
  fill(&this->Array,TypeParam(1.2));
  EXPECT_FLOAT_EQ(TypeParam(1.2),real(&this->Array)->at(125));
}

TYPED_TEST(hoNDArray_elemwise_TestReal,imagTest){
  fill(&this->Array,TypeParam(1.2));
  EXPECT_FLOAT_EQ(TypeParam(0.0),imag(&this->Array)->at(125));
}

TYPED_TEST(hoNDArray_elemwise_TestReal,conjTest){
  fill(&this->Array,TypeParam(1.2));
  EXPECT_FLOAT_EQ(TypeParam(1.2),real(&this->Array)->at(125));
  EXPECT_FLOAT_EQ(TypeParam(0.0),imag(&this->Array)->at(125));
}

TYPED_TEST_SUITE(hoNDArray_elemwise_TestCplx, cplxImplementations);

TYPED_TEST(hoNDArray_elemwise_TestCplx,fillTest){
  fill(&this->Array,TypeParam(1.1,2.2));
  EXPECT_FLOAT_EQ(1.1,real(TypeParam(this->Array.get_data_ptr()[52323])));
  EXPECT_FLOAT_EQ(2.2,imag(TypeParam(this->Array.get_data_ptr()[52323])));
}

TYPED_TEST(hoNDArray_elemwise_TestCplx,clearTest){
  fill(&this->Array,TypeParam(1,1));
  clear(&this->Array);
  EXPECT_FLOAT_EQ(0,real(TypeParam(this->Array.get_data_ptr()[325])));
  EXPECT_FLOAT_EQ(0,imag(TypeParam(this->Array.get_data_ptr()[325])));
}

TYPED_TEST(hoNDArray_elemwise_TestCplx,absTest){
  fill(&this->Array,TypeParam(-5.5,7.7));
  EXPECT_FLOAT_EQ(std::sqrt(5.5*5.5+7.7*7.7),abs(&this->Array)->get_data_ptr()[32113]);
}

TYPED_TEST(hoNDArray_elemwise_TestCplx,absSquareTest){
  fill(&this->Array,TypeParam(-5.5,7.7));
  EXPECT_FLOAT_EQ(5.5*5.5+7.7*7.7,abs_square(&this->Array)->get_data_ptr()[32113]);
}

TYPED_TEST(hoNDArray_elemwise_TestCplx,sqrtTest){
  fill(&this->Array,TypeParam(17.9,3.5));
  EXPECT_NEAR(real(sqrt(TypeParam(17.9,3.5))),real(sqrt(&this->Array)->get_data_ptr()[2131]),0.00001);
  EXPECT_NEAR(imag(sqrt(TypeParam(17.9,3.5))),imag(sqrt(&this->Array)->get_data_ptr()[2131]),0.00001);
  fill(&this->Array,TypeParam(3.14,4.13));
  sqrt_inplace(&this->Array);
  EXPECT_NEAR(real(sqrt(TypeParam(3.14,4.13))),real(this->Array.get_data_ptr()[120000]),0.00001);
  EXPECT_NEAR(imag(sqrt(TypeParam(3.14,4.13))),imag(this->Array.get_data_ptr()[120000]),0.00001);
}

TYPED_TEST(hoNDArray_elemwise_TestCplx,squareTest){
  fill(&this->Array,TypeParam(1.7,7.1));
  EXPECT_FLOAT_EQ(real(TypeParam(1.7,7.1)*TypeParam(1.7,7.1)),real(square(&this->Array)->get_data_ptr()[22123]));
  EXPECT_FLOAT_EQ(imag(TypeParam(1.7,7.1)*TypeParam(1.7,7.1)),imag(square(&this->Array)->get_data_ptr()[22123]));
  fill(&this->Array,TypeParam(31.4,4.31));
  square_inplace(&this->Array);
  EXPECT_FLOAT_EQ(real(TypeParam(31.4,4.31)*TypeParam(31.4,4.31)),real(this->Array.get_data_ptr()[51234]));
  EXPECT_FLOAT_EQ(imag(TypeParam(31.4,4.31)*TypeParam(31.4,4.31)),imag(this->Array.get_data_ptr()[51234]));
}

TYPED_TEST(hoNDArray_elemwise_TestCplx,reciprocalTest){
  fill(&this->Array,TypeParam(1.9,2.7));
  EXPECT_FLOAT_EQ(real(TypeParam(1,0)/TypeParam(1.9,2.7)),real(reciprocal(&this->Array)->get_data_ptr()[11232]));
  EXPECT_FLOAT_EQ(imag(TypeParam(1,0)/TypeParam(1.9,2.7)),imag(reciprocal(&this->Array)->get_data_ptr()[11232]));
  fill(&this->Array,TypeParam(1.14,4.32));
  reciprocal_inplace(&this->Array);
  EXPECT_FLOAT_EQ(real(TypeParam(1,0)/TypeParam(1.14,4.32)),real(this->Array.get_data_ptr()[10]));
  EXPECT_FLOAT_EQ(imag(TypeParam(1,0)/TypeParam(1.14,4.32)),imag(this->Array.get_data_ptr()[10]));
}

TYPED_TEST(hoNDArray_elemwise_TestCplx,reciprocal_sqrtTest){
  fill(&this->Array,TypeParam(1.9,2.7));
  EXPECT_FLOAT_EQ(real(TypeParam(1,0)/sqrt(TypeParam(1.9,2.7))),real(reciprocal_sqrt(&this->Array)->get_data_ptr()[12543]));
  EXPECT_FLOAT_EQ(imag(TypeParam(1,0)/sqrt(TypeParam(1.9,2.7))),imag(reciprocal_sqrt(&this->Array)->get_data_ptr()[12543]));
  fill(&this->Array,TypeParam(1.14,4.32));
  reciprocal_sqrt_inplace(&this->Array);
  EXPECT_FLOAT_EQ(real(TypeParam(1,0)/sqrt(TypeParam(1.14,4.32))),real(this->Array.get_data_ptr()[10000]));
  EXPECT_FLOAT_EQ(imag(TypeParam(1,0)/sqrt(TypeParam(1.14,4.32))),imag(this->Array.get_data_ptr()[10000]));
}

TYPED_TEST(hoNDArray_elemwise_TestCplx,realImagTest){
  fill(&this->Array,TypeParam(3.4,4.2));
  EXPECT_FLOAT_EQ(3.4,real(&this->Array)->get_data_ptr()[33425]);
  EXPECT_NEAR(4.2,imag(&this->Array)->get_data_ptr()[45], 0.000001);
}

TYPED_TEST(hoNDArray_elemwise_TestCplx,conjTest){
  fill(&this->Array,TypeParam(3.4,4.2));
  EXPECT_FLOAT_EQ(3.4,real(conj(&this->Array)->at(33425)));
  EXPECT_NEAR(-4.2,imag(conj(&this->Array)->at(45)), 0.000001);
}

TYPED_TEST(hoNDArray_elemwise_TestCplx,normalizeTest){
  fill(&this->Array,TypeParam(50,50));
  this->Array.get_data_ptr()[23]=TypeParam(-200,-200);
  normalize(&this->Array,110);
  EXPECT_FLOAT_EQ(real(TypeParam(50,50)*real(TypeParam(110,110))/abs(TypeParam(-200,-200))),real(&this->Array)->get_data_ptr()[12345]);
  EXPECT_FLOAT_EQ(imag(TypeParam(50,50)*real(TypeParam(110,110))/abs(TypeParam(-200,-200))),imag(&this->Array)->get_data_ptr()[12345]);
}

TYPED_TEST(hoNDArray_elemwise_TestCplx,clampTest){
  fill(&this->Array,TypeParam(-5.7, -4.6));
  this->Array.get_data_ptr()[354222] = TypeParam(101.3,203.4);
  clamp(&this->Array,real(TypeParam(4.9,0)),real(TypeParam(100.0,0)));
  EXPECT_FLOAT_EQ(real(TypeParam(4.9,0)),real(&this->Array)->get_data_ptr()[3435]);
  EXPECT_FLOAT_EQ(real(TypeParam(100.0,0)),real(&this->Array)->get_data_ptr()[354222]);
  EXPECT_FLOAT_EQ(imag(TypeParam(4.9,0)),imag(&this->Array)->get_data_ptr()[3435]);
  EXPECT_FLOAT_EQ(imag(TypeParam(100.0,0)),imag(&this->Array)->get_data_ptr()[354222]);
}

TYPED_TEST(hoNDArray_elemwise_TestCplx,clamp_minTest){
  fill(&this->Array,TypeParam(-5.7, -4.6));
  this->Array.get_data_ptr()[91] = TypeParam(-101.3, -203.4);
  clamp_min(&this->Array, real(TypeParam(-10.6,0)));
  EXPECT_FLOAT_EQ(real(TypeParam(-5.7,0)),real(&this->Array)->get_data_ptr()[28]);
  EXPECT_FLOAT_EQ(real(TypeParam(-10.6,0)),real(&this->Array)->get_data_ptr()[91]);
  EXPECT_FLOAT_EQ(imag(TypeParam(-5.7,0)),imag(&this->Array)->get_data_ptr()[28]);
  EXPECT_FLOAT_EQ(imag(TypeParam(-10.6,0)),imag(&this->Array)->get_data_ptr()[91]);
}

TYPED_TEST(hoNDArray_elemwise_TestCplx,clamp_maxTest){
  fill(&this->Array,TypeParam(5.7, 4.6));
  this->Array.get_data_ptr()[91] = TypeParam(101.3, 203.4);
  clamp_max(&this->Array,real(TypeParam(10.6,0)));
  EXPECT_FLOAT_EQ(real(TypeParam(5.7,0)),real(&this->Array)->get_data_ptr()[28]);
  EXPECT_FLOAT_EQ(real(TypeParam(10.6,0)),real(&this->Array)->get_data_ptr()[91]);
  EXPECT_FLOAT_EQ(imag(TypeParam(5.7,0)),imag(&this->Array)->get_data_ptr()[28]);
  EXPECT_FLOAT_EQ(imag(TypeParam(10.6,0)),imag(&this->Array)->get_data_ptr()[91]);
}

TYPED_TEST(hoNDArray_elemwise_TestCplx,shrink1Test){
  fill(&this->Array,TypeParam(1.2,1.4));
  shrink1(&this->Array,0.75);
  EXPECT_FLOAT_EQ(real(TypeParam(1.2,1.4)/abs(TypeParam(1.2,1.4)))*std::max(abs(TypeParam(1.2,1.4))-0.75,0.0),real(&this->Array)->get_data_ptr()[125]);
  EXPECT_FLOAT_EQ(imag(TypeParam(1.2,1.4)/abs(TypeParam(1.2,1.4)))*std::max(abs(TypeParam(1.2,1.4))-0.75,0.0),imag(&this->Array)->get_data_ptr()[125]);
  fill(&this->Array,TypeParam(1,1));
  shrink1(&this->Array,2.0);
  EXPECT_FLOAT_EQ(0.0,real(&this->Array)->get_data_ptr()[125]);
  EXPECT_FLOAT_EQ(0.0,imag(&this->Array)->get_data_ptr()[23125]);
}


TYPED_TEST(hoNDArray_elemwise_TestCplx,axpyTest){
    fill(&this->Array,TypeParam(3,4));
    fill(&this->Array2,TypeParam(2,5));
    axpy(TypeParam(2),this->Array2,this->Array);


    EXPECT_FLOAT_EQ(14,imag(&this->Array)->get_data_ptr()[23125]);
    EXPECT_FLOAT_EQ(7,real(&this->Array)->get_data_ptr()[23123]);



}

TYPED_TEST_SUITE(hoNDArray_elemwise_TestCplx4, cplxImplementations);

TYPED_TEST(hoNDArray_elemwise_TestCplx4,shrinkdTest){
  fill(&this->Array,TypeParam(1.2,1.4));
  fill(&this->Array2,real(TypeParam(4.0,4.0)));
  shrinkd(&this->Array,&this->Array2,1.0);
  EXPECT_FLOAT_EQ(real(TypeParam(1.2,1.4)/real(TypeParam(4.0,4.0)))*std::max(4.0-1.0,0.0),real(&this->Array)->get_data_ptr()[125]);
  EXPECT_FLOAT_EQ(imag(TypeParam(1.2,1.4)/imag(TypeParam(4.0,4.0)))*std::max(4.0-1.0,0.0),imag(&this->Array)->get_data_ptr()[125]);
  shrinkd(&this->Array,&this->Array2,8.0);
  EXPECT_FLOAT_EQ(0.0,real(&this->Array)->get_data_ptr()[125]);
  EXPECT_FLOAT_EQ(0.0,imag(&this->Array)->get_data_ptr()[23125]);
}

TYPED_TEST_SUITE(hoNDArray_elemwise_TestCplx2, stdCplxImplementations);

TYPED_TEST(hoNDArray_elemwise_TestCplx2,realToCplxTest){
  fill(&this->Array,TypeParam(3.4,4.2));
  EXPECT_FLOAT_EQ(3.4,real(real_to_complex<TypeParam>(real(&this->Array).get())->get_data_ptr()[33425]));
  EXPECT_FLOAT_EQ(0.0,imag(real_to_complex<TypeParam>(real(&this->Array).get())->get_data_ptr()[33425]));
}

TYPED_TEST(hoNDArray_elemwise_TestCplx2,multiplyConj){
fill(&this->Array,TypeParam(3,1));
fill(&this->Array2,TypeParam(3,2));

multiplyConj(this->Array,this->Array2,this->Array);

EXPECT_FLOAT_EQ(11,real(this->Array[33425]));
EXPECT_FLOAT_EQ(-3,imag(this->Array[33425]));
}

TYPED_TEST_SUITE(hoNDArray_elemwise_TestCplx3, cplxtImplementations);

TYPED_TEST(hoNDArray_elemwise_TestCplx3,realToCplxTest){
  fill(&this->Array,TypeParam(3.4,4.2));
  EXPECT_FLOAT_EQ(3.4,real(real_to_complex<TypeParam>(real(&this->Array).get())->get_data_ptr()[33425]));
  EXPECT_FLOAT_EQ(0.0,imag(real_to_complex<TypeParam>(real(&this->Array).get())->get_data_ptr()[33425]));
}

template <typename T> class hoNDArray_operators_TestReal : public ::testing::Test {
protected:
  virtual void SetUp() {
    size_t vdims[] = {37, 49, 23, 19}; //Using prime numbers for setup because they are messy
    size_t vdims2[] = {37, 49}; //Smaller dimensionality to test batch mode
    dims = std::vector<size_t>(vdims,vdims+sizeof(vdims)/sizeof(size_t));
    dims2 = std::vector<size_t>(vdims2,vdims2+sizeof(vdims2)/sizeof(size_t));
    Array = hoNDArray<T>(dims);
    Array2 = hoNDArray<T>(dims2);
  }
  std::vector<size_t> dims;
  std::vector<size_t> dims2;
  hoNDArray<T> Array;
  hoNDArray<T> Array2;
};

template <typename T> class hoNDArray_operators_TestCplx : public ::testing::Test {
protected:
  virtual void SetUp() {
    dims = {37, 49, 23, 19};
    dims2 = {37,49};
    Array = hoNDArray<T>(dims);
    Array2 = hoNDArray<T>(dims2);
  }
  std::vector<size_t> dims;
  std::vector<size_t> dims2;
  hoNDArray<T> Array;
  hoNDArray<T> Array2;
};

typedef Types<float, double> realImplementations;
typedef Types<std::complex<float>, std::complex<double>, float_complext, double_complext> cplxImplementations;

TYPED_TEST_SUITE(hoNDArray_operators_TestReal, realImplementations);

TYPED_TEST(hoNDArray_operators_TestReal,equalsAddTest1){
  TypeParam v1 = TypeParam(46865.35435);
  TypeParam v2 = TypeParam(13784.34);
  unsigned int idx = 73243;
  fill(&this->Array,v1);
  fill(&this->Array2,v2);
  this->Array += this->Array2;
  auto val = this->Array[idx];
  EXPECT_FLOAT_EQ(v1+v2,this->Array[0]);
  EXPECT_FLOAT_EQ(v1+v2,this->Array[idx]);
}

TYPED_TEST(hoNDArray_operators_TestReal,equalsAddTest2){
  TypeParam v1 = TypeParam(98.4);
  TypeParam v2 = TypeParam(2.2);
  unsigned int idx = 12295;
  fill(&this->Array,v1);
  this->Array += v2;
  EXPECT_FLOAT_EQ(v1+v2,this->Array.get_data_ptr()[idx]);
}

TYPED_TEST(hoNDArray_operators_TestReal,equalsSubtractTest1){
  TypeParam v1 = TypeParam(98475334.34);
  TypeParam v2 = TypeParam(2452.234);
  unsigned int idx = 124999;
  fill(&this->Array,v1);
  fill(&this->Array2,v2);
  this->Array -= this->Array2;
  EXPECT_FLOAT_EQ(v1-v2,this->Array.get_data_ptr()[idx]);
}

TYPED_TEST(hoNDArray_operators_TestReal,equalsSubtractTest2){
  TypeParam v1 = TypeParam(4.4);
  TypeParam v2 = TypeParam(9212.21);
  unsigned int idx = 122131;
  fill(&this->Array,v1);
  this->Array -= v2;
  EXPECT_FLOAT_EQ(v1-v2,this->Array.get_data_ptr()[idx]);
}

TYPED_TEST(hoNDArray_operators_TestReal,equalsMultiplyTest1){
  TypeParam v1 = TypeParam(342.145);
  TypeParam v2 = TypeParam(43545.43);
  unsigned int idx = 12344;
  fill(&this->Array,v1);
  fill(&this->Array2,v2);
  this->Array *= this->Array2;
  EXPECT_FLOAT_EQ(v1*v2,this->Array.get_data_ptr()[idx]);
}

TYPED_TEST(hoNDArray_operators_TestReal,equalsMultiplyTest2){
  TypeParam v1 = TypeParam(43534.443);
  TypeParam v2 = TypeParam(92.842);
  unsigned int idx = 96735;
  fill(&this->Array,v1);
  this->Array *= v2;
  EXPECT_FLOAT_EQ(v1*v2,this->Array.get_data_ptr()[idx]);
}

TYPED_TEST(hoNDArray_operators_TestReal,equalsDivideTest1){
  TypeParam v1 = TypeParam(644.24);
  TypeParam v2 = TypeParam(38564.64);
  unsigned int idx = 98322;
  fill(&this->Array,v1);
  fill(&this->Array2,v2);
  this->Array /= this->Array2;
  EXPECT_FLOAT_EQ(v1/v2,this->Array.get_data_ptr()[idx]);
}

TYPED_TEST(hoNDArray_operators_TestReal,equalsDivideTest2){
  TypeParam v1 = TypeParam(56342.24);
  TypeParam v2 = TypeParam(23434.34);
  unsigned int idx = 12591;
  fill(&this->Array,v1);
  this->Array /= v2;
  EXPECT_FLOAT_EQ(v1/v2,this->Array.get_data_ptr()[idx]);
}

TYPED_TEST_SUITE(hoNDArray_operators_TestCplx, cplxImplementations);

TYPED_TEST(hoNDArray_operators_TestCplx,equalsAddTest1){
  TypeParam v1 = TypeParam(46865.35435, 534544.534523);
  TypeParam v2 = TypeParam(13784.34, 54543543.1243);
  unsigned int idx = 73243;
  fill(&this->Array,v1);
  fill(&this->Array2,v2);
  this->Array += this->Array2;
  EXPECT_FLOAT_EQ(real(v1+v2),real(this->Array.get_data_ptr()[idx]));
  EXPECT_FLOAT_EQ(imag(v1+v2),imag(this->Array.get_data_ptr()[idx]));
}

TYPED_TEST(hoNDArray_operators_TestCplx,equalsAddTest2){
  TypeParam v1 = TypeParam(98.4, 45.34);
  TypeParam v2 = TypeParam(2.2,3.23);
  unsigned int idx = 12925;
  fill(&this->Array,v1);
  this->Array += v2;
  EXPECT_FLOAT_EQ(real(v1+v2),real(this->Array.get_data_ptr()[idx]));
  EXPECT_FLOAT_EQ(imag(v1+v2),imag(this->Array.get_data_ptr()[idx]));
}

TYPED_TEST(hoNDArray_operators_TestCplx,equalsAddTest3){
  TypeParam v1 = TypeParam(98.4, 45.34);
  TypeParam v2 = TypeParam(2.2,0.0);
  unsigned int idx = 12295;
  fill(&this->Array,v1);
  this->Array += real(v2);
  EXPECT_FLOAT_EQ(real(v1+v2),real(this->Array.get_data_ptr()[idx]));
  EXPECT_FLOAT_EQ(imag(v1+v2),imag(this->Array.get_data_ptr()[idx]));
}

TYPED_TEST(hoNDArray_operators_TestCplx,equalsSubtractTest1){
  TypeParam v1 = TypeParam(46865.35435, 534544.534523);
  TypeParam v2 = TypeParam(13784.34, 54543543.1243);
  unsigned int idx = 73243;
  fill(&this->Array,v1);
  fill(&this->Array2,v2);
  this->Array -= this->Array2;
  EXPECT_FLOAT_EQ(real(v1-v2),real(this->Array.get_data_ptr()[idx]));
  EXPECT_FLOAT_EQ(imag(v1-v2),imag(this->Array.get_data_ptr()[idx]));
}

TYPED_TEST(hoNDArray_operators_TestCplx,equalsSubtractTest2){
  TypeParam v1 = TypeParam(98.4, 45.34);
  TypeParam v2 = TypeParam(2.2,3.23);
  unsigned int idx = 12925;
  fill(&this->Array,v1);
  this->Array -= v2;
  EXPECT_FLOAT_EQ(real(v1-v2),real(this->Array.get_data_ptr()[idx]));
  EXPECT_FLOAT_EQ(imag(v1-v2),imag(this->Array.get_data_ptr()[idx]));
}

TYPED_TEST(hoNDArray_operators_TestCplx,equalsSubtractTest3){
  TypeParam v1 = TypeParam(98.4, 45.34);
  TypeParam v2 = TypeParam(2.2,0.0);
  unsigned int idx = 12925;
  fill(&this->Array,v1);
  this->Array -= real(v2);
  EXPECT_FLOAT_EQ(real(v1-v2),real(this->Array.get_data_ptr()[idx]));
  EXPECT_FLOAT_EQ(imag(v1-v2),imag(this->Array.get_data_ptr()[idx]));
}

TYPED_TEST(hoNDArray_operators_TestCplx,equalsMultiplyTest1){
  TypeParam v1 = TypeParam(46865.35435, 534544.534523);
  TypeParam v2 = TypeParam(13784.34, 54543543.1243);
  unsigned int idx = 73243;
  fill(&this->Array,v1);
  fill(&this->Array2,v2);
  this->Array *= this->Array2;
  EXPECT_FLOAT_EQ(real(v1*v2),real(this->Array.get_data_ptr()[idx]));
  EXPECT_FLOAT_EQ(imag(v1*v2),imag(this->Array.get_data_ptr()[idx]));
}

TYPED_TEST(hoNDArray_operators_TestCplx,equalsMultiplyTest2){
  TypeParam v1 = TypeParam(98.4, 45.34);
  TypeParam v2 = TypeParam(2.2,3.23);
  unsigned int idx = 12925;
  fill(&this->Array,v1);
  this->Array *= v2;
  EXPECT_FLOAT_EQ(real(v1*v2),real(this->Array.get_data_ptr()[idx]));
  EXPECT_FLOAT_EQ(imag(v1*v2),imag(this->Array.get_data_ptr()[idx]));
}

TYPED_TEST(hoNDArray_operators_TestCplx,equalsMultiplyTest3){
  TypeParam v1 = TypeParam(98.4, 45.34);
  TypeParam v2 = TypeParam(2.2,0.0);
  unsigned int idx = 12295;
  fill(&this->Array,v1);
  this->Array *= real(v2);
  EXPECT_FLOAT_EQ(real(v1*v2),real(this->Array.get_data_ptr()[idx]));
  EXPECT_FLOAT_EQ(imag(v1*v2),imag(this->Array.get_data_ptr()[idx]));
}

TYPED_TEST(hoNDArray_operators_TestCplx,equalsDivideTest1){
  TypeParam v1 = TypeParam(46865.35435, 534544.534523);
  TypeParam v2 = TypeParam(13784.34, 54543543.1243);
  unsigned int idx = 73243;
  fill(&this->Array,v1);
  fill(&this->Array2,v2);
  this->Array /= this->Array2;
  EXPECT_FLOAT_EQ(real(v1/v2),real(this->Array.get_data_ptr()[idx]));
  EXPECT_FLOAT_EQ(imag(v1/v2),imag(this->Array.get_data_ptr()[idx]));
}

TYPED_TEST(hoNDArray_operators_TestCplx,equalsDivideTest2){
  TypeParam v1 = TypeParam(98.4, 45.34);
  TypeParam v2 = TypeParam(2.2,3.23);
  unsigned int idx = 12295;
  fill(&this->Array,v1);
  this->Array /= v2;
  EXPECT_FLOAT_EQ(real(v1/v2),real(this->Array.get_data_ptr()[idx]));
  EXPECT_FLOAT_EQ(imag(v1/v2),imag(this->Array.get_data_ptr()[idx]));
}

TYPED_TEST(hoNDArray_operators_TestCplx,equalsDivideTest3){
  TypeParam v1 = TypeParam(98.4, 45.34);
  TypeParam v2 = TypeParam(2.2,0.0);
  unsigned int idx = 12295;
  fill(&this->Array,v1);
  this->Array /= real(v2);
  EXPECT_FLOAT_EQ(real(v1/v2),real(this->Array.get_data_ptr()[idx]));
  EXPECT_FLOAT_EQ(imag(v1/v2),imag(this->Array.get_data_ptr()[idx]));
}

TYPED_TEST(hoNDArray_operators_TestCplx,add){
    auto copy = this->Array;
    fill(&this->Array,TypeParam(3,4));
    fill(&this->Array2,TypeParam(2,7));
    add(this->Array,this->Array2,copy);

    EXPECT_FLOAT_EQ(5,real(copy[122]));
    EXPECT_FLOAT_EQ(11,imag(copy[122]));
}
TYPED_TEST(hoNDArray_operators_TestCplx,subtract){
    auto copy = this->Array;
    fill(&this->Array,TypeParam(3,4));
    fill(&this->Array2,TypeParam(2,7));
    subtract(this->Array,this->Array2,copy);

    EXPECT_FLOAT_EQ(1,real(copy[122]));
    EXPECT_FLOAT_EQ(-3,imag(copy[122]));
}
TYPED_TEST(hoNDArray_operators_TestCplx,multiply){
    auto copy = this->Array;
    fill(&this->Array,TypeParam(3,4));
    fill(&this->Array2,TypeParam(2,7));
    multiply(this->Array,this->Array2,copy);

    EXPECT_FLOAT_EQ(-22,real(copy[122]));
    EXPECT_FLOAT_EQ(29,imag(copy[122]));
}
TYPED_TEST(hoNDArray_operators_TestCplx,divide){
    auto copy = this->Array;
    fill(&this->Array,TypeParam(3,4));
    fill(&this->Array2,TypeParam(2,7));
    divide(this->Array,this->Array2,copy);

    EXPECT_FLOAT_EQ(0.6415094339622641,real(copy[122]));
    EXPECT_FLOAT_EQ(-0.24528301886792456,imag(copy[122]));
}
TYPED_TEST(hoNDArray_operators_TestCplx,conjugate){
    auto copy = this->Array;
    fill(&this->Array,TypeParam(3,4));
    conjugate(this->Array,copy);

    EXPECT_FLOAT_EQ(3,real(copy[122]));
    EXPECT_FLOAT_EQ(-4,imag(copy[122]));
}

TYPED_TEST(hoNDArray_operators_TestCplx, inv) {
    auto copy = this->Array;
    fill(&this->Array, TypeParam(3, 4));
    inv(this->Array, copy);

    EXPECT_FLOAT_EQ(0.12, real(copy[122]));
    EXPECT_FLOAT_EQ(-0.16, imag(copy[122]));
}

TYPED_TEST(hoNDArray_operators_TestCplx, sum_over_dimension) {
    fill(&this->Array,TypeParam(1,2));
    sum_over_dimension(this->Array,this->Array2,2);

    auto dims2 = this->Array2.dimensions();
    EXPECT_EQ(this->Array2.get_number_of_dimensions(),this->Array.get_number_of_dimensions());
    EXPECT_EQ(dims2[2],1);

    EXPECT_EQ(dims2[0],this->dims[0]);
    EXPECT_EQ(dims2[1],this->dims[1]);
    EXPECT_EQ(dims2[3],this->dims[3]);

    EXPECT_FLOAT_EQ(1*this->dims[2], real(this->Array2[122]));
    EXPECT_FLOAT_EQ(2*this->dims[2], imag(this->Array2[122]));
}
TYPED_TEST(hoNDArray_operators_TestCplx, nrm2) {
    fill(&this->Array,TypeParam(1,2));
    sum_over_dimension(this->Array,this->Array2,2);

    auto dims2 = this->Array2.dimensions();
    EXPECT_EQ(this->Array2.get_number_of_dimensions(),this->Array.get_number_of_dimensions());
    EXPECT_EQ(dims2[2],1);

    EXPECT_EQ(dims2[0],this->dims[0]);
    EXPECT_EQ(dims2[1],this->dims[1]);
    EXPECT_EQ(dims2[3],this->dims[3]);

    EXPECT_FLOAT_EQ(1*this->dims[2], real(this->Array2[122]));
    EXPECT_FLOAT_EQ(2*this->dims[2], imag(this->Array2[122]));
}
