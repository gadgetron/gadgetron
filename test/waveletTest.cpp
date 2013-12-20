#include "cuNDArray_blas.h"
#include "cuNDArray_elemwise.h"
#include "cuHaarWaveletOperator.h"
#include "cuNDArray_utils.h"
#include <gtest/gtest.h>
#include <vector>


using namespace Gadgetron;
using testing::Types;



template <typename T> class cuWaveletTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    size_t vdims[] = {2,4}; //Using prime numbers for setup because they are messy
    dims = std::vector<size_t>(vdims,vdims+sizeof(vdims)/sizeof(unsigned int));
    Array = cuNDArray<T>(&dims);
    op.set_domain_dimensions(&dims);
  }
  std::vector<size_t> dims;
  cuNDArray<T> Array;
  cuHaarWaveletOperator<T,2> op;

};

typedef Types<float,double,float_complext,double_complext> implementations;

TYPED_TEST_CASE(cuWaveletTest,implementations);
/*
TYPED_TEST(cuWaveletTest,twobytwo){

	thrust::device_ptr<TypeParam> ptr = this->Array.get_device_ptr();

	for (int i = 0; i < this->Array.get_number_of_elements(); i++) ptr[i] = TypeParam(i+1);

	cuNDArray<TypeParam> out(this->Array);
	clear(&out);
	this->op.mult_M(&this->Array,&out,false);
	TypeParam res = out.get_device_ptr()[0];
	EXPECT_FLOAT_EQ(real(res),real(mean(&this->Array)));

	res = out.get_device_ptr()[1];
	EXPECT_FLOAT_EQ(real(res),real(TypeParam(-2)));

	res = out.get_device_ptr()[2];
	EXPECT_FLOAT_EQ(real(res),real(TypeParam(-4)));

	res = out.get_device_ptr()[3];
		EXPECT_FLOAT_EQ(real(res),real(TypeParam(0)));
}


TYPED_TEST(cuWaveletTest,inv_twobytwo){

	thrust::device_ptr<TypeParam> ptr = this->Array.get_device_ptr();

	for (int i = 0; i < this->Array.get_number_of_elements(); i++) ptr[i] = TypeParam(i+1);

	cuNDArray<TypeParam> out(this->Array);
	clear(&out);
	this->op.mult_MH(&this->Array,&out,false);
	TypeParam res = out.get_device_ptr()[0];
	EXPECT_FLOAT_EQ(real(res),real(TypeParam(2.5)));

	res = out.get_device_ptr()[1];
	EXPECT_FLOAT_EQ(real(res),real(TypeParam(1)));

	res = out.get_device_ptr()[2];
	EXPECT_FLOAT_EQ(real(res),real(TypeParam(0.5)));

	res = out.get_device_ptr()[3];
		EXPECT_FLOAT_EQ(real(res),real(TypeParam(0)));
}
*/

TYPED_TEST(cuWaveletTest,MHM_twobyfour){

	thrust::device_ptr<TypeParam> ptr = this->Array.get_device_ptr();

	for (int i = 0; i < this->Array.get_number_of_elements(); i++) ptr[i] = TypeParam(i+1);

	cuNDArray<TypeParam> out(this->Array);
	clear(&out);
	cuNDArray<TypeParam> result(this->Array);
	clear(&result);
	this->op.mult_M(&this->Array,&out,false);
	this->op.mult_MH(&out,&result,false);

	EXPECT_FLOAT_EQ(real(this->Array.at(0)),real(result.at(0)));
	EXPECT_FLOAT_EQ(real(this->Array.at(1)),real(result.at(1)));
	EXPECT_FLOAT_EQ(real(this->Array.at(2)),real(result.at(2)));
	EXPECT_FLOAT_EQ(real(this->Array.at(3)),real(result.at(3)));
	EXPECT_FLOAT_EQ(real(this->Array.at(4)),real(result.at(4)));
	EXPECT_FLOAT_EQ(real(this->Array.at(5)),real(result.at(5)));
	EXPECT_FLOAT_EQ(real(this->Array.at(6)),real(result.at(6)));
	EXPECT_FLOAT_EQ(real(this->Array.at(7)),real(result.at(7)));

}

TYPED_TEST(cuWaveletTest,M_twobyfour){

	thrust::device_ptr<TypeParam> ptr = this->Array.get_device_ptr();

	for (int i = 0; i < this->Array.get_number_of_elements(); i++) ptr[i] = TypeParam(i+1);

	cuNDArray<TypeParam> result(this->Array);
	clear(&result);
	this->op.mult_M(&this->Array,&result,false);


	EXPECT_FLOAT_EQ(4.5,real(result.at(0)));
	EXPECT_FLOAT_EQ(-4,real(result.at(1)));
	EXPECT_FLOAT_EQ(-1,real(result.at(2)));
	EXPECT_FLOAT_EQ(-1,real(result.at(3)));
	EXPECT_FLOAT_EQ(-2,real(result.at(4)));
	EXPECT_FLOAT_EQ(-2,real(result.at(5)));
	EXPECT_FLOAT_EQ(0,real(result.at(6)));
	EXPECT_FLOAT_EQ(0,real(result.at(7)));

}
