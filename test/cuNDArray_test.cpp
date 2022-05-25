/*
 * cuNDArray_test.cpp
 *
 *  Created on: Mar 1, 2013
 *      Author: Dae
 */


#include "gtest/gtest.h"
#include "cuNDArray.h"
#include <vector>

using namespace Gadgetron;
using testing::Types;


template <typename T> class cuNDArray_Test : public ::testing::Test {
	protected:
	 virtual void SetUp() {
		 unsigned int vdims[] = {37, 49, 23, 19}; //Using prime numbers for setup because they are messy
		 dims= std::vector<unsigned int>(vdims,vdims+sizeof(vdims)/sizeof(unsigned int));
		 Array =cuNDArray<T>(&dims);
		 Array2 =cuNDArray<T>(&dims);


	}
	 std::vector<unsigned int> dims;
	 cuNDArray<T> Array;
	 cuNDArray<T> Array2;

};

typedef Types<float,double,float_complext,double_complext> Implementations;

TYPED_TEST_SUITE(cuNDArray_Test, Implementations);

TYPED_TEST(cuNDArray_Test,fillTest){
	this->Array.fill(TypeParam(1));
	TypeParam res = this->Array.get_device_ptr()[5];
	EXPECT_FLOAT_EQ(1,real(res));
	this->Array.fill(TypeParam(27));
	res = this->Array.get_device_ptr()[42];
	EXPECT_FLOAT_EQ(27,real(res));
}


TYPED_TEST(cuNDArray_Test,clearTest){
	this->Array.fill(TypeParam(1));
	TypeParam res = this->Array.get_device_ptr()[5];
	EXPECT_FLOAT_EQ(1,real(res));
	this->Array.clear();
	res = this->Array.get_device_ptr()[5];
	EXPECT_FLOAT_EQ(0,real(res));
}

TYPED_TEST(cuNDArray_Test,equalsMultiplyTest){
	this->Array.fill(TypeParam(2));
	this->Array2.fill(TypeParam(4));
	this->Array *= this->Array2;
	TypeParam res = this->Array.get_device_ptr()[105];
	EXPECT_FLOAT_EQ(8,real(res));

}

TYPED_TEST(cuNDArray_Test,absTest){
	this->Array.fill(TypeParam(2.2));
	this->Array.abs();
	TypeParam res = this->Array.get_device_ptr()[121];
	EXPECT_FLOAT_EQ(real(res),2.2);
	this->Array.fill(TypeParam(-2.2));
	this->Array.abs();
	res = this->Array.get_device_ptr()[121];
	EXPECT_FLOAT_EQ(real(res),2.2);
}


TYPED_TEST(cuNDArray_Test,sqrtTest){
	this->Array.fill(TypeParam(12.1));
	this->Array.sqrt();
	TypeParam res = this->Array.get_device_ptr()[121];
	EXPECT_FLOAT_EQ(real(res),3.478505426);

}
