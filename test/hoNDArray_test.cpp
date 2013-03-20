/*
 * hoNDArray_test.cpp
 *
 *  Created on: Mar 1, 2013
 *      Author: Dae
 */

#include "hoNDArray.h"
#include "complext.h"

#include <gtest/gtest.h>
#include <complex>
#include <vector>

using namespace Gadgetron;
using testing::Types;

template <typename T> class hoNDArray_Test : public ::testing::Test {
	protected:
	 virtual void SetUp() {
		 unsigned int vdims[] = {37, 49, 23, 19}; //Using prime numbers for setup because they are messy
		 dims = std::vector<unsigned int>(vdims,vdims+sizeof(vdims)/sizeof(unsigned int));
		 Array = hoNDArray<T>(&dims);
		 Array2 = hoNDArray<T>(&dims);
	}
	 std::vector<unsigned int> dims;
	 hoNDArray<T> Array;
	 hoNDArray<T> Array2;

};

typedef Types<float, double, std::complex<float>, std::complex<double>, float_complext, double_complext> Implementations;

TYPED_TEST_CASE(hoNDArray_Test, Implementations);

TYPED_TEST(hoNDArray_Test,fillTest){
	this->Array.fill(TypeParam(1));
	EXPECT_FLOAT_EQ(1,real(this->Array.get_data_ptr()[5]));
	this->Array.fill(TypeParam(27));
	EXPECT_FLOAT_EQ(27,real(this->Array.get_data_ptr()[42]));
}

TYPED_TEST(hoNDArray_Test,clearTest){
	this->Array.fill(TypeParam(1));
	EXPECT_FLOAT_EQ(1,real(this->Array.get_data_ptr()[5]));
	this->Array.clear();
	EXPECT_FLOAT_EQ(0,real(this->Array.get_data_ptr()[5]));
}

TYPED_TEST(hoNDArray_Test,equalsMultiplyTest){
	this->Array.fill(TypeParam(2));
	this->Array2.fill(TypeParam(4));
	this->Array *= this->Array2;
	EXPECT_FLOAT_EQ(8,real(this->Array.get_data_ptr()[105]));
}


