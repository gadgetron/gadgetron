/*
 * hoCuGTBLAS_test.cpp
 *
 *  Created on: Feb 28, 2013
 *      Author: Dae
 */
#include "gtest/gtest.h"
#include "hoCuNDArray_blas.h"
#include <vector>
using namespace Gadgetron;
using testing::Types;
template <typename T> class hoCuGTBLAS_Test : public ::testing::Test {
	protected:
	 virtual void SetUp() {
		 unsigned int vdims[] = {37, 49, 23, 19}; //Using prime numbers for setup because they are messy
		 dims= std::vector<unsigned int>(vdims,vdims+sizeof(vdims)/sizeof(unsigned int));
		 Array =hoCuNDArray<T>(&dims);
		 Array2 =hoCuNDArray<T>(&dims);


	}
	 std::vector<unsigned int> dims;
	 hoCuNDArray<T> Array;
	 hoCuNDArray<T> Array2;

};

typedef Types<float,double,float_complext,double_complext> Implementations;

TYPED_TEST_SUITE(hoCuGTBLAS_Test, Implementations);


TYPED_TEST(hoCuGTBLAS_Test,dotTest){
	this->Array.fill(TypeParam(1));
	EXPECT_FLOAT_EQ(this->Array.get_number_of_elements(),real(dot(&this->Array,&this->Array)));

	this->Array2.fill(TypeParam(2));
	EXPECT_FLOAT_EQ(this->Array.get_number_of_elements()*2,real(dot(&this->Array,&this->Array2)));
}

TYPED_TEST(hoCuGTBLAS_Test,axpyTest){
	this->Array.fill(TypeParam(71));
	this->Array2.fill(TypeParam(97));
	axpy(TypeParam(11),&this->Array,&this->Array2);

	TypeParam val = this->Array2.get_data_ptr()[10];
	EXPECT_FLOAT_EQ(878,real(val));

}

TYPED_TEST(hoCuGTBLAS_Test,nrm2Test){
	this->Array.fill(TypeParam(1));
	EXPECT_FLOAT_EQ(std::sqrt((double)this->Array.get_number_of_elements()),nrm2(&this->Array));
	this->Array.fill(TypeParam(3));
	EXPECT_FLOAT_EQ(std::sqrt(3.0*3.0*this->Array.get_number_of_elements()),nrm2(&this->Array));
}

TYPED_TEST(hoCuGTBLAS_Test,asumTest){
	this->Array.fill(TypeParam(1));
	EXPECT_FLOAT_EQ(this->Array.get_number_of_elements(),real(asum(&this->Array)));
	this->Array.fill(TypeParam(-3));
	EXPECT_FLOAT_EQ(this->Array.get_number_of_elements()*3,real(asum(&this->Array)));
}

TYPED_TEST(hoCuGTBLAS_Test,aminTest){
	this->Array.fill(TypeParam(100));
	this->Array.get_data_ptr()[23]=TypeParam(-50);
	EXPECT_EQ(23,amin(&this->Array));
	this->Array.get_data_ptr()[48]=TypeParam(2);
	EXPECT_EQ(48,amin(&this->Array));

}
TYPED_TEST(hoCuGTBLAS_Test,amaxTest){
	this->Array.fill(TypeParam(1));
	this->Array.get_data_ptr()[23]=TypeParam(2);
	EXPECT_EQ(23,amax(&this->Array));
	this->Array.get_data_ptr()[48]=TypeParam(-50);
	EXPECT_EQ(48,amax(&this->Array));

}
