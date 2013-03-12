/*
 * hoNDArray_utils_test.cpp
 *
 *  Created on: Feb 28, 2013
 *      Author: Dae
 */
#include "gtest/gtest.h"
#include "hoNDArray_utils.h"
#include <vector>
using namespace Gadgetron;
using testing::Types;
template <typename T> class hoNDArray_utils_Test : public ::testing::Test {
	protected:
	 virtual void SetUp() {
		 unsigned int vdims[] = {37, 49, 23, 19}; //Using prime numbers for setup because they are messy
		 dims= std::vector<unsigned int>(vdims,vdims+sizeof(vdims)/sizeof(unsigned int));
		 Array =hoNDArray<T>(&dims);
		 Array2 =hoNDArray<T>(&dims);


	}
	 std::vector<unsigned int> dims;
	 hoNDArray<T> Array;
	 hoNDArray<T> Array2;



};


template <typename T> class hoNDArray_utils_Test2 : public ::testing::Test {
	protected:
	 virtual void SetUp() {
		 unsigned int vdims[] = {37, 49, 23, 19}; //Using prime numbers for setup because they are messy
		 dims= std::vector<unsigned int>(vdims,vdims+sizeof(vdims)/sizeof(unsigned int));
		 Array =hoNDArray<T>(&dims);
		 Array2 =hoNDArray<T>(&dims);


	}
	 std::vector<unsigned int> dims;
	 hoNDArray<T> Array;
	 hoNDArray<T> Array2;



};
typedef Types<float,double,float_complext,double_complext> Implementations;

typedef Types<float,double> RealImplementations;

TYPED_TEST_CASE(hoNDArray_utils_Test, Implementations);

TYPED_TEST_CASE(hoNDArray_utils_Test2, RealImplementations);

TYPED_TEST(hoNDArray_utils_Test,absTest){
	typedef typename realType<TypeParam>::type REAL;
	this->Array.fill(TypeParam(-5));

	boost::shared_ptr<hoNDArray<REAL> > res = abs(&this->Array);
	EXPECT_FLOAT_EQ(REAL(5),res->get_data_ptr()[28]);
}


TYPED_TEST(hoNDArray_utils_Test2,clampTest){

	this->Array.fill(TypeParam(-5));
	this->Array.get_data_ptr()[91] = TypeParam(101);
	clamp(&this->Array,TypeParam(5),TypeParam(100));
	EXPECT_FLOAT_EQ(5,this->Array.get_data_ptr()[28]);
	EXPECT_FLOAT_EQ(100,this->Array.get_data_ptr()[91]);
}

TYPED_TEST(hoNDArray_utils_Test,sgnTest){
	typedef typename realType<TypeParam>::type REAL;
	this->Array.fill(TypeParam(-5));
	this->Array.get_data_ptr()[91] = TypeParam(101);
	this->Array.get_data_ptr()[191] = TypeParam(0);
	inplace_sgn(&this->Array);
	EXPECT_FLOAT_EQ(REAL(-1),real(this->Array.get_data_ptr()[28]));
	EXPECT_FLOAT_EQ(REAL(1),real(this->Array.get_data_ptr()[91]));
	EXPECT_FLOAT_EQ(REAL(0),real(this->Array.get_data_ptr()[191]));
}




