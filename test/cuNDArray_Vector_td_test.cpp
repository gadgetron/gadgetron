/*
 * cuGTBLAS_test.cpp
 *
 *  Created on: Feb 28, 2013
 *      Author: Dae
 */
#include "gtest/gtest.h"


#include <vector>
#include "complext.h"
#include "cuNDArray.h"
#include "vector_td_utilities.h"

using namespace Gadgetron;
using testing::Types;
template <typename T> class cuNDArray_vector_td_Test : public ::testing::Test {
	protected:
	 virtual void SetUp() {
		 unsigned int vdims[] = {37}; //Using prime numbers for setup because they are messy
		 dims= std::vector<unsigned int>(vdims,vdims+sizeof(vdims)/sizeof(unsigned int));
		 cuData = cuNDArray<vector_td<T,3> >(dims);
		 cuData.clear();
	}
	 cuNDArray<vector_td<T,3> > cuData;
	 std::vector<unsigned int> dims;


};

//typedef Types<float,double,float_complext,double_complext> Implementations;
typedef Types<float> Implementations;

TYPED_TEST_SUITE(cuNDArray_vector_td_Test, Implementations);

TYPED_TEST(cuNDArray_vector_td_Test,absTest){
	this->cuData.fill(vector_td<TypeParam,3>(-2));
	this->cuData.abs();
	vector_td<TypeParam,3> expected(2);
	vector_td<TypeParam,3> result = this->cuData.get_device_ptr()[2];
	EXPECT_EQ(expected,result);
}

TYPED_TEST(cuNDArray_vector_td_Test,sqrtTest){
	this->cuData.fill(vector_td<TypeParam,3>(12.1));
	this->cuData.sqrt();
	vector_td<TypeParam,3> expected(TypeParam(3.478505426));
	vector_td<TypeParam,3> result = this->cuData.get_device_ptr()[2];
	EXPECT_FLOAT_EQ(result[1],expected[1]);
}
