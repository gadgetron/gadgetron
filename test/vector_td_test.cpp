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
#include "cuVector_td_test_kernels.h"
#include <sstream>
using namespace Gadgetron;
using testing::Types;
template <typename T> class vector_td_Test : public ::testing::Test {
	protected:
	 virtual void SetUp() {
		 unsigned int vdims[] = {37}; //Using prime numbers for setup because they are messy
		 dims= std::vector<unsigned int>(vdims,vdims+sizeof(vdims)/sizeof(unsigned int));
		 cuData = cuNDArray<vector_td<T,3> >(&dims);
		 cuData.clear();
	}
	 cuNDArray<vector_td<T,3> > cuData;
	 std::vector<unsigned int> dims;


};

//typedef Types<float,double,float_complext,double_complext> Implementations;
typedef Types<float,double> Implementations;

TYPED_TEST_CASE(vector_td_Test, Implementations);


TYPED_TEST(vector_td_Test,absTest){
	this->cuData.fill(vector_td<TypeParam,3>(-2));

	test_abs(&this->cuData);
	vector_td<TypeParam,3> expected(2);
	vector_td<TypeParam,3> result = this->cuData.get_device_ptr()[2];
	EXPECT_EQ(expected,result);
}

TYPED_TEST(vector_td_Test,normTest){
	this->cuData.fill(vector_td<TypeParam,3>(12.1));

	thrust::device_vector<TypeParam> out = test_norm(&this->cuData);

	EXPECT_FLOAT_EQ(real(20.957814772),out[3]);
}



TYPED_TEST(vector_td_Test,minTest){
	this->cuData.fill(vector_td<TypeParam,3>(2.2,1.1,5.3));

	thrust::device_vector<TypeParam> out = test_min(&this->cuData);

	EXPECT_FLOAT_EQ(TypeParam(1.1),out[5]);
}

TYPED_TEST(vector_td_Test,maxTest){
	this->cuData.fill(vector_td<TypeParam,3>(2.2,1.1,5.3));

	thrust::device_vector<TypeParam> out = test_max(&this->cuData);

	EXPECT_FLOAT_EQ(TypeParam(5.3),out[5]);
}

TEST(vector_td,parseTest){
std::string base ="[23,22,25]";
std::stringstream ss(base);

vector_td<float,3> vec;
vector_td<float,3> res(23,22,25);
ss >> vec;

EXPECT_FALSE(ss.fail());
EXPECT_EQ(res,vec);

}


TEST(vector_td,parseEqualTest){
	vector_td<float,3> res(23,22,25);
	std::stringstream ss;
	ss << res;

	vector_td<float,3> vec;

	ss >> vec;

	EXPECT_FALSE(ss.fail());
	EXPECT_EQ(res,vec);

}
