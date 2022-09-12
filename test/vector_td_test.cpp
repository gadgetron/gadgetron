/*
 * cuGTBLAS_test.cpp
 *
 *  Created on: Feb 28, 2013
 *      Author: Dae
 */
#include "gtest/gtest.h"


#include <vector>
#include "complext.h"
#include "cuNDArray_elemwise.h"
#include "vector_td_utilities.h"
#include "vector_td_io.h"
#include "cuVector_td_test_kernels.h"
#include <sstream>
using namespace Gadgetron;
using testing::Types;
template <typename T> class vector_td_Test : public ::testing::Test {
	protected:
	 virtual void SetUp() {
		 size_t vdims[] = {37}; //Using prime numbers for setup because they are messy
		 dims= std::vector<size_t>(vdims,vdims+sizeof(vdims)/sizeof(size_t));
		 cuData = cuNDArray<vector_td<T,3> >(&dims);
		 cuData2 = cuNDArray<vector_td<T,3> >(&dims);
	}
	 cuNDArray<vector_td<T,3> > cuData;
	 cuNDArray<vector_td<T,3> > cuData2;
	 std::vector<size_t> dims;


};

//typedef Types<float,double,float_complext,double_complext> Implementations;
typedef Types<float> Implementations;

TYPED_TEST_SUITE(vector_td_Test, Implementations);


TYPED_TEST(vector_td_Test,absTest){

	vector_fill(&this->cuData,vector_td<TypeParam,3>(-2));

	test_abs(&this->cuData);
	vector_td<TypeParam,3> expected(2);
	vector_td<TypeParam,3> result = this->cuData.get_device_ptr()[2];
	EXPECT_EQ(expected,result);
}

TYPED_TEST(vector_td_Test,normTest){
	vector_fill(&this->cuData,vector_td<TypeParam,3>(12.1));

	thrust::device_vector<TypeParam> out = test_norm(&this->cuData);

	EXPECT_FLOAT_EQ(real(20.957814772),out[3]);
}



TYPED_TEST(vector_td_Test,minTest){
	vector_fill(&this->cuData,vector_td<TypeParam,3>(2.2,1.1,5.3));

	thrust::device_vector<TypeParam> out = test_min(&this->cuData);

	EXPECT_FLOAT_EQ(TypeParam(1.1),out[5]);
}

TYPED_TEST(vector_td_Test,maxTest){
	vector_fill(&this->cuData,vector_td<TypeParam,3>(2.2,1.1,5.3));

	thrust::device_vector<TypeParam> out = test_max(&this->cuData);

	EXPECT_FLOAT_EQ(TypeParam(5.3),out[5]);
}


TYPED_TEST(vector_td_Test,aminTest){
	vector_fill(&this->cuData,vector_td<TypeParam,3>(2.2,1.1,5.3));
	vector_fill(&this->cuData2,vector_td<TypeParam,3>(20.2,0.11,5.3));

	boost::shared_ptr<cuNDArray<vector_td<TypeParam,3> > > out = test_amin(&this->cuData,&this->cuData2);
	vector_td<TypeParam,3> expected(2.2,0.11,5.3);
	boost::shared_ptr<hoNDArray<vector_td<TypeParam,3> > > host = out->to_host();
	EXPECT_EQ(expected,host->begin()[35]);
}

TYPED_TEST(vector_td_Test,amin2Test){
	vector_fill(&this->cuData,vector_td<TypeParam,3>(2.2,1.1,5.3));

	boost::shared_ptr<cuNDArray<vector_td<TypeParam,3> > > out = test_amin2(&this->cuData,TypeParam(4));
	vector_td<TypeParam,3> expected(2.2,1.1,4);
	boost::shared_ptr<hoNDArray<vector_td<TypeParam,3> > > host = out->to_host();
	EXPECT_EQ(expected,host->begin()[35]);
}

TYPED_TEST(vector_td_Test,amaxTest){
	vector_fill(&this->cuData,vector_td<TypeParam,3>(2.2,1.1,5.3));
	vector_fill(&this->cuData2,vector_td<TypeParam,3>(20.2,0.11,5.3));

	boost::shared_ptr<cuNDArray<vector_td<TypeParam,3> > > out = test_amax(&this->cuData,&this->cuData2);
	vector_td<TypeParam,3> expected(20.2,1.1,5.3);
	boost::shared_ptr<hoNDArray<vector_td<TypeParam,3> > > host = out->to_host();
	EXPECT_EQ(expected,host->begin()[23]);
}

TYPED_TEST(vector_td_Test,amax2Test){
	vector_fill(&this->cuData,vector_td<TypeParam,3>(2.2,1.1,5.3));

	boost::shared_ptr<cuNDArray<vector_td<TypeParam,3> > > out = test_amax2(&this->cuData,TypeParam(4));
	vector_td<TypeParam,3> expected(4,4,5.3);
	boost::shared_ptr<hoNDArray<vector_td<TypeParam,3> > > host = out->to_host();
	EXPECT_EQ(expected,host->begin()[26]);
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
