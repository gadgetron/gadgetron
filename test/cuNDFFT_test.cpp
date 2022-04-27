#include "cuNDFFT.h"
#include "cuNDArray_math.h"
#include "complext.h"
#include <gtest/gtest.h>
#include <boost/random.hpp>

using namespace Gadgetron;
using testing::Types;

template<typename REAL> class cuNDFFT_test : public ::testing::Test {
protected:
	virtual void SetUp(){
		boost::random::mt19937 rng;
		boost::random::uniform_real_distribution<REAL> uni(0,1);
		std::vector<size_t > dimensions(3,128);

		hoNDArray<complext<REAL> > tmp(dimensions);
		complext<REAL>* data = tmp.get_data_ptr();

		for (size_t i = 0; i < tmp.get_number_of_elements(); i++)
			data[i] = complext<REAL>(uni(rng),uni(rng));

		Array = cuNDArray<complext<REAL> >(tmp);
		Array2 = Array;
	}

	cuNDArray<complext<REAL> > Array;

	cuNDArray<complext<REAL> > Array2;

};
typedef Types<float, double> realImplementations;
TYPED_TEST_SUITE(cuNDFFT_test, realImplementations);

TYPED_TEST(cuNDFFT_test,fftNrm2Test){
	cuNDFFT<TypeParam>::instance()->fft(&this->Array);

	EXPECT_NEAR(nrm2(&this->Array2),nrm2(&this->Array),nrm2(&this->Array)*1e-3);

}

TYPED_TEST(cuNDFFT_test,ifftNrm2Test){
	cuNDFFT<TypeParam>::instance()->ifft(&this->Array);

	EXPECT_NEAR(nrm2(&this->Array2),nrm2(&this->Array),nrm2(&this->Array)*1e-3);

}
