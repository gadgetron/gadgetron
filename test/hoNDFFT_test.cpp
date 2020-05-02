#include "hoNDFFT.h"
#include "hoNDArray_math.h"
#include "complext.h"
#include <gtest/gtest.h>
#include <boost/random.hpp>

using namespace Gadgetron;
using testing::Types;

template<typename REAL> class hoNDFFT_test : public ::testing::Test {
protected:
	virtual void SetUp(){
		boost::random::mt19937 rng;
		boost::random::uniform_real_distribution<REAL> uni(0,1);
		std::vector<size_t > dimensions(3,128);

		Array = hoNDArray<complext<REAL> >(dimensions);
		complext<REAL>* data = Array.get_data_ptr();

		for (size_t i = 0; i < Array.get_number_of_elements(); i++)
			data[i] = complext<REAL>(uni(rng),uni(rng));

		Array2 = Array;
	}

	hoNDArray<complext<REAL> > Array;

	hoNDArray<complext<REAL> > Array2;

};
typedef Types<float, double> realImplementations;
TYPED_TEST_CASE(hoNDFFT_test, realImplementations);

TYPED_TEST(hoNDFFT_test,fftNrm2Test){
	hoNDFFT<TypeParam>::instance()->fft(&this->Array);

	EXPECT_NEAR(nrm2(&this->Array2),nrm2(&this->Array),nrm2(&this->Array)*1e-2);

}


TYPED_TEST(hoNDFFT_test,ifftNrm2Test){
	hoNDFFT<TypeParam>::instance()->ifft(&this->Array);

	EXPECT_NEAR(nrm2(&this->Array2),nrm2(&this->Array),nrm2(&this->Array)*1e-2);

}
TYPED_TEST(hoNDFFT_test,fft1Nrm2Test){
    hoNDFFT<TypeParam>::instance()->fft(&this->Array,0);

    EXPECT_NEAR(nrm2(&this->Array2),nrm2(&this->Array),nrm2(&this->Array)*1e-2);

}
