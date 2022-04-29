#include "hoNDFFT.h"
#include "hoNDArray_math.h"
#include "complext.h"
#include <gtest/gtest.h>
#include <boost/random.hpp>
#include <random>

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
TYPED_TEST_SUITE(hoNDFFT_test, realImplementations);

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

TEST(FFTshiftTest,shift1d_inplace){

    auto array = hoNDArray<std::complex<float>>(24);
    array.fill(0.0f);
    array(2) = 1;
    hoNDFFT<float>::instance()->fftshift1D(array);
    EXPECT_EQ(array(14),1.0f);
    EXPECT_EQ(array(2),0.0f);

}


TEST(FFTshiftTest,shift1d){

    auto array = hoNDArray<std::complex<float>>(24);
    array.fill(0.0f);
    array(2) = 1;
    auto output = array;
    output.fill(0.0f);
    hoNDFFT<float>::instance()->fftshift1D(array,output);
    EXPECT_EQ(output(14),1.0f);
    EXPECT_EQ(output(2),0.0f);

}

TEST(FFTshiftTest,shift2d_inplace){

    auto array = hoNDArray<std::complex<float>>(24,8);
    array.fill(0.0f);
    array(2,3) = 1;
    hoNDFFT<float>::instance()->fftshift2D(array);
    EXPECT_EQ(array(14,7),1.0f);
    EXPECT_EQ(array(2,3),0.0f);

}


TEST(FFTshiftTest,shift2d){

    auto array = hoNDArray<std::complex<float>>(24,8);
    array.fill(0.0f);
    array(2,3) = 1;
    auto output = array;
    output.fill(0.0f);
    hoNDFFT<float>::instance()->fftshift2D(array,output);
    EXPECT_EQ(output(14,7),1.0f);
    EXPECT_EQ(output(2,3),0.0f);

}

TEST(FFTshiftTest,shift3d_inplace){

    auto array = hoNDArray<std::complex<float>>(24,8, 26,3);
    array.fill(0.0f);
    array(2,3, 4,0) = 1;
    hoNDFFT<float>::instance()->fftshift3D(array);
    EXPECT_EQ(array(14,7, 17,0),1.0f);
    EXPECT_EQ(array(2,3, 4,0),0.0f);

}


TEST(FFTshiftTest,shift3d){

    auto array = hoNDArray<std::complex<float>>(24,8, 26,3);
    array.fill(0.0f);
    array(2,3, 4,0) = 1;
    auto output = array;
    output.fill(0.0f);
    hoNDFFT<float>::instance()->fftshift3D(array,output);
    EXPECT_EQ(output(14,7,17,0),1.0f);
    EXPECT_EQ(output(2,3, 4,0),0.0f);

}

template<class... INDICES>
static auto make_random_array(INDICES... indices){
    auto array = hoNDArray<std::complex<float>>(indices...);
    std::default_random_engine e1;
    std::uniform_real_distribution<float> dist{};
    for (auto& val : array) val = {dist(e1),dist(e1)};

    return array;
}

TEST(FFTshiftTest, shift3d_random){
    auto array = make_random_array(24,8,26,3);

    const auto array_copy = array;

    hoNDFFT<float>::instance()->fftshift3D(array);
    hoNDFFT<float>::instance()->ifftshift3D(array);

    EXPECT_EQ(array,array_copy);

    auto buffer = array;
    hoNDFFT<float>::instance()->fftshift3D(array,buffer);
    hoNDFFT<float>::instance()->ifftshift3D(buffer,array);

    EXPECT_EQ(array,array_copy);



    hoNDFFT<float>::instance()->ifftshift3D(array, buffer);
    hoNDFFT<float>::instance()->fftshift3D(buffer);

    EXPECT_EQ(buffer,array_copy);

}


