#include "hoNDHarrWavelet.h"
#include "hoNDRedundantWavelet.h"
#include "hoNDArray_math.h"
#include <gtest/gtest.h>
#include <boost/random.hpp>

using namespace Gadgetron;
using testing::Types;

template<typename REAL> class hoNDWavelet_test : public ::testing::Test
{
protected:
    virtual void SetUp()
    {
        boost::random::mt19937 rng;
        boost::random::uniform_real_distribution<REAL> uni(0,1);
        std::vector<size_t > dimensions(3,128);

        Array = hoNDArray< std::complex<REAL> >(dimensions);
        std::complex<REAL>* data = Array.get_data_ptr();

        for (size_t i = 0; i < Array.get_number_of_elements(); i++)
            data[i] = std::complex<REAL>(uni(rng), uni(rng));
    }

    hoNDArray< std::complex<REAL> > Array;
};

typedef Types<float, double> realImplementations;
TYPED_TEST_SUITE(hoNDWavelet_test, realImplementations);

TYPED_TEST(hoNDWavelet_test, hoNDHarrWaveletTest1D)
{
    Gadgetron::hoNDHarrWavelet< std::complex<TypeParam> > wav;

    hoNDArray< std::complex<TypeParam> > r, rr, diff;

    size_t WavDim = 1;
    size_t level = 3;

    wav.transform(this->Array, r, WavDim, level, true);
    wav.transform(r, rr, WavDim, level, false);

    Gadgetron::subtract(this->Array, rr, diff);

    TypeParam v =  Gadgetron::nrm2(diff);

    EXPECT_NEAR(v, 0, 0.001);
}

TYPED_TEST(hoNDWavelet_test, hoNDHarrWaveletTest2D)
{
    Gadgetron::hoNDHarrWavelet< std::complex<TypeParam> > wav;

    hoNDArray< std::complex<TypeParam> > r, rr, diff;

    size_t WavDim = 2;
    size_t level = 2;

    wav.transform(this->Array, r, WavDim, level, true);
    wav.transform(r, rr, WavDim, level, false);

    Gadgetron::subtract(this->Array, rr, diff);

    TypeParam v = Gadgetron::nrm2(diff);

    EXPECT_NEAR(v, 0, 0.001);
}

TYPED_TEST(hoNDWavelet_test, hoNDHarrWaveletTest3D)
{
    Gadgetron::hoNDHarrWavelet< std::complex<TypeParam> > wav;

    hoNDArray< std::complex<TypeParam> > r, rr, diff;

    size_t WavDim = 3;
    size_t level = 1;

    wav.transform(this->Array, r, WavDim, level, true);
    wav.transform(r, rr, WavDim, level, false);

    Gadgetron::subtract(this->Array, rr, diff);

    TypeParam v = Gadgetron::nrm2(diff);

    EXPECT_NEAR(v, 0, 0.001);
}

TYPED_TEST(hoNDWavelet_test, hoNDRedundantWaveletTest1D)
{
    Gadgetron::hoNDRedundantWavelet< std::complex<TypeParam> > wav;
    wav.compute_wavelet_filter("db2");

    hoNDArray< std::complex<TypeParam> > r, rr, diff;

    size_t WavDim = 1;
    size_t level = 2;

    wav.transform(this->Array, r, WavDim, level, true);
    wav.transform(r, rr, WavDim, level, false);

    Gadgetron::subtract(this->Array, rr, diff);

    TypeParam v = Gadgetron::nrm2(diff);

    EXPECT_NEAR(v, 0, 0.001);
}

TYPED_TEST(hoNDWavelet_test, hoNDRedundantWaveletTest2D)
{
    Gadgetron::hoNDRedundantWavelet< std::complex<TypeParam> > wav;
    wav.compute_wavelet_filter("db2");

    hoNDArray< std::complex<TypeParam> > r, rr, diff;

    size_t WavDim = 2;
    size_t level = 3;

    wav.transform(this->Array, r, WavDim, level, true);
    wav.transform(r, rr, WavDim, level, false);

    Gadgetron::subtract(this->Array, rr, diff);

    TypeParam v = Gadgetron::nrm2(diff);

    EXPECT_NEAR(v, 0, 0.001);
}

TYPED_TEST(hoNDWavelet_test, hoNDRedundantWaveletTest3D)
{
    Gadgetron::hoNDRedundantWavelet< std::complex<TypeParam> > wav;
    wav.compute_wavelet_filter("db2");

    hoNDArray< std::complex<TypeParam> > r, rr, diff;

    size_t WavDim = 3;
    size_t level = 2;

    wav.transform(this->Array, r, WavDim, level, true);
    wav.transform(r, rr, WavDim, level, false);

    Gadgetron::subtract(this->Array, rr, diff);

    TypeParam v = Gadgetron::nrm2(diff);

    EXPECT_NEAR(v, 0, 0.001);
}

