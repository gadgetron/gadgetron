#include <numeric>
#include <random>
#include <iostream>
#include <tuple>
#include <gtest/gtest.h>
#include <string>

#include "NHLBICompression.h"

#define TESTSAMPLES 10000000

using namespace NHLBI;

void fill_random(std::vector<float> &v, float mean = 0.0f, float sigma = 1.0f, int seed = 42)
{
    std::default_random_engine generator(seed);
    std::normal_distribution<float> distribution(mean, sigma);

    std::generate(v.begin(), v.end(), [&distribution, &generator]() { return distribution(generator); });
}

double mean(std::vector<float> &v)
{
    return std::accumulate(v.begin(), v.end(), 0.0, [](double &a, float &b){ return a + b;}) / v.size();
}

double variance(std::vector<float> &v)
{
    auto m = mean(v);
    return std::accumulate(v.begin(), v.end(), 0.0, [&m](double &a, float &b) { return a + (b - m)*(b - m); }) / (v.size() -1);
}

// We are gonna parametrize the tests with seed, mean, sigma
class NHLBICompression : public testing::TestWithParam< std::tuple<int, float, float> > {
};

// Internal test for the helper functions
TEST_P(NHLBICompression, NormDist)
{
    auto [ seed, dmean, dsigma ] = GetParam();

    std::vector<float> noise(TESTSAMPLES);
    fill_random(noise, dmean, dsigma, seed);

    // We are not testing against a value relative to the mean, since the zero mean case would fail
    EXPECT_NEAR(mean(noise), dmean, 0.05);

    // Variance should be within 5%
    auto expected_variance = dsigma*dsigma;
    EXPECT_NEAR(variance(noise), dsigma*dsigma, expected_variance*0.05);
}

// Test for compression error statistics
TEST_P(NHLBICompression, Tolerance)
{
    auto [ seed, signal_mean, signal_sigma ] = GetParam();
    
    std::vector<float> signal(TESTSAMPLES);
    std::vector<float> signal_decompressed(TESTSAMPLES);
    std::vector<float> compression_error(TESTSAMPLES);

    const float signal_variance = signal_sigma*signal_sigma;

    // We will make a fairly big compression error here to make the compression error noticeable relative to noise
    const float tolerance = 0.5*signal_sigma; 

    fill_random(signal, signal_mean, signal_sigma, seed);

    std::unique_ptr<CompressedFloatBuffer> signal_compressed(CompressedFloatBuffer::createCompressedBuffer());
    signal_compressed->compress(signal, tolerance);

    float actual_tolerance = signal_compressed->getTolerance();

    EXPECT_LE(actual_tolerance, tolerance);

    // Let the width of our error distribution be W = 2*tolerance,
    // then variance = W^2/12
    // https://en.wikipedia.org/wiki/Continuous_uniform_distribution
    const float expected_error_variance = ((actual_tolerance*2)*(actual_tolerance*2))/12;

    signal_compressed->decompress(signal_decompressed.data());

    for (size_t i = 0; i < signal.size(); i++)
    {
        compression_error[i] = signal_decompressed[i]-signal[i];
    }

    auto error_variance = variance(compression_error);
    auto decompressed_variance = variance(signal_decompressed);

    // Max error should not exceed the tolerance
    auto max_abs_error = std::abs(*std::max_element(compression_error.begin(), compression_error.end(), [](float a, float b) { return std::abs(a) < std::abs(b); }));
    EXPECT_LE(max_abs_error, actual_tolerance);

    // There should be no error bias
    // This may start to fail on too few samples
    EXPECT_NEAR(mean(compression_error), 0.0f, 1e-3*signal_sigma);

    // The variance of the error should follow a uniform distribution
    EXPECT_NEAR(error_variance, expected_error_variance, expected_error_variance*0.01);

    // The variance of a sum is the sum of the variances
    auto expected_decompressed_variance = signal_variance + expected_error_variance;

    // The resulting decompression variance should not differ by more than 5% of error variance
    EXPECT_NEAR(decompressed_variance, expected_decompressed_variance, expected_error_variance*0.05);

    // Now see how many bits were used for compression and verify that
    // if we use one less bit that the error will exceed our tolerance
    size_t reduced_precision = signal_compressed->getPrecision() - 1;
    std::unique_ptr<CompressedFloatBuffer> signal_over_compressed(CompressedFloatBuffer::createCompressedBuffer());
    signal_over_compressed->compress(signal, -1, reduced_precision);
    signal_over_compressed->decompress(signal_decompressed.data());

    for (size_t i = 0; i < signal.size(); i++)
    {
        compression_error[i] = signal_decompressed[i]-signal[i];
    }
    
    EXPECT_TRUE(std::any_of(compression_error.begin(), compression_error.end(), [tolerance](float err) { return err > tolerance; }));
}

float nhlbi_compression_roundtrip_scalar(std::vector<float>& input, float* output, uint8_t precision_bits)
{
    std::unique_ptr<CompressedFloatBuffer> compressor(CompressedFloatBuffer::createCompressedBuffer(InstructionSet::Scalar));

    compressor->compress(input, -1.0f, precision_bits);

    // this will fail if the AVX2 instruction set if not available on this CPU.
    EXPECT_EQ(compressor->getInstructionSet(), InstructionSet::Scalar);

    compressor->decompress(output);

    return compressor->getTolerance();
}

float nhlbi_compression_roundtrip_sse(std::vector<float>&input, float* output, uint8_t precision_bits)
{
    std::unique_ptr<CompressedFloatBuffer> compressor(CompressedFloatBuffer::createCompressedBuffer(InstructionSet::Sse41));

    compressor->compress(input, -1.0f, precision_bits);

    // this will fail if the AVX2 instruction set if not available on this CPU.
    EXPECT_EQ(compressor->getInstructionSet(), InstructionSet::Sse41);

    compressor->decompress(output);

    return compressor->getTolerance();
}

float nhlbi_compression_roundtrip_avx2(std::vector<float>&input, float* output, uint8_t precision_bits)
{
    std::unique_ptr<CompressedFloatBuffer> compressor(CompressedFloatBuffer::createCompressedBuffer(InstructionSet::Avx2));

    compressor->compress(input, -1.0f, precision_bits);

    // this will fail if the AVX2 instruction set if not available on this CPU.
    EXPECT_EQ(compressor->getInstructionSet(), InstructionSet::Avx2);

    compressor->decompress(output);

    return compressor->getTolerance();
}

TEST_P(NHLBICompression, Roundtrip)
{
    auto [ seed, signal_mean, signal_sigma ] = GetParam();
    
    // Running compression/decompression ~100 times 
    // So reducing the data size
    constexpr size_t elements = TESTSAMPLES / 100;

    std::vector<float> signal(elements);
    std::vector<float> output_scalar(elements);
    std::vector<float> output_sse(elements);
    std::vector<float> output_avx2(elements);

    fill_random(signal, signal_mean, signal_sigma, seed);

    for (uint8_t precision_bits = 4; precision_bits < 32; precision_bits++)
    {
        nhlbi_compression_roundtrip_scalar(signal, output_scalar.data(), precision_bits);

        nhlbi_compression_roundtrip_sse(signal, output_sse.data(), precision_bits);

        nhlbi_compression_roundtrip_avx2(signal, output_avx2.data(), precision_bits);

        SCOPED_TRACE(std::to_string(precision_bits));
        EXPECT_EQ(std::memcmp(output_scalar.data(), output_sse.data(), elements), 0);
        EXPECT_EQ(std::memcmp(output_scalar.data(), output_avx2.data(), elements), 0);
    }
}

INSTANTIATE_TEST_SUITE_P(
    TestMatrixForCompression, 
    NHLBICompression, 
    ::testing::Combine(
        ::testing::Values(42, 287), // Seeds
        ::testing::Values(0.0f, 1000.0f), // Means
        ::testing::Range(2.0f, 10.0f))); // Sigmas

class NHLBICompressionBits : public testing::TestWithParam< std::tuple<int, float, float, int, std::string> > 
{    
public:
    std::map<std::string, std::function<float(std::vector<float>&, float*, uint8_t)>> compressionFunctions =
    {
        { "scalar", nhlbi_compression_roundtrip_scalar },
        { "sse", nhlbi_compression_roundtrip_sse },
        { "avx2", nhlbi_compression_roundtrip_avx2 },
    };
};

// Test for compression equivalence between implementations when precision bits is specified.
TEST_P(NHLBICompressionBits, PrecisionTolerance)
{
    auto [ seed, signal_mean, signal_sigma, precision_bits, compressionFunctionName ] = GetParam();
    
    // So reducing the data size
    constexpr size_t elements = TESTSAMPLES / 10;

    std::vector<float> signal(elements);
    std::vector<float> output(elements);

    fill_random(signal, signal_mean, signal_sigma, seed);

    auto compress = compressionFunctions[compressionFunctionName];
    float tolerance = compress(signal, output.data(), precision_bits);

    for (size_t i = 0; i < signal.size(); i++)
    {
        SCOPED_TRACE(std::to_string(i));
        ASSERT_NEAR(signal[i], output[i], tolerance);
    }
}

INSTANTIATE_TEST_SUITE_P(
    TestMatrixForCompression, 
    NHLBICompressionBits, 
    ::testing::Combine(
        ::testing::Values(10), // Seeds
        ::testing::Values(0.0f, 1000.0f), // Means
        ::testing::Values(1.0f), // Sigmas
        ::testing::Values(4, 5, 7, 8, 9, 15, 16, 17, 21, 25, 31), // bits of precision
        ::testing::Values(
            "scalar"
            ,"sse" 
            ,"avx2"
            ))); 

// Verify that the actual error does not exceed requested tolerance
TEST(NHLBICompression, ToleranceNotExceeded)
{
    std::vector<float> input(TESTSAMPLES);
    std::vector<float> compression_error(TESTSAMPLES);

    // Fill input starting at 100 and decrementing by epsilon each index.
    input[0] = 100;
    for(size_t i = 1; i < input.size(); i++)
    {
        input[i] = nextafterf(input[i - 1], 0);
    }

    const float tolerance = input[0] * 1e-6;

    std::unique_ptr<CompressedFloatBuffer> compressor(CompressedFloatBuffer::createCompressedBuffer(InstructionSet::Avx2));
    compressor->compress(input, tolerance);

    EXPECT_LE(compressor->getTolerance(), tolerance);

    // compute error for each value
    for (size_t i = 0; i < input.size(); i++)
    {
        compression_error[i] = compressor->getValue(i) - input[i];
    }

    // count how many times the error exceeds tolerance
    int count_exceeded = std::count_if(compression_error.begin(), compression_error.end(), [tolerance](float x) { return std::abs(x) > tolerance; });

    EXPECT_EQ(count_exceeded, 0);
}
// Verifies constructor argument validation.
TEST(NHLBICompression, ArgumentValidation)
{
    std::unique_ptr<CompressedFloatBuffer> compressor(CompressedFloatBuffer::createCompressedBuffer());

    std::vector<float> v;
    ASSERT_THROW(compressor->compress(v, -1, 32), std::runtime_error);
    ASSERT_THROW(compressor->compress(v, -1, 0), std::runtime_error);
    
    v = { 1e10 };
    ASSERT_THROW(compressor->compress(v, 1e-5), std::runtime_error);
}

TEST(NHLBICompression, Deserialize)
{
    std::vector<float> signal(TESTSAMPLES);
    fill_random(signal);
    std::unique_ptr<CompressedFloatBuffer> compressor(CompressedFloatBuffer::createCompressedBuffer());
    compressor->compress(signal, -1, 12);
    auto serialized = compressor->serialize();

    std::unique_ptr<CompressedFloatBuffer> decompressor(CompressedFloatBuffer::createCompressedBuffer());
    decompressor->deserialize(serialized);
    EXPECT_EQ(compressor->size(), decompressor->size());
    for(size_t i = 0; i < compressor->size(); i++)
    {
        ASSERT_FLOAT_EQ(compressor->getValue(i), decompressor->getValue(i));
    }
}