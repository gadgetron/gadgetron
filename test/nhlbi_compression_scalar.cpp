#undef __AVX2__
#undef __SSE4_1__

#define NHLBI_NAMESPACE Nhlbi_scalar

#include "NHLBICompression.h"
#include <gtest/gtest.h>

using namespace std;

float nhlbi_compression_roundtrip_scalar(vector<float>& input, float* output, uint8_t precision_bits)
{
    NHLBI_NAMESPACE::CompressedBufferFloat compressor(input, -1.0f, precision_bits);

    EXPECT_EQ(compressor.getInstructionSet(), NHLBI_NAMESPACE::CompressedBufferFloat::InstructionSet::Scalar);

    compressor.decompress(output);

    return compressor.getTolerance();
}