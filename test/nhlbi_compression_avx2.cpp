#define NHLBI_NAMESPACE Nhlbi_avx

#include "NHLBICompression.h"
#include <gtest/gtest.h>

using namespace std;

float nhlbi_compression_roundtrip_avx2(vector<float>& input, float* output, uint8_t precision_bits)
{
    NHLBI_NAMESPACE::CompressedBufferFloat compressor(input, -1.0f, precision_bits);

    // this will fail if the AVX2 instruction set if not available on this CPU.
    EXPECT_EQ(compressor.getInstructionSet(), NHLBI_NAMESPACE::CompressedBufferFloat::InstructionSet::Avx2);

    compressor.decompress(output);

    return compressor.getTolerance();
}