#undef __AVX2__

#define NHLBI_NAMESPACE Nhlbi_sse41

#include "NHLBICompression.h"
#include <gtest/gtest.h>

using namespace std;

float nhlbi_compression_roundtrip_sse(vector<float>& input, float* output, uint8_t precision_bits)
{
    NHLBI_NAMESPACE::CompressedBufferFloat compressor(input, -1.0f, precision_bits);

    EXPECT_EQ(compressor.getInstructionSet(), NHLBI_NAMESPACE::CompressedBufferFloat::InstructionSet::Sse41);

    compressor.decompress(output);
    
    return compressor.getTolerance();
}