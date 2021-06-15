#include "NHLBICompression.h"

#include <immintrin.h>
using namespace NHLBI;


InstructionSet CompressedFloatBufferSse41::getInstructionSet()
{
    return InstructionSet::Sse41;
}

static inline __m128i scale_and_cast_ps(const __m128 _float_val, const __m128 _scale, const __m128 _sign_mask, __m128i& _sign)
{
    const __m128 _half = _mm_set1_ps(0.5f);

    _sign = _mm_srai_epi32(
        _mm_castps_si128(
            _mm_andnot_ps(
                _sign_mask,
                _float_val)),
        31);

    __m128i _int_val = _mm_cvtps_epi32(
        _mm_round_ps(
            _mm_add_ps(
                _mm_mul_ps(
                    _mm_and_ps(
                        _float_val,
                        _sign_mask),
                    _scale),
                _half),
            _MM_FROUND_TO_ZERO));

    return _int_val;
}

static inline __m128i compact_epi32(const __m128i _int_val, const __m128i _sign, const __m128i _bitmask)
{
    const __m128i _one = _mm_set1_epi32(1);

    __m128i _compact_val = _mm_and_si128(
        _mm_blendv_epi8(
            _int_val,
            _mm_add_epi32(
                _mm_xor_si128(
                    _int_val,
                    _bitmask),
                _one),
            _sign),
        _bitmask);

    return _compact_val;
}

void CompressedFloatBufferSse41::compress(std::vector<float>& d, float tolerance, uint8_t precision_bits)
{
    if (tolerance <= 0 && (precision_bits < 1 || precision_bits > 31))
    {
        throw std::runtime_error("Argument 'precision_bits' must be between 1 and 31 if 'tolerance' is less than 0.");
    }

    elements_ = d.size();

    float* dptr = d.data();
    size_t i = 0;

    const __m128 _sign_mask = _mm_castsi128_ps(_mm_set1_epi32(0x7fffffff));

    // Initialize
    float max_val;

    size_t nIterations = elements_ / 16;

    if (nIterations > 0)
    {
        // Initialize
        // Using four accumulator registers to:
        // 1. Reduce data dependencies in the CPU pipeline
        // 2. Unroll the loop
        __m128 _max0 = _mm_and_ps(_mm_loadu_ps(dptr + 0), _sign_mask);
        __m128 _max1 = _mm_and_ps(_mm_loadu_ps(dptr + 4), _sign_mask);
        __m128 _max2 = _mm_and_ps(_mm_loadu_ps(dptr + 8), _sign_mask);
        __m128 _max3 = _mm_and_ps(_mm_loadu_ps(dptr + 12), _sign_mask);

        i += 16;

        // Vector body
        for (size_t n = 1; n < nIterations; i += 16, n++) {
            _max0 = _mm_max_ps(_max0,
                _mm_and_ps(_mm_loadu_ps(dptr + i + 0), _sign_mask));

            _max1 = _mm_max_ps(_max1,
                _mm_and_ps(_mm_loadu_ps(dptr + i + 4), _sign_mask));

            _max2 = _mm_max_ps(_max2,
                _mm_and_ps(_mm_loadu_ps(dptr + i + 8), _sign_mask));

            _max3 = _mm_max_ps(_max3,
                _mm_and_ps(_mm_loadu_ps(dptr + i + 12), _sign_mask));
        }

        // Reduce
        _max0 = _mm_max_ps(_max0, _max1);
        _max2 = _mm_max_ps(_max2, _max3);

        __m128 _max = _mm_max_ps(_max0, _max2);

        _max = _mm_max_ps(_max,
            _mm_shuffle_ps(_max, _max, 0xb1));

        _max = _mm_max_ps(_max,
            _mm_shuffle_ps(_max, _max, 0x7e));

        max_val = _mm_cvtss_f32(_max);
    }
    else
    {
        max_val = dptr[i++];
    }

    // Short tail
    for (; i < elements_; i++) {
        float float_val = std::abs(d[i]);

        if (max_val < float_val) {
            max_val = float_val;
        }
    }

    initialize(max_val, tolerance, precision_bits);

    uint64_t* bptr = reinterpret_cast<uint64_t*>(comp_.data());

    //Create mask with ones corresponding to current bits
    const uint32_t bitmask = ((1 << bits_) - 1);

    uint64_t b = 0;
    size_t upshift = 0;

    i = 0;

    if (bits_ <= 16) {
        size_t bits_x_4 = bits_ * 4;

        const __m128 _scale = _mm_set1_ps(scale_);
        const __m128i _bitmask = _mm_set1_epi32(bitmask);

        // Vector body
        for (size_t n = 0; n < nIterations; i += 16, n++) {
            // Load floating point data
            // Using four sets of registers to:
            // 1. Reduce data dependencies in the CPU pipeline
            // 2. Unroll the loop
            // 3. Keep enough data for the vector pack

            // The data currently will look like this:
            // _float_val0 = af0,af1,af2,af3
            // _float_val1 = bf0,bf1,bf2,bf3
            // _float_val2 = cf0,cf1,cf2,cf3
            // _float_val3 = df0,df1,df2,df3

            __m128 _float_val0 = _mm_loadu_ps(dptr + i + 0);
            __m128 _float_val1 = _mm_loadu_ps(dptr + i + 4);
            __m128 _float_val2 = _mm_loadu_ps(dptr + i + 8);
            __m128 _float_val3 = _mm_loadu_ps(dptr + i + 12);

            // Transposing values allows using vector operations 
            // Instead of packing values horizontally
            // It is convenient to transpose early
            // Since there are more convenient instructions for transposing floats

            _MM_TRANSPOSE4_PS(_float_val0, _float_val1, _float_val2, _float_val3);

            // Record float signs
            __m128i _sign0, _sign1, _sign2, _sign3;

            // Convert to integers
            __m128i _int_val0 = scale_and_cast_ps(_float_val0, _scale, _sign_mask, _sign0);
            __m128i _int_val1 = scale_and_cast_ps(_float_val1, _scale, _sign_mask, _sign1);
            __m128i _int_val2 = scale_and_cast_ps(_float_val2, _scale, _sign_mask, _sign2);
            __m128i _int_val3 = scale_and_cast_ps(_float_val3, _scale, _sign_mask, _sign3);

            // Compact
            __m128i _compact_val0 = compact_epi32(_int_val0, _sign0, _bitmask);
            __m128i _compact_val1 = compact_epi32(_int_val1, _sign1, _bitmask);
            __m128i _compact_val2 = compact_epi32(_int_val2, _sign2, _bitmask);
            __m128i _compact_val3 = compact_epi32(_int_val3, _sign3, _bitmask);

            // Now we have this:
            // ai0,bi0,ci0,di0
            // ai1,bi1,ci1,di1
            // ai2,bi2,ci2,di2
            // ai3,bi3,ci3,di3

            // And we want to write the following to bptr:
            // ai0ai1ai2ai3bi0bi1bi2bi3ci0ci1ci2ci3

            // The first step is to pack each column together. Packing means ORing
            // shifted values together. But this will require more than 128 bits, so we need to split  
            // split the data out into two sets of 128-bit registers

            // lo:      hi:
            // ai0,bi0  ci0,di0
            // ai1,bi1  ci1,di1
            // ai2,bi2  ci2,di2
            // ai3,bi3  ci3,di3

            __m128i _extended_lo0 = _mm_cvtepi32_epi64(_compact_val0);
            __m128i _extended_hi0 = _mm_cvtepi32_epi64(_mm_srli_si128(_compact_val0, 8));
            __m128i _extended_lo1 = _mm_cvtepi32_epi64(_compact_val1);
            __m128i _extended_hi1 = _mm_cvtepi32_epi64(_mm_srli_si128(_compact_val1, 8));
            __m128i _extended_lo2 = _mm_cvtepi32_epi64(_compact_val2);
            __m128i _extended_hi2 = _mm_cvtepi32_epi64(_mm_srli_si128(_compact_val2, 8));
            __m128i _extended_lo3 = _mm_cvtepi32_epi64(_compact_val3);
            __m128i _extended_hi3 = _mm_cvtepi32_epi64(_mm_srli_si128(_compact_val3, 8));

            // Now pack the columns:

            // lo:
            // ai0|(ai1 << 1 * bits_)|(ai2 << 2 * bits_)|(ai3 << 3 * bits_),bi0|(ai1 << 1 * bits_)|(bi2 << 2 * bits_)|(bi3 << 3 * bits_)

            // hi:
            // ci0|(ci1 << 1 * bits_)|(ci2 << 2 * bits_)|(ci3 << 3 * bits_),di0|(di1 << 1 * bits_)|(di2 << 2 * bits_)|(di3 << 3 * bits_)

            __m128i _pack_lo = _mm_or_si128(
                _extended_lo0,
                _mm_slli_epi64(
                    _mm_or_si128(
                        _extended_lo1,
                        _mm_slli_epi64(
                            _mm_or_si128(
                                _extended_lo2,
                                _mm_slli_epi64(_extended_lo3, bits_)),
                            bits_)),
                    bits_));

            __m128i _pack_hi = _mm_or_si128(
                _extended_hi0,
                _mm_slli_epi64(
                    _mm_or_si128(
                        _extended_hi1,
                        _mm_slli_epi64(
                            _mm_or_si128(
                                _extended_hi2,
                                _mm_slli_epi64(_extended_hi3, bits_)),
                            bits_)),
                    bits_));

            // And finally write the packed values to the bitstream.
            pack(bptr, b, upshift, bits_x_4, static_cast<uint64_t>(_mm_extract_epi64(_pack_lo, 0)));
            pack(bptr, b, upshift, bits_x_4, static_cast<uint64_t>(_mm_extract_epi64(_pack_lo, 1)));
            pack(bptr, b, upshift, bits_x_4, static_cast<uint64_t>(_mm_extract_epi64(_pack_hi, 0)));
            pack(bptr, b, upshift, bits_x_4, static_cast<uint64_t>(_mm_extract_epi64(_pack_hi, 1)));
        }
    }


    // Short tail
    for (; i < elements_; i++) {
        //Convert number to compact integeter representation
        int32_t int_val = static_cast<int32_t>(std::round(d[i] * scale_));

        compact_and_pack(bptr, b, upshift, bits_, bitmask, int_val);
    }

    // Writing out the last (incomplete) qword
    if (upshift > 0) {
        *bptr = b;
    }
}
