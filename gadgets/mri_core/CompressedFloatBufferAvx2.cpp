#include "NHLBICompression.h"

#include <immintrin.h>

using namespace NHLBI;

InstructionSet CompressedFloatBufferAvx2::getInstructionSet()
{
    return InstructionSet::Avx2;
}

// Vector scale, round, and convert to integer
static inline __m256i scale_and_cast_ps(const __m256 _float_val, const __m256 _scale, const __m256 _sign_mask, __m256i& _sign)
{
    const __m256 _half = _mm256_set1_ps(0.5f);

    _sign = _mm256_srai_epi32(
        _mm256_castps_si256(
            _mm256_andnot_ps(
                _sign_mask,
                _float_val)),
        31);

    __m256i _int_val = _mm256_cvtps_epi32(
        _mm256_round_ps(
            _mm256_fmadd_ps(
                _mm256_and_ps(
                    _float_val,
                    _sign_mask),
                _scale,
                _half),
            _MM_FROUND_TO_ZERO));

    return _int_val;
}

// Convert number to compact integer representation
// Literal translation on compacting code
// Into vector instructions
static inline __m256i compact_epi32(const __m256i _int_val, const __m256i _sign, const __m256i _bitmask)
{
    const __m256i _one = _mm256_set1_epi32(1);

    __m256i _compact_val = _mm256_and_si256(
        _mm256_blendv_epi8(
            _int_val,
            _mm256_add_epi32(
                _mm256_xor_si256(
                    _int_val,
                    _bitmask),
                _one),
            _sign),
        _bitmask);

    return _compact_val;
}

void CompressedFloatBufferAvx2::compress(std::vector<float>& d, float tolerance, uint8_t precision_bits)
{
    if (tolerance <= 0 && (precision_bits < 1 || precision_bits > 31))
    {
        throw std::runtime_error("Argument 'precision_bits' must be between 1 and 31 if 'tolerance' is less than 0.");
    }

    elements_ = d.size();

    float* dptr = d.data();
    size_t i = 0;

    const __m256 _sign_mask = _mm256_castsi256_ps(_mm256_set1_epi32(0x7fffffff));

    // Initialize
    float max_val;

    size_t nIterations = elements_ / 32;

    if (nIterations > 0)
    {
        // Using four accumulator registers to:
        // 1. Reduce data dependencies in the CPU pipeline
        // 2. Unroll the loop
        __m256 _max0 = _mm256_and_ps(_mm256_loadu_ps(dptr + 0), _sign_mask);
        __m256 _max1 = _mm256_and_ps(_mm256_loadu_ps(dptr + 8), _sign_mask);
        __m256 _max2 = _mm256_and_ps(_mm256_loadu_ps(dptr + 16), _sign_mask);
        __m256 _max3 = _mm256_and_ps(_mm256_loadu_ps(dptr + 24), _sign_mask);

        i += 32;

        // Vector body
        for (size_t n = 1; n < nIterations; i += 32, n++) {
            _max0 = _mm256_max_ps(_max0,
                _mm256_and_ps(_mm256_loadu_ps(dptr + i + 0), _sign_mask));

            _max1 = _mm256_max_ps(_max1,
                _mm256_and_ps(_mm256_loadu_ps(dptr + i + 8), _sign_mask));

            _max2 = _mm256_max_ps(_max2,
                _mm256_and_ps(_mm256_loadu_ps(dptr + i + 16), _sign_mask));

            _max3 = _mm256_max_ps(_max3,
                _mm256_and_ps(_mm256_loadu_ps(dptr + i + 24), _sign_mask));
        }

        // Reduce
        _max0 = _mm256_max_ps(_max0, _max1);
        _max2 = _mm256_max_ps(_max2, _max3);

        _max0 = _mm256_max_ps(_max0, _max2);

        __m128 _max = _mm_max_ps(
            _mm256_castps256_ps128(_max0),
            _mm256_extractf128_ps(_max0, 1));

        _max = _mm_max_ps(_max,
            _mm_permute_ps(_max, 0xb1));

        _max = _mm_max_ps(_max,
            _mm_permute_ps(_max, 0x7e));

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

        const __m256 _scale = _mm256_set1_ps(scale_);
        const __m256i _bitmask = _mm256_set1_epi32(bitmask);

        // Vector body
        for (size_t n = 0; n < nIterations; i += 32, n++) {
            // Load floating point data
            // Using four sets of registers to:
            // 1. Reduce data dependencies in the CPU pipeline
            // 2. Unroll the loop
            // 3. Keep enough data for the vector pack

            // The data currently will look like this:
            // _float_val0 = af0,af1,af2,af3,af4,af5,af6,af7
            // _float_val1 = bf0,bf1,bf2,bf3,bf4,bf5,bf6,bf7
            // _float_val2 = cf0,cf1,cf2,cf3,cf4,cf5,cf6,cf7
            // _float_val3 = df0,df1,df2,df3,df4,df5,df6,df7

            __m256 _float_val0 = _mm256_loadu_ps(dptr + i + 0);
            __m256 _float_val1 = _mm256_loadu_ps(dptr + i + 8);
            __m256 _float_val2 = _mm256_loadu_ps(dptr + i + 16);
            __m256 _float_val3 = _mm256_loadu_ps(dptr + i + 24);

            // Transposing values allows using vector operations 
            // Instead of packing values horizontally
            // It is convenient to transpose early
            // Since there are more convenient instructions for transposing floats
            // Due to other quirks of the instruction set
            // Lower and upper halves of the AVX registers
            // Are transposed independently here 
            // As two separate 4x4 matrices
            __m256 t0 = _mm256_unpacklo_ps(_float_val0, _float_val1);
            __m256 t1 = _mm256_unpacklo_ps(_float_val2, _float_val3);
            __m256 t2 = _mm256_unpackhi_ps(_float_val0, _float_val1);
            __m256 t3 = _mm256_unpackhi_ps(_float_val2, _float_val3);

            _float_val0 = _mm256_shuffle_ps(t0, t1, 0x44);
            _float_val1 = _mm256_shuffle_ps(t0, t1, 0xEE);
            _float_val2 = _mm256_shuffle_ps(t2, t3, 0x44);
            _float_val3 = _mm256_shuffle_ps(t2, t3, 0xEE);

            // Now the data looks like this:
            // _float_val0 = af0,bf0,cf0,df0,af4,bf4,cf4,df4
            // _float_val1 = af1,bf1,cf1,df1,af5,bf5,cf5,df5
            // _float_val2 = af2,bf2,cf2,df2,af6,bf6,cf6,df6
            // _float_val3 = af3,bf3,cf3,df3,af7,bf7,cf7,df7

            // Record float signs
            __m256i _sign0, _sign1, _sign2, _sign3;

            // Convert to integers
            __m256i _int_val0 = scale_and_cast_ps(_float_val0, _scale, _sign_mask, _sign0);
            __m256i _int_val1 = scale_and_cast_ps(_float_val1, _scale, _sign_mask, _sign1);
            __m256i _int_val2 = scale_and_cast_ps(_float_val2, _scale, _sign_mask, _sign2);
            __m256i _int_val3 = scale_and_cast_ps(_float_val3, _scale, _sign_mask, _sign3);

            // Compact
            __m256i _compact_val0 = compact_epi32(_int_val0, _sign0, _bitmask);
            __m256i _compact_val1 = compact_epi32(_int_val1, _sign1, _bitmask);
            __m256i _compact_val2 = compact_epi32(_int_val2, _sign2, _bitmask);
            __m256i _compact_val3 = compact_epi32(_int_val3, _sign3, _bitmask);

            // Now the data looks like this:
            // ai0,bi0,ci0,di0,ai4,bi4,ci4,di4
            // ai1,bi1,ci1,di1,ai5,bi5,ci5,di5
            // ai2,bi2,ci2,di2,ai6,bi6,ci6,di6
            // ai3,bi3,ci3,di3,ai7,bi7,ci7,di7

            // And we want to write the following to bptr:
            // ai0bi0ci0di0ai4bi4ci4di4ai1bi1ci1di1ai5bi5ci5di5ai2bi2ci2di2ai6bi6ci6di6ai3bi3ci3di3ai7bi7ci7di7

            // The first step is to pack each column together. Packing means ORing
            // shifted values together. But this will require more than 256 bits, so we need to split  
            // split the data out into two sets of 256-bit registers

            // lo:              hi:  
            // ai0,bi0,ci0,di0  ai4,bi4,ci4,di4
            // ai1,bi1,ci1,di1  ai5,bi5,ci5,di5
            // ai2,bi2,ci2,di2  ai6,bi6,ci6,di6
            // ai3,bi3,ci3,di3  ai7,bi7,ci7,di7

            __m256i _extended_lo0 = _mm256_cvtepi32_epi64(_mm256_castsi256_si128(_compact_val0));
            __m256i _extended_hi0 = _mm256_cvtepi32_epi64(_mm256_extracti128_si256(_compact_val0, 1));
            __m256i _extended_lo1 = _mm256_cvtepi32_epi64(_mm256_castsi256_si128(_compact_val1));
            __m256i _extended_hi1 = _mm256_cvtepi32_epi64(_mm256_extracti128_si256(_compact_val1, 1));
            __m256i _extended_lo2 = _mm256_cvtepi32_epi64(_mm256_castsi256_si128(_compact_val2));
            __m256i _extended_hi2 = _mm256_cvtepi32_epi64(_mm256_extracti128_si256(_compact_val2, 1));
            __m256i _extended_lo3 = _mm256_cvtepi32_epi64(_mm256_castsi256_si128(_compact_val3));
            __m256i _extended_hi3 = _mm256_cvtepi32_epi64(_mm256_extracti128_si256(_compact_val3, 1));

            // Now pack the columns:

            // lo:
            // ai0|(ai1 << 1 * bits_)|(ai2 << 2 * bits_)|(ai3 << 3 * bits_),bi0|(bi1 << 1 * bits_)|(bi2 << 2 * bits_)|(bi3 << 3 * bits_),ci0|(ci1 << 1 * bits_)|(ci2 << 2 * bits_)|(ci3 << 3 * bits_),di0|(di1 << 1 * bits_)|(di2 << 2 * bits_)|(di3 << 3 * bits_)

            // hi:
            // ai4|(ai5 << 5 * bits_)|(ai6 << 6 * bits_)|(ai7 << 7 * bits_),bi4|(bi5 << 5 * bits_)|(bi6 << 6 * bits_)|(bi7 << 7 * bits_),ci4|(ci5 << 5 * bits_)|(ci6 << 6 * bits_)|(ci7 << 7 * bits_),di4|(di5 << 5 * bits_)|(di6 << 6 * bits_)|(di7 << 7 * bits_)

            __m256i _pack_lo = _mm256_or_si256(
                _extended_lo0,
                _mm256_slli_epi64(
                    _mm256_or_si256(
                        _extended_lo1,
                        _mm256_slli_epi64(
                            _mm256_or_si256(
                                _extended_lo2,
                                _mm256_slli_epi64(_extended_lo3, bits_)),
                            bits_)),
                    bits_));

            __m256i _pack_hi = _mm256_or_si256(
                _extended_hi0,
                _mm256_slli_epi64(
                    _mm256_or_si256(
                        _extended_hi1,
                        _mm256_slli_epi64(
                            _mm256_or_si256(
                                _extended_hi2,
                                _mm256_slli_epi64(_extended_hi3, bits_)),
                            bits_)),
                    bits_));

            // And finally write the packed values to the bitstream.
            // To restore the correct order after transposition,
            // the chunks from the lower and upper halves are now interspersed
            pack(bptr, b, upshift, bits_x_4, static_cast<uint64_t>(_mm256_extract_epi64(_pack_lo, 0)));
            pack(bptr, b, upshift, bits_x_4, static_cast<uint64_t>(_mm256_extract_epi64(_pack_hi, 0)));
            pack(bptr, b, upshift, bits_x_4, static_cast<uint64_t>(_mm256_extract_epi64(_pack_lo, 1)));
            pack(bptr, b, upshift, bits_x_4, static_cast<uint64_t>(_mm256_extract_epi64(_pack_hi, 1)));
            pack(bptr, b, upshift, bits_x_4, static_cast<uint64_t>(_mm256_extract_epi64(_pack_lo, 2)));
            pack(bptr, b, upshift, bits_x_4, static_cast<uint64_t>(_mm256_extract_epi64(_pack_hi, 2)));
            pack(bptr, b, upshift, bits_x_4, static_cast<uint64_t>(_mm256_extract_epi64(_pack_lo, 3)));
            pack(bptr, b, upshift, bits_x_4, static_cast<uint64_t>(_mm256_extract_epi64(_pack_hi, 3)));
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

void CompressedFloatBufferAvx2::decompress(float* dptr)
{
    size_t i = 0;

    // Literal translation of getValue
    // into vector instructions.
    // Processes 8 floats per iteration.

    const long long int* bptr = reinterpret_cast<long long int*>(comp_.data());
    const uint64_t bitmask = ((1 << bits_) - 1);

    const __m256 _scale = _mm256_set1_ps(scale_);
    const __m256i _zero = _mm256_setzero_si256();
    const __m256i _one = _mm256_set1_epi32(1);
    const __m256i _bitmask_32 = _mm256_set1_epi32(bitmask);
    const __m256i _bitmask_64 = _mm256_set1_epi64x(bitmask);
    const __m256i _bits = _mm256_set1_epi64x((int64_t)bits_);
    const __m256i _sign_mask = _mm256_set1_epi32(1 << (bits_ - 1));
    const __m256i _bits_x_8 = _mm256_slli_epi64(_bits, 3);

    __m256i _bits_x_idx_0 = _mm256_mullo_epi32(
        _mm256_set_epi64x(6, 4, 2, 0),
        _bits);

    __m256i _bits_x_idx_1 = _mm256_mullo_epi32(
        _mm256_set_epi64x(7, 5, 3, 1),
        _bits);

    size_t nIterations = elements_ / 8;

    // Vector body
    for (size_t n = 0; n < nIterations; i += 8, n++) {
        // Load/Gather
        __m256i _sb_0 = _mm256_srli_epi64(_bits_x_idx_0, 3);
        __m256i _sb_1 = _mm256_srli_epi64(_bits_x_idx_1, 3);

        __m256i _b_0 = _mm256_i64gather_epi64(bptr, _sb_0, 1);
        __m256i _b_1 = _mm256_i64gather_epi64(bptr, _sb_1, 1);

        //Shift back down, and mask other bits
        __m256i _upshift_0 = _mm256_sub_epi64(
            _bits_x_idx_0,
            _mm256_slli_epi64(_sb_0, 3));

        __m256i _upshift_1 = _mm256_sub_epi64(
            _bits_x_idx_1,
            _mm256_slli_epi64(_sb_1, 3));

        __m256i _compact_val_0 = _mm256_and_si256(
            _mm256_srlv_epi64(_b_0, _upshift_0),
            _bitmask_64);

        __m256i _compact_val_1 = _mm256_and_si256(
            _mm256_srlv_epi64(_b_1, _upshift_1),
            _bitmask_64);

        // All values are <32 bits now.
        // Interleave into a single register

        __m256i _compact_val_combined = _mm256_or_si256(
            _compact_val_0,
            _mm256_slli_si256(_compact_val_1, 4));

        //Convert back to binary
        __m256i _int_val = _mm256_blendv_epi8(
            _mm256_sub_epi32(
                _zero,
                _mm256_add_epi32(
                    _mm256_xor_si256(
                        _compact_val_combined,
                        _bitmask_32),
                    _one)),
            _compact_val_combined,
            _mm256_cmpeq_epi32(
                _mm256_and_si256(
                    _compact_val_combined,
                    _sign_mask),
                _zero));

        //Convert to float and scale back
        __m256 _float_val = _mm256_div_ps(
            _mm256_cvtepi32_ps(_int_val),
            _scale);

        // Store
        _mm256_storeu_ps(dptr + i, _float_val);

        _bits_x_idx_0 = _mm256_add_epi64(_bits_x_idx_0, _bits_x_8);
        _bits_x_idx_1 = _mm256_add_epi64(_bits_x_idx_1, _bits_x_8);
    }

    for (; i < elements_; i++) {
        dptr[i] = getValue(i);
    }
}