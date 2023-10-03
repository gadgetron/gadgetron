#ifndef NHLBICOMPRESSION_H
#define NHLBICOMPRESSION_H

#include <algorithm>
#include <cstdint>
#include <cassert>
#include <cstring>
#include <exception>
#include <fstream>
#include <iostream>
#include <cmath>
#include <string>
#include <sstream>
#include <vector>

namespace NHLBI {
    enum class InstructionSet
    {
        Native,
        Scalar,
        Sse41,
        Avx2
    };

    class CompressedFloatBuffer
    {
    protected:
        CompressedFloatBuffer()
        {
            tolerance_ = 0.0;
            elements_ = 0;
            bits_ = 0;
            max_val_ = 0.0;
            scale_ = 0.0;
        }

    public:
        virtual ~CompressedFloatBuffer() = default;

        static CompressedFloatBuffer* createCompressedBuffer(InstructionSet instructionSet = InstructionSet::Native);

        virtual InstructionSet getInstructionSet();

        size_t size()
        {
            return elements_;
        }

        size_t getPrecision()
        {
            return bits_;
        }

        float getTolerance()
        {
            return tolerance_;
        }

        float getCompressionRatio()
        {
            return (1.0f * elements_ * sizeof(float)) / comp_.size();
        }

        std::vector<uint8_t> serialize();

        void deserialize(std::vector<uint8_t>& buffer);

        virtual void compress(std::vector<float>& d, float tolerance = -1.0f, uint8_t precision_bits = 16);

        virtual void decompress(float* dptr);

        float operator[](size_t idx)
        {
            return getValue(idx);
        }

        float getValue(size_t idx)
        {
            size_t sb = (idx * bits_) / 8;
            uint64_t* bptr = reinterpret_cast<uint64_t*>(&comp_[sb]);

            size_t upshift = idx * bits_ - sb * 8;

            //Create mask with ones corresponding to current bits
            const uint64_t bitmask = ((1ULL << bits_) - 1) << upshift;

            //Mask other bits and shift back down
            uint64_t compact_val = (*bptr & bitmask) >> upshift;

            //Convert back to binary
            int64_t int_val = uncompact_int(compact_val);

            //Scale back and return
            return int_val / scale_;
        }

    protected:
        uint8_t bits_;
        size_t elements_;
        float tolerance_;
        float max_val_;
        float scale_;
        std::vector<uint8_t> comp_;

        void initialize(float max_val, float tolerance, uint8_t precision_bits);

        // Packing bits into qwords
        // Assuming the packed array is filled sequentially
        inline void pack(uint64_t*& bptr, uint64_t& b, size_t& upshift, const size_t bits, const uint64_t compact_val)
        {
            b |= compact_val << upshift;

            upshift += bits;

            if (upshift >= 64) {
                *bptr++ = b;

                upshift -= 64;
                if (upshift == 0) {
                    b = 0;
                }
                else {
                    size_t downshift = bits - upshift;

                    b = compact_val >> downshift;
                }
            }
        }

        // Convert number to compact integer representation
        // And pack it into qwords
        inline void compact_and_pack(uint64_t*& bptr, uint64_t& b, size_t& upshift, const size_t bits, const uint32_t bitmask, const int32_t int_val)
        {
            uint32_t abs_val = static_cast<uint32_t>(std::abs(int_val));

            if (int_val < 0) {
                abs_val ^= bitmask;
                abs_val += 1;
            }
            abs_val &= bitmask;

            uint64_t compact_val = static_cast<uint64_t>(abs_val);

            pack(bptr, b, upshift, bits, compact_val);
        }

        inline int32_t uncompact_int(uint32_t cbin)
        {
            if (cbin & (1 << (bits_ - 1))) {
                const uint32_t bitmask = ((1 << bits_) - 1);
                int32_t out = static_cast<int32_t>((cbin ^ bitmask) + 1);
                out = -out;
                return out;
            }

            return static_cast<int32_t>(cbin);
        }


        static float getScaleForBits(size_t precision_bits, float max_val)
        {
            uint64_t max_int = (1 << (precision_bits - 1)) - 1;

            // ensure that max_val * scale will never exceed max_int because of floating-point rounding
            int i = 0;
            float scale;
            do {
                scale = (max_int - i) / max_val;
                i++;
            } while (static_cast<uint64_t>(std::round(max_val * scale)) > max_int);

            if (scale == 0) {
                throw std::runtime_error("Please choose a higher value as the 'precision_bits' argument.");
            }

            return scale;
        }
    };

    class CompressedFloatBufferSse41 : public CompressedFloatBuffer
    {
    public:
        virtual InstructionSet getInstructionSet();

        virtual void compress(std::vector<float>& d, float tolerance = -1.0f, uint8_t precision_bits = 16);
    };

    class CompressedFloatBufferAvx2 : public CompressedFloatBuffer
    {
    public:
        virtual InstructionSet getInstructionSet();

        virtual void compress(std::vector<float>& d, float tolerance = -1.0f, uint8_t precision_bits = 16);

        virtual void decompress(float* dptr);
    };
}

#endif //NHLBICOMPRESSION
