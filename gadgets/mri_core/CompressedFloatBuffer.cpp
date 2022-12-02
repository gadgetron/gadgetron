#include "NHLBICompression.h"

#include "cpuisa.h"

using namespace NHLBI;

CompressedFloatBuffer* CompressedFloatBuffer::createCompressedBuffer(InstructionSet instructionSet)
{
    switch (instructionSet)
    {
    case InstructionSet::Native:
    case InstructionSet::Avx2:
        if (CPU_supports_AVX2())
        {
            return new CompressedFloatBufferAvx2;
        }
        // FALLTHROUGH
    case InstructionSet::Sse41:
        if (CPU_supports_SSE41())
        {
            return new CompressedFloatBufferSse41;
        }
        break;
    }

    return new CompressedFloatBuffer;
}

#pragma pack(push, 1)
struct CompressionHeader
{
    uint64_t elements_;
    float scale_;
    uint8_t bits_;
};
#pragma pack(pop)

std::vector<uint8_t> CompressedFloatBuffer::serialize()
{
    // extra padding was allocated at the end because we use uint64_t to write to the buffer.
    size_t buffer_size = comp_.size() - sizeof(uint64_t) + sizeof(uint8_t);

    std::vector<uint8_t> out(buffer_size + sizeof(CompressionHeader), 0);
    CompressionHeader h;
    h.elements_ = this->elements_;
    h.scale_ = this->scale_;
    h.bits_ = static_cast<uint8_t>(this->bits_);
    memcpy(&out[0], &h, sizeof(CompressionHeader));
    memcpy(&out[sizeof(CompressionHeader)], &comp_[0], buffer_size);
    return out;
}

void CompressedFloatBuffer::deserialize(std::vector<uint8_t>& buffer)
{
    if (buffer.size() <= sizeof(CompressionHeader)) {
        throw std::runtime_error("Invalid buffer size");
    }

    CompressionHeader h;
    memcpy(&h, &buffer[0], sizeof(CompressionHeader));

    size_t bytes_needed = static_cast<size_t>(std::ceil((h.bits_ * h.elements_) / 8.0f));
    if (bytes_needed != (buffer.size() - sizeof(CompressionHeader))) {
        throw std::runtime_error("Incorrect number of bytes in buffer");
    }

    this->bits_ = h.bits_;
    this->elements_ = h.elements_;
    this->scale_ = h.scale_;
    this->tolerance_ = 0.5f / h.scale_;
    this->comp_.resize(bytes_needed, 0);

    memcpy(&comp_[0], &buffer[sizeof(CompressionHeader)], bytes_needed);
}

InstructionSet CompressedFloatBuffer::getInstructionSet()
{
    return InstructionSet::Scalar;
}

void NHLBI::CompressedFloatBuffer::initialize(float max_val, float tolerance, uint8_t precision_bits)
{
    max_val_ = max_val;

    float epsilon = max_val_ - std::nextafterf(max_val_, 0);

    if (tolerance > 0) {
        // Determine how many bits we need to use to achieve this tolerance.
        // Take floating-point precision loss into account to ensure the actual
        // error does not exceed the tolerance.

        float scale = 0.5 / (tolerance - epsilon);
        uint64_t max_int = static_cast<uint64_t>(std::ceil(scale * max_val_));
        precision_bits = 0;
        while (max_int) {
            precision_bits++;
            max_int = max_int >> 1;
        }
        precision_bits++; //Signed

        if (precision_bits > 31) {
            throw std::runtime_error("Tolerance is too small relative to the largest item in the buffer.");
        }

        // Now we'll maximize the scale we can achieve with these bits
        // (which minimizes the error)
    }

    bits_ = precision_bits;

    scale_ = getScaleForBits(precision_bits, max_val_);

    // The actual tolerance (error) has contributions from integer rounding
    // and the floating-point precision loss.

    tolerance_ = 0.5f / scale_ + epsilon;

    size_t bytes_needed = static_cast<size_t>(std::ceil((bits_ * elements_) / 8.0f))
        + sizeof(uint64_t) - sizeof(uint8_t); // we write using uint64_t, so ensure there is sufficient space for the last write.

    comp_.resize(bytes_needed, 0);
}

void CompressedFloatBuffer::compress(std::vector<float>& d, float tolerance, uint8_t precision_bits)
{
    if (tolerance <= 0 && (precision_bits < 1 || precision_bits > 31))
    {
        throw std::runtime_error("Argument 'precision_bits' must be between 1 and 31 if 'tolerance' is less than 0.");
    }

    elements_ = d.size();

    float* dptr = d.data();

    float max_val = *dptr;

    for (size_t i = 1; i < elements_; i++) {
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

    for (size_t i = 0; i < elements_; i++) {
        //Convert number to compact integeter representation
        int32_t int_val = static_cast<int32_t>(std::round(d[i] * scale_));

        compact_and_pack(bptr, b, upshift, bits_, bitmask, int_val);
    }

    // Writing out the last (incomplete) qword
    if (upshift > 0) {
        *bptr = b;
    }
}

void CompressedFloatBuffer::decompress(float* dptr)
{
    size_t i = 0;

    for (; i < elements_; i++) {
        dptr[i] = getValue(i);
    }
}

