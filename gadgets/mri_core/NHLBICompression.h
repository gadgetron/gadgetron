#ifndef NHLBICOMPRESSION_H
#define NHLBICOMPRESSION_H

#include <cstdint>
#include <exception>
#include <fstream>
#include <iostream>
#include <cmath>
#include <string>
#include <sstream>
#include <vector>

#pragma pack(push, 1)
struct CompressionHeader
{
    uint64_t elements_;
    float scale_;
    uint8_t bits_;
};
#pragma pack(pop)

template <typename T> class CompressedBuffer
{

public:
    CompressedBuffer()
    {
        tolerance_ = 0.0;
        elements_ = 0;
        bits_ = 0;
        max_val_ = 0.0;
        scale_ = 0.0;
    }
    
    CompressedBuffer(std::vector<T>& d, T tolerance = -1.0, uint8_t precision_bits = 32)
    {
        auto comp_func = [](T a, T b) { return std::abs(a) < std::abs(b); };
        max_val_ = *std::max_element(d.begin(), d.end(), comp_func);

        if (tolerance > 0) {
            tolerance_ = tolerance;
            scale_ = 0.5/tolerance_;
            uint64_t max_int = static_cast<uint64_t>(std::ceil(std::abs(scale_*max_val_+1)));
            bits_ = 0;
            while (max_int) {
                bits_++;
                max_int = max_int>>1;
            }
            bits_++; //Signed
        } else {
            bits_ = precision_bits;
            uint64_t max_int = (1<<(bits_-1))-1;
            scale_ = (max_int-1)/max_val_;
            tolerance_ = 0.5/scale_;
        }

        elements_ = d.size();
        size_t bytes_needed = static_cast<size_t>(std::ceil((bits_*elements_)/8.0f));
        comp_.resize(bytes_needed, 0);

        for (size_t i = 0; i < d.size(); i++) {
            setValue(i,d[i]);
        }
    }

    float operator[](size_t idx)
    {
        return getValue(idx);
    }

    size_t size()
    {
        return elements_;
    }

    size_t getPrecision()
    {
        return bits_;
    }

    T getCompressionRatio()
    {
        return (1.0*elements_*sizeof(T))/comp_.size();
    }

    std::vector<uint8_t> serialize()
    {
        std::vector<uint8_t> out(comp_.size()+sizeof(CompressionHeader),0);
        CompressionHeader h;
        h.elements_ = this->elements_;
        h.scale_ = this->scale_;
        h.bits_ = static_cast<uint8_t>(this->bits_);
        memcpy(&out[0],&h, sizeof(CompressionHeader));
        memcpy(&out[sizeof(CompressionHeader)], &comp_[0], comp_.size());
        return out;
    }

    void deserialize(std::vector<uint8_t>& buffer)
    {
        if (buffer.size() <= sizeof(CompressionHeader)) {
            throw std::runtime_error("Invalid buffer size");
        }

        CompressionHeader h;
        memcpy(&h, &buffer[0], sizeof(CompressionHeader));
        
        size_t bytes_needed = static_cast<size_t>(std::ceil((h.bits_*h.elements_)/8.0f));
        if (bytes_needed != (buffer.size()-sizeof(CompressionHeader))) {
            throw std::runtime_error("Incorrect number of bytes in buffer");
        }

        this->bits_ = h.bits_;
        this->elements_ = h.elements_;
        this->scale_ = h.scale_;
        this->tolerance_ = 0.5/h.scale_;
        this->comp_.resize(bytes_needed,0);

        memcpy(&comp_[0], &buffer[sizeof(CompressionHeader)], bytes_needed);
    }

private:
    size_t bits_;
    size_t elements_;
    T tolerance_;
    T max_val_;
    T scale_;
    std::vector<uint8_t> comp_;

    void setValue(size_t idx, T v)
    {
        size_t sb = (idx*bits_)/8;

        uint64_t* bptr = reinterpret_cast<uint64_t*>(&comp_[sb]);
        
        size_t upshift = idx*bits_-sb*8;

        //Create mask with ones corresponding to current bits
        const uint64_t bitmask = ((1<<bits_)-1)<<upshift;

        //Convert number to compact integeter representation
        int64_t int_val = static_cast<int64_t>(std::round(v*scale_));
        uint64_t compact_val = compact_int(int_val);
        *bptr = ((*bptr) & (~bitmask)) | (compact_val << upshift);         
    }

    float getValue(size_t idx)
    {
        size_t sb = (idx*bits_)/8;
        uint64_t* bptr = reinterpret_cast<uint64_t*>(&comp_[sb]);
        
        size_t upshift = idx*bits_-sb*8;

        //Create mask with ones corresponding to current bits
        const uint64_t bitmask = ((1<<bits_)-1)<<upshift;

        //Mask other bits and shift back down
        uint64_t compact_val =  (*bptr & bitmask)>>upshift;

        //Convert back to binary
        int64_t int_val = uncompact_int(compact_val);
        
        //Scale back and return
        return int_val / scale_;
    }

    uint64_t compact_int(int64_t bin)
    {

        uint64_t abs_val = static_cast<uint64_t>(std::abs(bin));
        
        if (bin < 0) {
            const uint64_t bitmask = ((1<<bits_)-1);
            abs_val ^= bitmask;
            abs_val += 1;
        }
        return abs_val;
    }

    int64_t uncompact_int(uint64_t cbin)
    {
        if (cbin & (1<<(bits_-1))) {
            const uint64_t bitmask = ((1<<bits_)-1);
            int64_t out = static_cast<int64_t>((cbin ^ bitmask)+1);
            out = -out;
            return out;
        }

        return static_cast<int64_t>(cbin);
    }
};


#endif //NHLBICOMPRESSION


