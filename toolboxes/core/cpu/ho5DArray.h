#pragma once

#include "hoNDArray.h"

namespace Gadgetron{

template <class T> class ho5DArray : public hoNDArray<T>
{
public:

    typedef hoNDArray<T> BaseClass;

    ho5DArray();
    ho5DArray(size_t sx, size_t sy, size_t sz, size_t ss, size_t sp);
    explicit ho5DArray(const std::vector<size_t>& dimensions);
    ho5DArray(const std::vector<size_t>& dimensions, T* data, bool delete_data_on_destruct = false);
    ho5DArray(size_t sx, size_t sy, size_t sz, size_t ss, size_t sp, T* data, bool delete_data_on_destruct = false);

    virtual ~ho5DArray();

    ho5DArray(const ho5DArray<T>& a);
    ho5DArray<T>& operator=(const ho5DArray<T>& rhs);

    virtual void create(const std::vector<size_t>& dimensions);
    virtual void create(const std::vector<size_t>& dimensions, T* data, bool delete_data_on_destruct = false);

    virtual bool createArray(size_t sx, size_t sy, size_t sz, size_t ss, size_t sp);
    virtual bool createArray(size_t sx, size_t sy, size_t sz, size_t ss, size_t sp, T* data, bool delete_data_on_destruct = false);

    T& operator()(size_t x, size_t y, size_t z, size_t s, size_t p);
    const T& operator()(size_t x , size_t y, size_t z, size_t s, size_t p) const;

    virtual void print(std::ostream& os) const;

protected:

    using BaseClass::dimensions_;
    using BaseClass::offsetFactors_;
    using BaseClass::data_;
    using BaseClass::elements_;
    using BaseClass::delete_data_on_destruct_;

    bool init_accesser();
    bool release_accesser();

    T***** accesser_;
};

}

#include "ho5DArray.hxx"
