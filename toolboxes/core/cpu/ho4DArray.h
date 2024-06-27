#pragma once

#include "hoNDArray.h"

namespace Gadgetron{

template <class T> class ho4DArray : public hoNDArray<T>
{
public:

    typedef hoNDArray<T> BaseClass;

    ho4DArray();
    ho4DArray(size_t sx, size_t sy, size_t sz, size_t ss);
    explicit ho4DArray(const std::vector<size_t>& dimensions);
    ho4DArray(const std::vector<size_t>& dimensions, T* data, bool delete_data_on_destruct = false);
    ho4DArray(size_t sx, size_t sy, size_t sz, size_t ss, T* data, bool delete_data_on_destruct = false);

    virtual ~ho4DArray();

    ho4DArray(const ho4DArray<T>& a);
    ho4DArray<T>& operator=(const ho4DArray<T>& rhs);

    virtual void create(const std::vector<size_t>& dimensions);
    virtual void create(const std::vector<size_t>& dimensions, T* data, bool delete_data_on_destruct = false);

    virtual bool createArray(size_t sx, size_t sy, size_t sz, size_t ss);
    virtual bool createArray(size_t sx, size_t sy, size_t sz, size_t ss, T* data, bool delete_data_on_destruct = false);

    T& operator()(size_t x, size_t y, size_t z, size_t s);
    const T& operator()(size_t x , size_t y, size_t z, size_t s) const;

    virtual void print(std::ostream& os) const;

protected:

    using BaseClass::dimensions_;
    using BaseClass::offsetFactors_;
    using BaseClass::data_;
    using BaseClass::elements_;
    using BaseClass::delete_data_on_destruct_;

    bool init_accesser();
    bool release_accesser();

    T**** accesser_;
};

}

#include "ho4DArray.hxx"
