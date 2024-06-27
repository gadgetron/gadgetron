#pragma once

#include "hoNDArray.h"

namespace Gadgetron{

template <class T> class ho3DArray : public hoNDArray<T>
{
public:

    typedef hoNDArray<T> BaseClass;

    ho3DArray();
    ho3DArray(size_t sx, size_t sy, size_t sz);
    explicit ho3DArray(const std::vector<size_t>& dimensions);
    ho3DArray(const std::vector<size_t>& dimensions, T* data, bool delete_data_on_destruct = false);
    ho3DArray(size_t sx, size_t sy, size_t sz, T* data, bool delete_data_on_destruct = false);
    ho3DArray(boost::shared_ptr< std::vector<size_t> > dimensions);
    ho3DArray(boost::shared_ptr< std::vector<size_t> > dimensions, T* data, bool delete_data_on_destruct = false);

    virtual ~ho3DArray();

    ho3DArray(const ho3DArray<T>& a);
    ho3DArray<T>& operator=(const ho3DArray<T>& rhs);

    virtual void create(const std::vector<size_t>& dimensions);
    virtual void create(const std::vector<size_t>& dimensions, T* data, bool delete_data_on_destruct = false);

    virtual bool createArray(size_t sx, size_t sy, size_t sz);
    virtual bool createArray(size_t sx, size_t sy, size_t sz, T* data, bool delete_data_on_destruct = false);

    T& operator()(size_t x, size_t y, size_t z);
    const T& operator()(size_t x, size_t y, size_t z) const;

    virtual void print(std::ostream& os) const;

protected:

    using BaseClass::dimensions_;
    using BaseClass::offsetFactors_;
    using BaseClass::data_;
    using BaseClass::elements_;
    using BaseClass::delete_data_on_destruct_;

    bool init_accesser();
    bool release_accesser();

    T*** accesser_;
};

}

#include "ho3DArray.hxx"
