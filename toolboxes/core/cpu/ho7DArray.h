#pragma once

#include "hoNDArray.h"

namespace Gadgetron{

template <class T> class ho7DArray : public hoNDArray<T>
{
public:

    typedef hoNDArray<T> BaseClass;

    ho7DArray();
    ho7DArray(size_t sx, size_t sy, size_t sz, size_t ss, size_t sp, size_t sr, size_t sa);
    explicit ho7DArray(const std::vector<size_t>& dimensions);
    ho7DArray(const std::vector<size_t>& dimensions, T* data, bool delete_data_on_destruct = false);
    ho7DArray(size_t sx, size_t sy, size_t sz, size_t ss, size_t sp, size_t sr, size_t sa, T* data, bool delete_data_on_destruct = false);
    ho7DArray(boost::shared_ptr< std::vector<size_t> > dimensions);
    ho7DArray(boost::shared_ptr< std::vector<size_t> > dimensions, T* data, bool delete_data_on_destruct = false);

    virtual ~ho7DArray();

    ho7DArray(const ho7DArray<T>& a);
    ho7DArray<T>& operator=(const ho7DArray<T>& rhs);

    virtual void create(const std::vector<size_t>& dimensions);
    virtual void create(const std::vector<size_t>& dimensions, T* data, bool delete_data_on_destruct = false);

    virtual bool createArray(size_t sx, size_t sy, size_t sz, size_t ss, size_t sp, size_t sr, size_t sa);
    virtual bool createArray(size_t sx, size_t sy, size_t sz, size_t ss, size_t sp, size_t sr, size_t sa, T* data, bool delete_data_on_destruct = false);

    T& operator()(size_t x, size_t y, size_t z, size_t s, size_t p, size_t r, size_t a);
    const T& operator()(size_t x , size_t y, size_t z, size_t s, size_t p, size_t r, size_t a) const;

    virtual void print(std::ostream& os) const;

protected:

    using BaseClass::dimensions_;
    using BaseClass::offsetFactors_;
    using BaseClass::data_;
    using BaseClass::elements_;
    using BaseClass::delete_data_on_destruct_;

    bool init_accesser();
    bool release_accesser();

    T******* accesser_;
};

}

#include "ho7DArray.hxx"
