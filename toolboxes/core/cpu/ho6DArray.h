#pragma once

#include "hoNDArray.h"

namespace Gadgetron{

template <class T> class ho6DArray : public hoNDArray<T>
{
public:

    typedef hoNDArray<T> BaseClass;

    ho6DArray();
    ho6DArray(size_t sx, size_t sy, size_t sz, size_t ss, size_t sp, size_t sr);
    explicit ho6DArray(std::vector<size_t> *dimensions);
    ho6DArray(std::vector<size_t> *dimensions, T* data, bool delete_data_on_destruct = false);
    ho6DArray(size_t sx, size_t sy, size_t sz, size_t ss, size_t sp, size_t sr, T* data, bool delete_data_on_destruct = false);
    ho6DArray(boost::shared_ptr< std::vector<size_t> > dimensions);
    ho6DArray(boost::shared_ptr< std::vector<size_t> > dimensions, T* data, bool delete_data_on_destruct = false);

    virtual ~ho6DArray();

    ho6DArray(const ho6DArray<T>& a);
    ho6DArray<T>& operator=(const ho6DArray<T>& rhs);

    virtual void create(std::vector<size_t>& dimensions);
    virtual void create(std::vector<size_t> *dimensions);
    virtual void create(std::vector<size_t> *dimensions, T* data, bool delete_data_on_destruct = false);

    virtual bool createArray(size_t sx, size_t sy, size_t sz, size_t ss, size_t sp, size_t sr);
    virtual bool createArray(size_t sx, size_t sy, size_t sz, size_t ss, size_t sp, size_t sr, T* data, bool delete_data_on_destruct = false);

    T& operator()(size_t x, size_t y, size_t z, size_t s, size_t p, size_t r);
    const T& operator()(size_t x , size_t y, size_t z, size_t s, size_t p, size_t r) const;

    virtual void print(std::ostream& os) const;

protected:

    using BaseClass::dimensions_;
    using BaseClass::offsetFactors_;
    using BaseClass::data_;
    using BaseClass::elements_;
    using BaseClass::delete_data_on_destruct_;

    bool init_accesser();
    bool release_accesser();

    T****** accesser_;
};

}

#include "ho6DArray.hxx"
