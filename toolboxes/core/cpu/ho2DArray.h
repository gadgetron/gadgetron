#pragma once

#include "hoNDArray.h"

namespace Gadgetron{

template <class T> class ho2DArray : public hoNDArray<T>
{
public:

    typedef hoNDArray<T> BaseClass;

    ho2DArray();
    ho2DArray(size_t sx, size_t sy);
    explicit ho2DArray(std::vector<size_t> *dimensions);
    ho2DArray(std::vector<size_t> *dimensions, T* data, bool delete_data_on_destruct = false);
    ho2DArray(size_t sx, size_t sy, T* data, bool delete_data_on_destruct = false);
    ho2DArray(boost::shared_ptr< std::vector<size_t> > dimensions);
    ho2DArray(boost::shared_ptr< std::vector<size_t> > dimensions, T* data, bool delete_data_on_destruct = false);

    virtual ~ho2DArray();

    ho2DArray(const ho2DArray<T>& a);
    ho2DArray<T>& operator=(const ho2DArray<T>& rhs);

    virtual void create(std::vector<size_t>& dimensions);
    virtual void create(std::vector<size_t> *dimensions);
    virtual void create(std::vector<size_t> *dimensions, T* data, bool delete_data_on_destruct = false);

    virtual bool createArray(size_t sx, size_t sy);
    virtual bool createArray(size_t sx, size_t sy, T* data, bool delete_data_on_destruct = false);

    T& operator()(size_t x , size_t y);
    const T& operator()(size_t x , size_t y) const;

    virtual void print(std::ostream& os) const;

protected:

    using BaseClass::dimensions_;
    using BaseClass::offsetFactors_;
    using BaseClass::data_;
    using BaseClass::elements_;
    using BaseClass::delete_data_on_destruct_;

    bool init_accesser();
    bool release_accesser();

    T** accesser_;
};

}

#include <ho2DArray.hxx>
