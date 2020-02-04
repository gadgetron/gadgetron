#pragma once

#include "hoNDArray.h"

namespace Gadgetron{

template <class T> class ho2DArray : public hoNDArray<T>
{
public:

    typedef hoNDArray<T> BaseClass;

    using BaseClass::create;

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

    T& operator()( const std::vector<size_t>& ind ) { return (*this)(ind[0], ind[1]); }
    const T& operator()( const std::vector<size_t>& ind ) const  { return (*this)(ind[0], ind[1]); }

    T& operator()( size_t x ) { return (*this)(x, 0); }
    const T& operator()( size_t x ) const { return (*this)(x, 0); }

    T& operator()( size_t x, size_t y, size_t z ) { return (*this)(x, y); }
    const T& operator()( size_t x, size_t y, size_t z ) const { return (*this)(x, y); }

    T& operator()( size_t x, size_t y, size_t z, size_t s ) { return (*this)(x, y); }
    const T& operator()( size_t x, size_t y, size_t z, size_t s ) const { return (*this)(x, y); }

    T& operator()( size_t x, size_t y, size_t z, size_t s, size_t p ) { return (*this)(x, y); }
    const T& operator()( size_t x, size_t y, size_t z, size_t s, size_t p ) const { return (*this)(x, y); }

    T& operator()( size_t x, size_t y, size_t z, size_t s, size_t p, size_t r ) { return (*this)(x, y); }
    const T& operator()( size_t x, size_t y, size_t z, size_t s, size_t p, size_t r ) const { return (*this)(x, y); }

    T& operator()( size_t x, size_t y, size_t z, size_t s, size_t p, size_t r, size_t a ) { return (*this)(x, y); }
    const T& operator()( size_t x, size_t y, size_t z, size_t s, size_t p, size_t r, size_t a ) const { return (*this)(x, y); }

    T& operator()( size_t x, size_t y, size_t z, size_t s, size_t p, size_t r, size_t a, size_t q ) { return (*this)(x, y); }
    const T& operator()( size_t x, size_t y, size_t z, size_t s, size_t p, size_t r, size_t a, size_t q ) const { return (*this)(x, y); }

    T& operator()( size_t x, size_t y, size_t z, size_t s, size_t p, size_t r, size_t a, size_t q, size_t u ) { return (*this)(x, y); }
    const T& operator()( size_t x, size_t y, size_t z, size_t s, size_t p, size_t r, size_t a, size_t q, size_t u ) const { return (*this)(x, y); }

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

#include "ho2DArray.hxx"
