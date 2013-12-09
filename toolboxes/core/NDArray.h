/** \file NDArray.h
\brief Abstract base class for all Gadgetron host and device arrays
*/

#ifndef NDARRAY_H
#define NDARRAY_H
#pragma once

#include <new>
#include <vector>
#include <iostream>
#include <stdexcept>
#include <boost/shared_ptr.hpp>
#include <boost/cast.hpp>
#include <new>          // std::runtime_error
#include "GadgetronException.h"
#include "GadgetronCommon.h"

namespace Gadgetron{

template <typename T> class NDArray
{
public:

    typedef T element_type;
    typedef T value_type;

    void* operator new (size_t bytes)
    {
        return ::new char[bytes];
    }

    void operator delete (void *ptr)
    {
        delete [] static_cast <char *> (ptr);
    } 

    void * operator new(size_t s, void * p)
    {
        return p;
    }

    NDArray () : data_(0), elements_(0), delete_data_on_destruct_(true)
    {
        dimensions_ = boost::shared_ptr< std::vector<unsigned long long> >( new std::vector<unsigned long long> );
    }

    virtual ~NDArray() {}

    virtual void create(std::vector<unsigned long long>& dimensions);

    virtual void create(std::vector<unsigned int> *dimensions);
    virtual void create(std::vector<unsigned long long> *dimensions);

    virtual void create(std::vector<unsigned int> *dimensions, T* data, bool delete_data_on_destruct = false);
    virtual void create(std::vector<unsigned long long> *dimensions, T* data, bool delete_data_on_destruct = false);

    virtual void create(boost::shared_ptr< std::vector<unsigned long long> > dimensions);
    virtual void create(boost::shared_ptr<std::vector<unsigned long long>  > dimensions, T* data, bool delete_data_on_destruct = false);

    /*virtual void create(unsigned long long len);
    virtual void create(unsigned long long sx, unsigned long long sy);
    virtual void create(unsigned long long sx, unsigned long long sy, unsigned long long sz);
    virtual void create(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long st);
    virtual void create(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long st, unsigned long long sp);
    virtual void create(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long st, unsigned long long sp, unsigned long long sq);
    virtual void create(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long st, unsigned long long sp, unsigned long long sq, unsigned long long sr);
    virtual void create(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long st, unsigned long long sp, unsigned long long sq, unsigned long long sr, unsigned long long ss);
    virtual void create(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long st, unsigned long long sp, unsigned long long sq, unsigned long long sr, unsigned long long ss, unsigned long long su);

    virtual void create(unsigned long long len, T* data, bool delete_data_on_destruct = false);
    virtual void create(unsigned long long sx, unsigned long long sy, T* data, bool delete_data_on_destruct = false);
    virtual void create(unsigned long long sx, unsigned long long sy, unsigned long long sz, T* data, bool delete_data_on_destruct = false);
    virtual void create(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long st, T* data, bool delete_data_on_destruct = false);
    virtual void create(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long st, unsigned long long sp, T* data, bool delete_data_on_destruct = false);
    virtual void create(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long st, unsigned long long sp, unsigned long long sq, T* data, bool delete_data_on_destruct = false);
    virtual void create(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long st, unsigned long long sp, unsigned long long sq, unsigned long long sr, T* data, bool delete_data_on_destruct = false);
    virtual void create(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long st, unsigned long long sp, unsigned long long sq, unsigned long long sr, unsigned long long ss, T* data, bool delete_data_on_destruct = false);
    virtual void create(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long st, unsigned long long sp, unsigned long long sq, unsigned long long sr, unsigned long long ss, unsigned long long su, T* data, bool delete_data_on_destruct = false);*/

    void squeeze();
    void reshape(std::vector<unsigned long long> *dims);
    void reshape( boost::shared_ptr< std::vector<unsigned long long> > dims );
    bool dimensions_equal(std::vector<unsigned long long> *d) const;

    template<class S> bool dimensions_equal(const NDArray<S> *a) const
    {
        boost::shared_ptr<std::vector<unsigned long long > > adims = a->get_dimensions();
        return ((this->dimensions_->size() == adims->size()) &&
            std::equal(this->dimensions_->begin(), this->dimensions_->end(), adims->begin()));
    }

    unsigned long long get_number_of_dimensions() const;
    unsigned long long get_size(unsigned long long dimension) const;
    boost::shared_ptr< std::vector<unsigned long long> > get_dimensions() const;

    T* get_data_ptr() const;

    unsigned long long get_number_of_elements() const;
    unsigned long long get_number_of_bytes() const;

    bool delete_data_on_destruct() const;
    void delete_data_on_destruct(bool d);

    long long calculate_offset(const std::vector<unsigned long long>& ind) const;
    size_t get_offset_factor(unsigned long long dim) const;
    size_t get_offset_factor_lastdim() const;

    boost::shared_ptr< std::vector<size_t> > get_offset_factor() const;
    void calculate_offset_factors(const std::vector<unsigned long long>& dimensions);

    std::vector<unsigned long long> calculate_index( unsigned long long offset ) const;
    void calculate_index( unsigned long long offset, std::vector<unsigned long long>& index ) const;

    void clear();

protected:

    virtual void allocate_memory() = 0;
    virtual void deallocate_memory() = 0;

protected:

    boost::shared_ptr< std::vector<unsigned long long> > dimensions_;
    boost::shared_ptr< std::vector<size_t> > offsetFactors_;
    T* data_;
    unsigned long long elements_; // to handle the big array
    bool delete_data_on_destruct_;
};

template <typename T> 
inline void NDArray<T>::create(std::vector<unsigned long long>& dimensions) 
{
    std::vector<unsigned long long> *tmp = new std::vector<unsigned long long>;
    *tmp = dimensions;
    dimensions_ = boost::shared_ptr< std::vector<unsigned long long> >(tmp);
    allocate_memory();
    calculate_offset_factors(*dimensions_);
}

template <typename T> 
inline void NDArray<T>::create(std::vector<unsigned int> *dimensions) 
{
    std::vector<unsigned long long> dims(dimensions->size(), 0);
    for ( unsigned int ii=0; ii<dimensions->size(); ii++ )
    {
        dims[ii] = (*dimensions)[ii];
    }

    this->create(&dims);
}

template <typename T> 
inline void NDArray<T>::create(std::vector<unsigned long long> *dimensions) 
{
    std::vector<unsigned long long> *tmp = new std::vector<unsigned long long>;
    *tmp = *dimensions;
    dimensions_ = boost::shared_ptr< std::vector<unsigned long long> >(tmp);
    allocate_memory();
    calculate_offset_factors(*dimensions_);
}

template <typename T> 
inline void NDArray<T>::create(std::vector<unsigned int> *dimensions, T* data, bool delete_data_on_destruct) 
{
    std::vector<unsigned long long> dims(dimensions->size(), 0);
    for ( unsigned int ii=0; ii<dimensions->size(); ii++ )
    {
        dims[ii] = (*dimensions)[ii];
    }

    this->create(&dims, data, delete_data_on_destruct);
}

template <typename T> 
void NDArray<T>::create(std::vector<unsigned long long> *dimensions, T* data, bool delete_data_on_destruct) 
{
    if (!data)
    {
        GADGET_THROW("NDArray<T>::create: 0x0 pointer provided");
    }

    std::vector<unsigned long long> *tmp = new std::vector<unsigned long long>;
    *tmp = *dimensions;
    dimensions_ = boost::shared_ptr< std::vector<unsigned long long> >(tmp);
    this->data_ = data;
    this->delete_data_on_destruct_ = delete_data_on_destruct;
    this->elements_ = 1;

    for (unsigned long long i = 0; i < this->dimensions_->size(); i++)
    {
        this->elements_ *= (*this->dimensions_)[i];
    }
    calculate_offset_factors(*dimensions_);
}

template <typename T> 
inline void NDArray<T>::create(boost::shared_ptr< std::vector<unsigned long long> > dimensions)
{
    this->create(dimensions.get());
}

template <typename T> 
inline void NDArray<T>::create(boost::shared_ptr<std::vector<unsigned long long>  > dimensions, 
    T* data, bool delete_data_on_destruct)
{
    this->create(dimensions.get(), data, delete_data_on_destruct);
}

//template <typename T> 
//inline void NDArray<T>::create(unsigned long long len)
//{
//    std::vector<unsigned long long> dim(1);
//    dim[0] = len;
//    this->create(dim);
//}
//
//template <typename T> 
//inline void NDArray<T>::create(unsigned long long sx, unsigned long long sy)
//{
//    std::vector<unsigned long long> dim(2);
//    dim[0] = sx;
//    dim[1] = sy;
//    this->create(dim);
//}
//
//template <typename T> 
//inline void NDArray<T>::create(unsigned long long sx, unsigned long long sy, unsigned long long sz)
//{
//    std::vector<unsigned long long> dim(3);
//    dim[0] = sx;
//    dim[1] = sy;
//    dim[2] = sz;
//    this->create(dim);
//}
//
//template <typename T> 
//inline void NDArray<T>::create(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long st)
//{
//    std::vector<unsigned long long> dim(4);
//    dim[0] = sx;
//    dim[1] = sy;
//    dim[2] = sz;
//    dim[3] = st;
//    this->create(dim);
//}
//
//template <typename T> 
//inline void NDArray<T>::create(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long st, unsigned long long sp)
//{
//    std::vector<unsigned long long> dim(5);
//    dim[0] = sx;
//    dim[1] = sy;
//    dim[2] = sz;
//    dim[3] = st;
//    dim[4] = sp;
//    this->create(dim);
//}
//
//template <typename T> 
//inline void NDArray<T>::create(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long st, unsigned long long sp, unsigned long long sq)
//{
//    std::vector<unsigned long long> dim(6);
//    dim[0] = sx;
//    dim[1] = sy;
//    dim[2] = sz;
//    dim[3] = st;
//    dim[4] = sp;
//    dim[5] = sq;
//    this->create(dim);
//}
//
//template <typename T> 
//inline void NDArray<T>::create(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long st, unsigned long long sp, unsigned long long sq, unsigned long long sr)
//{
//    std::vector<unsigned long long> dim(7);
//    dim[0] = sx;
//    dim[1] = sy;
//    dim[2] = sz;
//    dim[3] = st;
//    dim[4] = sp;
//    dim[5] = sq;
//    dim[6] = sr;
//    this->create(dim);
//}
//
//template <typename T> 
//inline void NDArray<T>::create(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long st, unsigned long long sp, unsigned long long sq, unsigned long long sr, unsigned long long ss)
//{
//    std::vector<unsigned long long> dim(8);
//    dim[0] = sx;
//    dim[1] = sy;
//    dim[2] = sz;
//    dim[3] = st;
//    dim[4] = sp;
//    dim[5] = sq;
//    dim[6] = sr;
//    dim[7] = ss;
//    this->create(dim);
//}
//
//template <typename T> 
//inline void NDArray<T>::create(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long st, unsigned long long sp, unsigned long long sq, unsigned long long sr, unsigned long long ss, unsigned long long su)
//{
//    std::vector<unsigned long long> dim(9);
//    dim[0] = sx;
//    dim[1] = sy;
//    dim[2] = sz;
//    dim[3] = st;
//    dim[4] = sp;
//    dim[5] = sq;
//    dim[6] = sr;
//    dim[7] = ss;
//    dim[8] = su;
//    this->create(dim);
//}
//
//template <typename T> 
//inline void NDArray<T>::create(unsigned long long len, T* data, bool delete_data_on_destruct)
//{
//    std::vector<unsigned long long> dim(1);
//    dim[0] = len;
//    this->create(&dim, data, delete_data_on_destruct);
//}
//
//template <typename T> 
//inline void NDArray<T>::create(unsigned long long sx, unsigned long long sy, T* data, bool delete_data_on_destruct)
//{
//    std::vector<unsigned long long> dim(2);
//    dim[0] = sx;
//    dim[1] = sy;
//    this->create(&dim, data, delete_data_on_destruct);
//}
//
//template <typename T> 
//inline void NDArray<T>::create(unsigned long long sx, unsigned long long sy, unsigned long long sz, T* data, bool delete_data_on_destruct)
//{
//    std::vector<unsigned long long> dim(3);
//    dim[0] = sx;
//    dim[1] = sy;
//    dim[2] = sz;
//    this->create(&dim, data, delete_data_on_destruct);
//}
//
//template <typename T> 
//inline void NDArray<T>::create(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long st, T* data, bool delete_data_on_destruct)
//{
//    std::vector<unsigned long long> dim(4);
//    dim[0] = sx;
//    dim[1] = sy;
//    dim[2] = sz;
//    dim[3] = st;
//    this->create(&dim, data, delete_data_on_destruct);
//}
//
//template <typename T> 
//inline void NDArray<T>::create(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long st, unsigned long long sp, T* data, bool delete_data_on_destruct)
//{
//    std::vector<unsigned long long> dim(5);
//    dim[0] = sx;
//    dim[1] = sy;
//    dim[2] = sz;
//    dim[3] = st;
//    dim[4] = sp;
//    this->create(&dim, data, delete_data_on_destruct);
//}
//
//template <typename T> 
//inline void NDArray<T>::create(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long st, unsigned long long sp, unsigned long long sq, T* data, bool delete_data_on_destruct)
//{
//    std::vector<unsigned long long> dim(6);
//    dim[0] = sx;
//    dim[1] = sy;
//    dim[2] = sz;
//    dim[3] = st;
//    dim[4] = sp;
//    dim[5] = sq;
//    this->create(&dim, data, delete_data_on_destruct);
//}
//
//template <typename T> 
//inline void NDArray<T>::create(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long st, unsigned long long sp, unsigned long long sq, unsigned long long sr, T* data, bool delete_data_on_destruct)
//{
//    std::vector<unsigned long long> dim(7);
//    dim[0] = sx;
//    dim[1] = sy;
//    dim[2] = sz;
//    dim[3] = st;
//    dim[4] = sp;
//    dim[5] = sq;
//    dim[6] = sr;
//    this->create(&dim, data, delete_data_on_destruct);
//}
//
//template <typename T> 
//inline void NDArray<T>::create(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long st, unsigned long long sp, unsigned long long sq, unsigned long long sr, unsigned long long ss, T* data, bool delete_data_on_destruct)
//{
//    std::vector<unsigned long long> dim(8);
//    dim[0] = sx;
//    dim[1] = sy;
//    dim[2] = sz;
//    dim[3] = st;
//    dim[4] = sp;
//    dim[5] = sq;
//    dim[6] = sr;
//    dim[7] = ss;
//    this->create(&dim, data, delete_data_on_destruct);
//}
//
//template <typename T> 
//inline void NDArray<T>::create(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long st, unsigned long long sp, unsigned long long sq, unsigned long long sr, unsigned long long ss, unsigned long long su, T* data, bool delete_data_on_destruct)
//{
//    std::vector<unsigned long long> dim(9);
//    dim[0] = sx;
//    dim[1] = sy;
//    dim[2] = sz;
//    dim[3] = st;
//    dim[4] = sp;
//    dim[5] = sq;
//    dim[6] = sr;
//    dim[7] = ss;
//    dim[8] = su;
//    this->create(&dim, data, delete_data_on_destruct);
//}

template <typename T> 
inline void NDArray<T>::squeeze()
{
    boost::shared_ptr< std::vector<unsigned long long> > new_dimensions( new std::vector<unsigned long long> ); 
    for (unsigned long long i = 0; i < dimensions_->size(); i++)
    {
        if ((*dimensions_)[i] != 1)
        {
            new_dimensions->push_back((*dimensions_)[i]);
        }
    }

    dimensions_ = new_dimensions;
    this->calculate_offset_factors(*dimensions_);
}

template <typename T> 
inline void NDArray<T>::reshape(std::vector<unsigned long long> *dims)
{
    unsigned long long new_elements = 1;
    for (unsigned long long i = 0; i < dims->size(); i++)
    {
        new_elements *= (*dims)[i];
    }

    if (new_elements != elements_)
    {
        std::runtime_error("NDArray<T>::reshape : Number of elements cannot change during reshape");
    }

    // Copy the input dimensions array
    std::vector<unsigned long long> *tmp = new std::vector<unsigned long long>;
    *tmp = *dims;
    dimensions_ = boost::shared_ptr< std::vector<unsigned long long> >(tmp);
    this->calculate_offset_factors(*dimensions_);
}

template <typename T> 
inline void NDArray<T>::reshape( boost::shared_ptr< std::vector<unsigned long long> > dims )
{
    this->reshape(dims.get());
}

template <typename T> 
inline bool NDArray<T>::dimensions_equal(std::vector<unsigned long long> *d) const
{
    return ((this->dimensions_->size() == d->size()) &&
        std::equal(this->dimensions_->begin(), this->dimensions_->end(), d->begin()));
}

template <typename T> 
inline unsigned long long NDArray<T>::get_number_of_dimensions() const
{
    return (unsigned long long)dimensions_->size();
}

template <typename T> 
inline unsigned long long NDArray<T>::get_size(unsigned long long dimension) const
{
    if (dimension >= dimensions_->size())
    {
        return 1;
    }
    else
    {
        return (*dimensions_)[dimension];
    }
}

template <typename T> 
inline boost::shared_ptr< std::vector<unsigned long long> > NDArray<T>::get_dimensions() const
{
    // Make copy to ensure that the receiver cannot alter the array dimensions
    std::vector<unsigned long long> *tmp = new std::vector<unsigned long long>;
    *tmp=*dimensions_;
    return boost::shared_ptr< std::vector<unsigned long long> >(tmp); 
}

template <typename T> 
inline T* NDArray<T>::get_data_ptr() const
{ 
    return data_;
}

template <typename T> 
inline unsigned long long NDArray<T>::get_number_of_elements() const
{
    return elements_;
}

template <typename T> 
inline unsigned long long NDArray<T>::get_number_of_bytes() const
{
    return elements_*sizeof(T);
}

template <typename T> 
inline bool NDArray<T>::delete_data_on_destruct() const
{
    return delete_data_on_destruct_;
}

template <typename T> 
inline void NDArray<T>::delete_data_on_destruct(bool d)
{
    delete_data_on_destruct_ = d;
}

template <typename T> 
inline long long NDArray<T>::calculate_offset(const std::vector<unsigned long long>& ind) const
{
    long long offset = ind[0];
    unsigned long long i;
    for( i = 1; i < dimensions_->size(); i++ )
        offset += ind[i] * (*offsetFactors_)[i];
    return offset;
}

template <typename T> 
inline size_t NDArray<T>::get_offset_factor(unsigned long long dim) const
{
    if ( dim >= (*dimensions_).size() )
    {
        std::runtime_error("NDArray<T>::reshape : Number of elements cannot change during reshape");
        return 0;
    }

    return (*offsetFactors_)[dim];
}

template <typename T> 
inline size_t NDArray<T>::get_offset_factor_lastdim() const
{
    return get_offset_factor(dimensions_->size()-1);
}

template <typename T> 
inline boost::shared_ptr< std::vector<size_t> > NDArray<T>::get_offset_factor() const
{
    std::vector<size_t> *tmp = new std::vector<size_t>;
    *tmp=*offsetFactors_;
    return boost::shared_ptr< std::vector<size_t> >(tmp); 
}

template <typename T> 
inline void NDArray<T>::calculate_offset_factors(const std::vector<unsigned long long>& dimensions)
{
    if ( offsetFactors_.get() == NULL )
    {
        std::vector<size_t> *tmp = new std::vector<size_t>;
        offsetFactors_ = boost::shared_ptr< std::vector<size_t> >(tmp);
    }
    offsetFactors_->resize(dimensions.size());
    for( unsigned long long i = 0; i < dimensions.size(); i++ )
    {
        unsigned long long k = 1;
        for( unsigned long long j = 0; j < i; j++ )
            k *= dimensions[j];
        (*offsetFactors_)[i] = k;
    }
}

template <typename T> 
inline std::vector<unsigned long long> NDArray<T>::calculate_index( unsigned long long offset ) const
{
    std::vector<unsigned long long> index(dimensions_->size());
    for( unsigned long long i = dimensions_->size()-1; i>=0; i-- )
    {
        index[i] = offset / (*offsetFactors_)[i];
        offset %= (*offsetFactors_)[i];
    }
    return index;
}

template <typename T> 
inline void NDArray<T>::calculate_index( unsigned long long offset, std::vector<unsigned long long>& index ) const
{
    index.resize(dimensions_->size(), 0);
    for( unsigned long long i = dimensions_->size()-1; i>=0; i-- )
    {
        index[i] = offset / (*offsetFactors_)[i];
        offset %= (*offsetFactors_)[i];
    }
    return;
}

template <typename T> 
void NDArray<T>::clear()
{
    if ( this->delete_data_on_destruct_ )
    {
        this->deallocate_memory();
    }

    this->data_ = 0;
    this->elements_ = 0;
    this->delete_data_on_destruct_ = true;

    if ( !this->dimensions_ )
    {
        this->dimensions_->clear();
        this->offsetFactors_->clear();
    }
    else
    {
        this->dimensions_ = boost::shared_ptr< std::vector<unsigned long long> >( new std::vector<unsigned long long> );
        this->offsetFactors_ = boost::shared_ptr< std::vector<size_t> >( new std::vector<size_t> );
    }
}

}

#endif //NDARRAY_H
