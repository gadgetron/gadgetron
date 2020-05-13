/** \file NDArray.h
\brief Abstract base class for all Gadgetron host and device arrays
*/

#pragma once

#include "GadgetronException.h"
#include "log.h"

#include <new>
#include <vector>
#include <iostream>
#include <stdexcept>
#include <array>
#include <algorithm>

#include <boost/cast.hpp>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>
#include <numeric>

namespace Gadgetron{

    template <typename T> class NDArray
    {
    public:

        typedef T element_type;
        typedef T value_type;

        NDArray () : data_(0), elements_(0), delete_data_on_destruct_(true)
        {
        }

        virtual ~NDArray() {}

        virtual void create(const std::vector<size_t> &dimensions);
        virtual void create(const std::vector<size_t> *dimensions);
        virtual void create(boost::shared_ptr< std::vector<size_t> > dimensions);

        virtual void create(const std::vector<size_t> &dimensions, T* data, bool delete_data_on_destruct = false);
        virtual void create(const std::vector<size_t> *dimensions, T* data, bool delete_data_on_destruct = false);
        virtual void create(boost::shared_ptr< std::vector<size_t> > dimensions, T* data, bool delete_data_on_destruct = false);

        void squeeze();

        void reshape(const std::vector<size_t> *dims);
        void reshape(const std::vector<size_t> & dims){ this->reshape(&dims);}
        void reshape(boost::shared_ptr< std::vector<size_t> > dims);

        /**
         * Reshapes the array to the given dimensions.
         * One of the dimensions can be -1, in which case that dimension will be calculated based on the other.
         * @param dims
         */
        void reshape(std::initializer_list<ssize_t> dims);

        bool dimensions_equal(const std::vector<size_t> *d) const;
        bool dimensions_equal(const std::vector<size_t>& d) const;

        template<class S> bool dimensions_equal(const NDArray<S> *a) const
        {
            std::vector<size_t> dim;
            a->get_dimensions(dim);

            if ( this->dimensions_.size() != dim.size() ) return false;

            size_t NDim = this->dimensions_.size();
            for ( size_t d=0; d<NDim; d++ )
            {
                if ( this->dimensions_[d] != dim[d] ) return false;
            }

            return true;
        }

        size_t get_number_of_dimensions() const;

        size_t get_size(size_t dimension) const;

        boost::shared_ptr< std::vector<size_t> > get_dimensions() const;
        void get_dimensions(std::vector<size_t>& dim) const;

        std::vector<size_t> const &dimensions() const;

        const T* get_data_ptr() const;
        T* get_data_ptr();

        const T* data() const;
        T* data();

        size_t size() const;
        size_t get_number_of_elements() const;

        bool empty() const;

        size_t get_number_of_bytes() const;

        bool delete_data_on_destruct() const;
        void delete_data_on_destruct(bool d);

        size_t calculate_offset(const std::vector<size_t>& ind) const;
        static size_t calculate_offset(const std::vector<size_t>& ind, const std::vector<size_t>& offsetFactors);

        size_t calculate_offset(size_t x, size_t y) const;
        size_t calculate_offset(size_t x, size_t y, size_t z) const;
        size_t calculate_offset(size_t x, size_t y, size_t z, size_t s) const;
        size_t calculate_offset(size_t x, size_t y, size_t z, size_t s, size_t p) const;
        size_t calculate_offset(size_t x, size_t y, size_t z, size_t s, size_t p, size_t r) const;
        size_t calculate_offset(size_t x, size_t y, size_t z, size_t s, size_t p, size_t r, size_t a) const;
        size_t calculate_offset(size_t x, size_t y, size_t z, size_t s, size_t p, size_t r, size_t a, size_t q) const;
        size_t calculate_offset(size_t x, size_t y, size_t z, size_t s, size_t p, size_t r, size_t a, size_t q, size_t u) const;

        size_t get_offset_factor(size_t dim) const;
        void get_offset_factor(std::vector<size_t>& offset) const;
        boost::shared_ptr< std::vector<size_t> > get_offset_factor() const;

        size_t get_offset_factor_lastdim() const;

        void calculate_offset_factors(const std::vector<size_t>& dimensions);
        static void calculate_offset_factors(const std::vector<size_t>& dimensions, std::vector<size_t>& offsetFactors);

        std::vector<size_t> calculate_index( size_t offset ) const;
        void calculate_index( size_t offset, std::vector<size_t>& index ) const;
        static void calculate_index( size_t offset, const std::vector<size_t>& offsetFactors, std::vector<size_t>& index );

        void clear();

    

        /// whether a point is within the array range
        bool point_in_range(const std::vector<size_t>& ind) const;
        bool point_in_range(size_t x) const;
        bool point_in_range(size_t x, size_t y) const;
        bool point_in_range(size_t x, size_t y, size_t z) const;
        bool point_in_range(size_t x, size_t y, size_t z, size_t s) const;
        bool point_in_range(size_t x, size_t y, size_t z, size_t s, size_t p) const;
        bool point_in_range(size_t x, size_t y, size_t z, size_t s, size_t p, size_t r) const;
        bool point_in_range(size_t x, size_t y, size_t z, size_t s, size_t p, size_t r, size_t a) const;
        bool point_in_range(size_t x, size_t y, size_t z, size_t s, size_t p, size_t r, size_t a, size_t q) const;
        bool point_in_range(size_t x, size_t y, size_t z, size_t s, size_t p, size_t r, size_t a, size_t q, size_t u) const;

    protected:

        virtual void allocate_memory() = 0;
        virtual void deallocate_memory() = 0;

    protected:

        std::vector<size_t> dimensions_;
        std::array<size_t, 12> offsetFactors_;
        T* data_;
        size_t elements_;
        bool delete_data_on_destruct_;
    };

    template <typename T> 
    inline void NDArray<T>::create(const std::vector<size_t> *dimensions)
    {
        if(!dimensions) throw std::runtime_error("NDArray<T>::create(): 0x0 pointer provided");
        dimensions_ = *dimensions;
        allocate_memory();
        calculate_offset_factors(dimensions_);
    }

    template <typename T> 
    inline void NDArray<T>::create(const std::vector<size_t>& dimensions)
    {
        dimensions_ = dimensions;
        allocate_memory();
        calculate_offset_factors(dimensions_);
    }

    template <typename T> 
    inline void NDArray<T>::create(boost::shared_ptr< std::vector<size_t> > dimensions)
    {
        this->create(dimensions.get());
    }

    template <typename T> 
    void NDArray<T>::create(const std::vector<size_t> *dimensions, T* data, bool delete_data_on_destruct)
    {
        if (!dimensions) throw std::runtime_error("NDArray<T>::create(): 0x0 pointer provided");
        if (!data) throw std::runtime_error("NDArray<T>::create(): 0x0 pointer provided");    
        dimensions_ = *dimensions;
        this->data_ = data;
        this->delete_data_on_destruct_ = delete_data_on_destruct;
        this->elements_ = 1;
        for (size_t i = 0; i < this->dimensions_.size(); i++){
            this->elements_ *= this->dimensions_[i];
        }
        calculate_offset_factors(dimensions_);
    }

    template <typename T> 
    void NDArray<T>::create(const std::vector<size_t> &dimensions, T* data, bool delete_data_on_destruct)
    {
        if (!data) throw std::runtime_error("NDArray<T>::create(): 0x0 pointer provided");    
        dimensions_ = dimensions;
        this->data_ = data;
        this->delete_data_on_destruct_ = delete_data_on_destruct;
        this->elements_ = 1;
        for (size_t i = 0; i < this->dimensions_.size(); i++){
            this->elements_ *= this->dimensions_[i];
        }
        calculate_offset_factors(dimensions_);
    }

    template <typename T> 
    inline void NDArray<T>::create(boost::shared_ptr<std::vector<size_t>  > dimensions, 
        T* data, bool delete_data_on_destruct)
    {
        this->create(dimensions.get(), data, delete_data_on_destruct);
    }

    template <typename T> 
    inline void NDArray<T>::squeeze()
    {
        std::vector<size_t> new_dimensions;
        for (size_t i = 0; i < dimensions_.size(); i++){
            if (dimensions_[i] != 1){
                new_dimensions.push_back(dimensions_[i]);
            }
        }
        dimensions_ = new_dimensions;
        this->calculate_offset_factors(dimensions_);
    }

    template <typename T> 
    inline void NDArray<T>::reshape(const std::vector<size_t> *dims)
    {
        size_t new_elements = 1;
        for (size_t i = 0; i < dims->size(); i++){
            new_elements *= (*dims)[i];
        }

        if (new_elements != elements_)
            throw std::runtime_error("NDArray<T>::reshape : Number of elements cannot change during reshape");    

        // Copy the input dimensions array
        dimensions_ = *dims;
        this->calculate_offset_factors(dimensions_);
    }

    template <typename T> void NDArray<T>::reshape(std::initializer_list<ssize_t> dims) {
        std::vector<ssize_t> dim_vec(dims);
        auto negatives = std::count(dims.begin(), dims.end(), -1);
        if (negatives > 1)
            throw std::runtime_error("Only a single reshape dimension can be negative");

        if (negatives == 1) {
            auto elements    = std::accumulate(dims.begin(), dims.end(), -1, std::multiplies<ssize_t>());
            auto neg_element = std::find(dim_vec.begin(), dim_vec.end(), -1);
            *neg_element     = this->elements_ / elements;
        }

        auto new_dims = std::vector<size_t>(dim_vec.begin(), dim_vec.end());
        this->reshape(new_dims);
    }
    template <typename T> 
    inline void NDArray<T>::reshape( boost::shared_ptr< std::vector<size_t> > dims )
    {
        this->reshape(dims.get());
    }

    template <typename T> 
    inline bool NDArray<T>::dimensions_equal(const std::vector<size_t>& d) const
    {
        if ( this->dimensions_.size() != d.size() ) return false;

        size_t NDim = this->dimensions_.size();
        for ( size_t ii=0; ii<NDim; ii++ )
        {
            if ( this->dimensions_[ii] != d[ii] ) return false;
        }

        return true;
    }

    template <typename T>
    inline bool NDArray<T>::dimensions_equal(const std::vector<size_t>* d) const
    {
        return this->dimensions_equal(*d);
    }
    template <typename T>
    inline size_t NDArray<T>::get_number_of_dimensions() const
    {
        return (size_t)dimensions_.size();
    }

    template <typename T> 
    inline size_t NDArray<T>::get_size(size_t dimension) const
    {
        if (dimension >= dimensions_.size()){
            return 1;
        }
        else{
            return dimensions_[dimension];
        }
    }

    template <typename T> 
    inline boost::shared_ptr< std::vector<size_t> > NDArray<T>::get_dimensions() const
    {
        // Make copy to ensure that the receiver cannot alter the array dimensions
        std::vector<size_t> *tmp = new std::vector<size_t>;
        *tmp=dimensions_;
        return boost::shared_ptr< std::vector<size_t> >(tmp); 
    }

    template <typename T> 
    inline void NDArray<T>::get_dimensions(std::vector<size_t>& dim) const
    {
        dim = dimensions_;
    }

    template<class T>
    inline std::vector<size_t> const &NDArray<T>::dimensions() const {
        return dimensions_;
    }


    template <typename T> 
    inline const T* NDArray<T>::get_data_ptr() const
    { 
        return data_;
    }
    template <typename T>
    inline T* NDArray<T>::get_data_ptr()
    {
        return data_;
    }

    template <typename T>
    const T* NDArray<T>::data() const
    {
        return data_;
    }

    template <typename T>
    T* NDArray<T>::data()
    {
        return data_;
    }


    template<class T>
    inline size_t NDArray<T>::size() const
    {
        return elements_;
    }

    template<class T>
    inline bool NDArray<T>::empty() const
    {
        return elements_ == 0;
    }

    template <class T>
    inline size_t NDArray<T>::get_number_of_elements() const {
        return size();
    }

    template <typename T> 
    inline size_t NDArray<T>::get_number_of_bytes() const
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
    size_t NDArray<T>::calculate_offset(const std::vector<size_t>& ind, const std::vector<size_t>& offsetFactors)
    {
        size_t offset = ind[0];

        for( size_t i = 1; i < ind.size(); i++ )
        {
            offset += ind[i] * offsetFactors[i];
        }

        return offset;
    }

    template <typename T> 
    inline size_t NDArray<T>::calculate_offset(const std::vector<size_t>& ind) const
    {
        size_t offset = ind[0];
        for( size_t i = 1; i < dimensions_.size(); i++ )
            offset += ind[i] * offsetFactors_[i];
        return offset;
    }

    template <typename T> 
    inline size_t NDArray<T>::calculate_offset(size_t x, size_t y) const
    {
//        GADGET_DEBUG_CHECK_THROW(dimensions_.size()==2);
        return x + y * offsetFactors_[1];
    }

    template <typename T> 
    inline size_t NDArray<T>::calculate_offset(size_t x, size_t y, size_t z) const
    {
//        GADGET_DEBUG_CHECK_THROW(dimensions_.size()==3);
        return x + y * offsetFactors_[1] + z * offsetFactors_[2];
    }

    template <typename T> 
    inline size_t NDArray<T>::calculate_offset(size_t x, size_t y, size_t z, size_t s) const
    {
//        GADGET_DEBUG_CHECK_THROW(dimensions_.size()==4);
        return x + y * offsetFactors_[1] + z * offsetFactors_[2] + s * offsetFactors_[3];
    }

    template <typename T> 
    inline size_t NDArray<T>::calculate_offset(size_t x, size_t y, size_t z, size_t s, size_t p) const
    {
//        GADGET_DEBUG_CHECK_THROW(dimensions_.size()==5);
        return x + y * offsetFactors_[1] + z * offsetFactors_[2] + s * offsetFactors_[3] + p * offsetFactors_[4];
    }

    template <typename T> 
    inline size_t NDArray<T>::calculate_offset(size_t x, size_t y, size_t z, size_t s, size_t p, size_t r) const
    {
//        GADGET_DEBUG_CHECK_THROW(dimensions_.size()==6);
        return x + y * offsetFactors_[1] + z * offsetFactors_[2] + s * offsetFactors_[3] + p * offsetFactors_[4] + r * offsetFactors_[5];
    }

    template <typename T> 
    inline size_t NDArray<T>::calculate_offset(size_t x, size_t y, size_t z, size_t s, size_t p, size_t r, size_t a) const
    {
//        GADGET_DEBUG_CHECK_THROW(dimensions_.size()==7);
        return x + y * offsetFactors_[1] + z * offsetFactors_[2] + s * offsetFactors_[3] + p * offsetFactors_[4] + r * offsetFactors_[5] + a * offsetFactors_[6];
    }

    template <typename T> 
    inline size_t NDArray<T>::calculate_offset(size_t x, size_t y, size_t z, size_t s, size_t p, size_t r, size_t a, size_t q) const
    {
//        GADGET_DEBUG_CHECK_THROW(dimensions_.size()==8);
        return x + y * offsetFactors_[1] + z * offsetFactors_[2] + s * offsetFactors_[3] + p * offsetFactors_[4] + r * offsetFactors_[5] + a * offsetFactors_[6] + q * offsetFactors_[7];
    }

    template <typename T> 
    inline size_t NDArray<T>::calculate_offset(size_t x, size_t y, size_t z, size_t s, size_t p, size_t r, size_t a, size_t q, size_t u) const
    {
//        GADGET_DEBUG_CHECK_THROW(dimensions_.size()==9);
        return x + y * offsetFactors_[1] + z * offsetFactors_[2] + s * offsetFactors_[3] + p * offsetFactors_[4] + r * offsetFactors_[5] + a * offsetFactors_[6] + q * offsetFactors_[7]+ u * offsetFactors_[8];
    }

    template <typename T> 
    inline size_t NDArray<T>::get_offset_factor(size_t dim) const
    {
        if ( dim >= dimensions_.size() )
            throw std::runtime_error("NDArray<T>::get_offset_factor : index out of range");
        return offsetFactors_[dim];
    }

    template <typename T> 
    inline void NDArray<T>::get_offset_factor(std::vector<size_t>& offset) const
    {
        offset=std::vector<size_t>(offsetFactors_.begin(),offsetFactors_.end());
    }

    template <typename T> 
    inline size_t NDArray<T>::get_offset_factor_lastdim() const
    {
        if( dimensions_.size() == 0 )
            throw std::runtime_error("NDArray<T>::get_offset_factor_lastdim : array is empty");

        return get_offset_factor(dimensions_.size()-1);
    }

    template <typename T> 
    inline boost::shared_ptr< std::vector<size_t> > NDArray<T>::get_offset_factor() const
    {

        return boost::make_shared<std::vector<size_t>>(offsetFactors_.begin(),offsetFactors_.end());
    }

    template <typename T> 
    void NDArray<T>::calculate_offset_factors(const std::vector<size_t>& dimensions, std::vector<size_t>& offsetFactors)
    {
        offsetFactors.resize(dimensions.size());
        for( size_t i = 0; i < dimensions.size(); i++ )
        {
            size_t k = 1;
            for( size_t j = 0; j < i; j++ )
            {
                k *= dimensions[j];
            }

            offsetFactors[i] = k;
        }
    }

    template <typename T> 
    inline void NDArray<T>::calculate_offset_factors(const std::vector<size_t>& dimensions)
    {
        std::fill(offsetFactors_.begin(),offsetFactors_.end(),1);
        size_t a = dimensions.size();
        size_t b = offsetFactors_.size();
        size_t offsets = a<b ? a : b;
        for( size_t i = 0; i < offsets; i++ ){
            size_t k = 1;
            for( size_t j = 0; j < i; j++ )
                k *= dimensions[j];
            offsetFactors_[i] = k;
        }
    }

    template <typename T> 
    inline std::vector<size_t> NDArray<T>::calculate_index( size_t offset ) const
    {
        if( dimensions_.size() == 0 )
            throw std::runtime_error("NDArray<T>::calculate_index : array is empty");

        std::vector<size_t> index(dimensions_.size());
        for( long long i = dimensions_.size()-1; i>=0; i-- ){
            index[i] = offset / offsetFactors_[i];
            offset %= offsetFactors_[i];
        }
        return index;
    }

    template <typename T> 
    inline void NDArray<T>::calculate_index( size_t offset, std::vector<size_t>& index ) const
    {
        if( dimensions_.size() == 0 )
            throw std::runtime_error("NDArray<T>::calculate_index : array is empty");

        index.resize(dimensions_.size(), 0);
        for( long long i = dimensions_.size()-1; i>=0; i-- ){
            index[i] = offset / offsetFactors_[i];
            offset %= offsetFactors_[i];
        }
    }

    template <typename T> 
    void NDArray<T>::calculate_index( size_t offset, const std::vector<size_t>& offsetFactors, std::vector<size_t>& index )
    {
        index.resize(offsetFactors.size(), 0);

        for( long long i = offsetFactors.size()-1; i>=0; i-- )
        {
            index[i] = offset / offsetFactors[i];
            offset %= offsetFactors[i];
        }
    }

    template <typename T> 
    void NDArray<T>::clear()
    {
        if ( this->delete_data_on_destruct_ ){
            this->deallocate_memory();
        } else{
            throw std::runtime_error("NDArray<T>::clear : trying to reallocate memory not owned by array.");
        }

        this->data_ = 0;
        this->elements_ = 0;
        this->dimensions_.clear();
    }


    template <typename T> 
    inline bool NDArray<T>::point_in_range(const std::vector<size_t>& ind) const
    {
        unsigned int D = dimensions_.size();
        if ( ind.size() != D ) return false;

        unsigned int ii;
        for ( ii=0; ii<D; ii++ )
        {
            if ( ind[ii]>=dimensions_[ii] )
            {
                return false;
            }
        }

        return true;
    }

    template <typename T> 
    inline bool NDArray<T>::point_in_range(size_t x) const
    {
        GADGET_DEBUG_CHECK_THROW(dimensions_.size()==1);
        return (x<dimensions_[0]);
    }

    template <typename T> 
    inline bool NDArray<T>::point_in_range(size_t x, size_t y) const
    {
        GADGET_DEBUG_CHECK_THROW(dimensions_.size()==2);
        return ((x<dimensions_[0]) && (y<dimensions_[1]));
    }

    template <typename T> 
    inline bool NDArray<T>::point_in_range(size_t x, size_t y, size_t z) const
    {
        GADGET_DEBUG_CHECK_THROW(dimensions_.size()==3);
        return ( (x<dimensions_[0]) && (y<dimensions_[1]) && (z<dimensions_[2]));
    }

    template <typename T> 
    inline bool NDArray<T>::point_in_range(size_t x, size_t y, size_t z, size_t s) const
    {
        GADGET_DEBUG_CHECK_THROW(dimensions_.size()==4);
        return ( (x<dimensions_[0]) && (y<dimensions_[1]) && (z<dimensions_[2]) && (s<dimensions_[3]));
    }

    template <typename T> 
    inline bool NDArray<T>::point_in_range(size_t x, size_t y, size_t z, size_t s, size_t p) const
    {
        GADGET_DEBUG_CHECK_THROW(dimensions_.size()==5);
        return ( (x<dimensions_[0]) && (y<dimensions_[1]) && (z<dimensions_[2]) && (s<dimensions_[3]) && (p<dimensions_[4]));
    }

    template <typename T> 
    inline bool NDArray<T>::point_in_range(size_t x, size_t y, size_t z, size_t s, size_t p, size_t r) const
    {
        GADGET_DEBUG_CHECK_THROW(dimensions_.size()==6);
        return ( (x<dimensions_[0]) && (y<dimensions_[1]) && (z<dimensions_[2]) && (s<dimensions_[3]) && (p<dimensions_[4]) && (r<dimensions_[5]));
    }

    template <typename T> 
    inline bool NDArray<T>::point_in_range(size_t x, size_t y, size_t z, size_t s, size_t p, size_t r, size_t a) const
    {
        GADGET_DEBUG_CHECK_THROW(dimensions_.size()==7);
        return ( (x<dimensions_[0]) && (y<dimensions_[1]) && (z<dimensions_[2]) && (s<dimensions_[3]) && (p<dimensions_[4]) && (r<dimensions_[5]) && (a<dimensions_[6]));
    }

    template <typename T> 
    inline bool NDArray<T>::point_in_range(size_t x, size_t y, size_t z, size_t s, size_t p, size_t r, size_t a, size_t q) const
    {
        GADGET_DEBUG_CHECK_THROW(dimensions_.size()==8);
        return ( (x<dimensions_[0]) && (y<dimensions_[1]) && (z<dimensions_[2]) && (s<dimensions_[3]) && (p<dimensions_[4]) && (r<dimensions_[5]) && (a<dimensions_[6]) && (q<dimensions_[7]));
    }

    template <typename T> 
    inline bool NDArray<T>::point_in_range(size_t x, size_t y, size_t z, size_t s, size_t p, size_t r, size_t a, size_t q, size_t u) const
    {
        GADGET_DEBUG_CHECK_THROW(dimensions_.size()==9);
        return ( (x<dimensions_[0]) && (y<dimensions_[1]) && (z<dimensions_[2]) && (s<dimensions_[3]) && (p<dimensions_[4]) && (r<dimensions_[5]) && (a<dimensions_[6]) && (q<dimensions_[7]) && (u<dimensions_[8]));
    }




}

