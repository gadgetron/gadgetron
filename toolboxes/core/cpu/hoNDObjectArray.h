/** \file   hoNDObjectArray.h
\brief  CPU-based N-dimensional array for object pointers
if delete_data_on_destruct == true, the object will be released; otherwise, only the object array memory is released
\author Hui Xue
*/

#pragma once

#include "hoNDArray.h"

namespace Gadgetron
{

    template <typename TObjectType> class hoNDObjectArray
    {
    public:

        typedef float coord_type;

        hoNDObjectArray();
        virtual ~hoNDObjectArray();

        hoNDObjectArray(const std::vector<size_t> *dimensions);
        hoNDObjectArray(const std::vector<size_t> &dimensions);
        hoNDObjectArray(boost::shared_ptr< std::vector<size_t> > dimensions);

        hoNDObjectArray(const hoNDObjectArray<TObjectType>* a);
        hoNDObjectArray(const hoNDObjectArray<TObjectType> &a);
        hoNDObjectArray<TObjectType>& operator=(const hoNDObjectArray<TObjectType>& rhs);

        void copyFrom(const hoNDObjectArray<TObjectType>& aArray);

        virtual void create(const std::vector<size_t> &dimensions);
        virtual void create(const std::vector<size_t> *dimensions);
        virtual void create(boost::shared_ptr< std::vector<size_t> > dimensions);

        size_t get_number_of_dimensions() const;
        size_t get_size(size_t dimension) const;

        bool dimensions_equal(const std::vector<size_t>& d) const;

        std::vector<size_t> get_dimensions() const;
        void get_dimensions(std::vector<size_t>& dim) const;

        size_t get_number_of_elements() const;

        size_t get_offset_factor_lastdim() const;

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
        std::vector<size_t> get_offset_factor() const;

        void calculate_offset_factors(const std::vector<size_t>& dimensions);
        static void calculate_offset_factors(const std::vector<size_t>& dimensions, std::vector<size_t>& offsetFactors);

        std::vector<size_t> calculate_index(size_t offset) const;
        void calculate_index(size_t offset, std::vector<size_t>& index) const;
        static void calculate_index(size_t offset, const std::vector<size_t>& offsetFactors, std::vector<size_t>& index);

        void clear();

        TObjectType& operator()(const std::vector<size_t>& ind);
        const TObjectType& operator()(const std::vector<size_t>& ind) const;

        TObjectType& operator()(size_t x);
        const TObjectType& operator()(size_t x) const;

        TObjectType& operator()(size_t x, size_t y);
        const TObjectType& operator()(size_t x, size_t y) const;

        TObjectType& operator()(size_t x, size_t y, size_t z);
        const TObjectType& operator()(size_t x, size_t y, size_t z) const;

        TObjectType& operator()(size_t x, size_t y, size_t z, size_t s);
        const TObjectType& operator()(size_t x, size_t y, size_t z, size_t s) const;

        TObjectType& operator()(size_t x, size_t y, size_t z, size_t s, size_t p);
        const TObjectType& operator()(size_t x, size_t y, size_t z, size_t s, size_t p) const;

        TObjectType& operator()(size_t x, size_t y, size_t z, size_t s, size_t p, size_t r);
        const TObjectType& operator()(size_t x, size_t y, size_t z, size_t s, size_t p, size_t r) const;

        TObjectType& operator()(size_t x, size_t y, size_t z, size_t s, size_t p, size_t r, size_t a);
        const TObjectType& operator()(size_t x, size_t y, size_t z, size_t s, size_t p, size_t r, size_t a) const;

        TObjectType& operator()(size_t x, size_t y, size_t z, size_t s, size_t p, size_t r, size_t a, size_t q);
        const TObjectType& operator()(size_t x, size_t y, size_t z, size_t s, size_t p, size_t r, size_t a, size_t q) const;

        TObjectType& operator()(size_t x, size_t y, size_t z, size_t s, size_t p, size_t r, size_t a, size_t q, size_t u);
        const TObjectType& operator()(size_t x, size_t y, size_t z, size_t s, size_t p, size_t r, size_t a, size_t q, size_t u) const;

        TObjectType& at(size_t idx);
        const TObjectType& at(size_t idx) const;

        TObjectType& operator[](size_t idx);
        const TObjectType& operator[](size_t idx) const;

        virtual bool serialize(char*& buf, size_t& len) const;
        virtual bool deserialize(char* buf, size_t& len);

    protected:

        std::vector<size_t> dimensions_;
        std::vector<size_t> offsetFactors_;
        size_t elements_;
        std::vector<TObjectType> data_;
    };

    template <typename TObjectType>
    hoNDObjectArray<TObjectType>::hoNDObjectArray() : elements_(0)
    {
    }

    template <typename TObjectType>
    hoNDObjectArray<TObjectType>::hoNDObjectArray(const std::vector<size_t> *dimensions)
    {
        this->create(dimensions);
    }

    template <typename TObjectType>
    hoNDObjectArray<TObjectType>::hoNDObjectArray(const std::vector<size_t> &dimensions)
    {
        this->create(dimensions);
    }

    template <typename TObjectType>
    hoNDObjectArray<TObjectType>::hoNDObjectArray(boost::shared_ptr< std::vector<size_t> > dimensions)
    {
        this->create(dimensions);
    }

    template <typename TObjectType>
    hoNDObjectArray<TObjectType>::~hoNDObjectArray()
    {
        this->clear();
    }

    template <typename TObjectType>
    hoNDObjectArray<TObjectType>::hoNDObjectArray(const hoNDObjectArray<TObjectType>* a)
    {
        this->copyFrom(*a);
    }

    template <typename TObjectType>
    hoNDObjectArray<TObjectType>::hoNDObjectArray(const hoNDObjectArray<TObjectType> &a)
    {
        this->copyFrom(a);
    }

    template <typename TObjectType>
    hoNDObjectArray<TObjectType>& hoNDObjectArray<TObjectType>::operator=(const hoNDObjectArray<TObjectType>& rhs)
    {
        if (&rhs == this) return *this;

        this->copyFrom(rhs);

        return *this;
    }

    template <typename TObjectType>
    void hoNDObjectArray<TObjectType>::create(const std::vector<size_t>& dimensions)
    {
        this->dimensions_ = dimensions;
        this->calculate_offset_factors(dimensions_);

        size_t n;

        this->elements_ = 1;
        for (n = 0; n < this->dimensions_.size(); n++)
        {
            this->elements_ *= this->dimensions_[n];
        }

        this->data_.clear();
        this->data_.resize(this->elements_);
    }

    template <typename TObjectType>
    void hoNDObjectArray<TObjectType>::create(const std::vector<size_t> *dimensions)
    {
        this->create(*dimensions);
    }

    template <typename TObjectType>
    void hoNDObjectArray<TObjectType>::create(boost::shared_ptr< std::vector<size_t> > dimensions)
    {
        this->create(*dimensions);
    }

    template <typename TObjectType>
    void hoNDObjectArray<TObjectType>::copyFrom(const hoNDObjectArray<TObjectType>& aArray)
    {
        try
        {
            if (!this->dimensions_equal(aArray.get_dimensions()))
            {
                this->create(aArray.get_dimensions());
            }

            long long i;
            for (i = 0; i<(long long)elements_; i++)
            {
                data_[i] = aArray(i);
            }
        }
        catch (...)
        {
            GADGET_THROW("Exceptions happened in hoNDObjectArray::copyFrom(...) ... ");
        }
    }

    template <typename TObjectType>
    inline size_t hoNDObjectArray<TObjectType>::get_number_of_dimensions() const
    {
        return (size_t)dimensions_.size();
    }

    template <typename TObjectType>
    inline size_t hoNDObjectArray<TObjectType>::get_size(size_t dimension) const
    {
        if (dimension >= dimensions_.size())
        {
            return 1;
        }
        else
        {
            return dimensions_[dimension];
        }
    }

    template <typename TObjectType>
    inline bool hoNDObjectArray<TObjectType>::dimensions_equal(const std::vector<size_t>& d) const
    {
        if (this->dimensions_.size() != d.size()) return false;

        size_t NDim = this->dimensions_.size();
        for (size_t ii = 0; ii<NDim; ii++)
        {
            if (this->dimensions_[ii] != d[ii]) return false;
        }

        return true;
    }

    template <typename TObjectType>
    inline std::vector<size_t> hoNDObjectArray<TObjectType>::get_dimensions() const
    {
        return this->dimensions_;
    }

    template <typename TObjectType>
    inline void hoNDObjectArray<TObjectType>::get_dimensions(std::vector<size_t>& dim) const
    {
        dim = this->dimensions_;
    }

    template <typename TObjectType>
    inline size_t hoNDObjectArray<TObjectType>::get_number_of_elements() const
    {
        return elements_;
    }

    template <typename TObjectType>
    size_t hoNDObjectArray<TObjectType>::calculate_offset(const std::vector<size_t>& ind, const std::vector<size_t>& offsetFactors)
    {
        size_t offset = ind[0];

        for (size_t i = 1; i < ind.size(); i++)
        {
            offset += ind[i] * offsetFactors[i];
        }

        return offset;
    }

    template <typename TObjectType>
    inline size_t hoNDObjectArray<TObjectType>::calculate_offset(const std::vector<size_t>& ind) const
    {
        size_t offset = ind[0];
        for (size_t i = 1; i < dimensions_.size(); i++)
            offset += ind[i] * offsetFactors_[i];
        return offset;
    }

    template <typename TObjectType>
    inline size_t hoNDObjectArray<TObjectType>::calculate_offset(size_t x, size_t y) const
    {
        return x + y * offsetFactors_[1];
    }

    template <typename TObjectType>
    inline size_t hoNDObjectArray<TObjectType>::calculate_offset(size_t x, size_t y, size_t z) const
    {
        return x + y * offsetFactors_[1] + z * offsetFactors_[2];
    }

    template <typename TObjectType>
    inline size_t hoNDObjectArray<TObjectType>::calculate_offset(size_t x, size_t y, size_t z, size_t s) const
    {
        return x + y * offsetFactors_[1] + z * offsetFactors_[2] + s * offsetFactors_[3];
    }

    template <typename TObjectType>
    inline size_t hoNDObjectArray<TObjectType>::calculate_offset(size_t x, size_t y, size_t z, size_t s, size_t p) const
    {
        return x + y * offsetFactors_[1] + z * offsetFactors_[2] + s * offsetFactors_[3] + p * offsetFactors_[4];
    }

    template <typename TObjectType>
    inline size_t hoNDObjectArray<TObjectType>::calculate_offset(size_t x, size_t y, size_t z, size_t s, size_t p, size_t r) const
    {
        return x + y * offsetFactors_[1] + z * offsetFactors_[2] + s * offsetFactors_[3] + p * offsetFactors_[4] + r * offsetFactors_[5];
    }

    template <typename TObjectType>
    inline size_t hoNDObjectArray<TObjectType>::calculate_offset(size_t x, size_t y, size_t z, size_t s, size_t p, size_t r, size_t a) const
    {
        return x + y * offsetFactors_[1] + z * offsetFactors_[2] + s * offsetFactors_[3] + p * offsetFactors_[4] + r * offsetFactors_[5] + a * offsetFactors_[6];
    }

    template <typename TObjectType>
    inline size_t hoNDObjectArray<TObjectType>::calculate_offset(size_t x, size_t y, size_t z, size_t s, size_t p, size_t r, size_t a, size_t q) const
    {
        return x + y * offsetFactors_[1] + z * offsetFactors_[2] + s * offsetFactors_[3] + p * offsetFactors_[4] + r * offsetFactors_[5] + a * offsetFactors_[6] + q * offsetFactors_[7];
    }

    template <typename TObjectType>
    inline size_t hoNDObjectArray<TObjectType>::calculate_offset(size_t x, size_t y, size_t z, size_t s, size_t p, size_t r, size_t a, size_t q, size_t u) const
    {
        return x + y * offsetFactors_[1] + z * offsetFactors_[2] + s * offsetFactors_[3] + p * offsetFactors_[4] + r * offsetFactors_[5] + a * offsetFactors_[6] + q * offsetFactors_[7] + u * offsetFactors_[8];
    }

    template <typename TObjectType>
    inline size_t hoNDObjectArray<TObjectType>::get_offset_factor(size_t dim) const
    {
        if (dim >= dimensions_.size())
            throw std::runtime_error("hoNDObjectArray<TObjectType>::get_offset_factor : index out of range");
        return offsetFactors_[dim];
    }

    template <typename TObjectType>
    inline void hoNDObjectArray<TObjectType>::get_offset_factor(std::vector<size_t>& offset) const
    {
        offset = std::vector<size_t>(offsetFactors_.begin(), offsetFactors_.end());
    }

    template <typename TObjectType>
    inline size_t hoNDObjectArray<TObjectType>::get_offset_factor_lastdim() const
    {
        if (dimensions_.size() == 0)
            throw std::runtime_error("hoNDObjectArray<TObjectType>::get_offset_factor_lastdim : array is empty");

        return get_offset_factor(dimensions_.size() - 1);
    }

    template <typename TObjectType>
    void hoNDObjectArray<TObjectType>::calculate_offset_factors(const std::vector<size_t>& dimensions, std::vector<size_t>& offsetFactors)
    {
        offsetFactors.resize(dimensions.size());
        for (size_t i = 0; i < dimensions.size(); i++)
        {
            size_t k = 1;
            for (size_t j = 0; j < i; j++)
            {
                k *= dimensions[j];
            }

            offsetFactors[i] = k;
        }
    }

    template <typename TObjectType>
    inline void hoNDObjectArray<TObjectType>::calculate_offset_factors(const std::vector<size_t>& dimensions)
    {
        this->offsetFactors_.resize(dimensions.size());
        std::fill(offsetFactors_.begin(), offsetFactors_.end(), 1);
        size_t a = dimensions.size();
        size_t b = offsetFactors_.size();
        size_t offsets = a<b ? a : b;
        for (size_t i = 0; i < offsets; i++)
        {
            size_t k = 1;
            for (size_t j = 0; j < i; j++)
                k *= dimensions[j];
            offsetFactors_[i] = k;
        }
    }

    template <typename TObjectType>
    inline std::vector<size_t> hoNDObjectArray<TObjectType>::calculate_index(size_t offset) const
    {
        if (dimensions_.size() == 0)
            throw std::runtime_error("hoNDObjectArray<TObjectType>::calculate_index : array is empty");

        std::vector<size_t> index(dimensions_.size());
        for (long long i = dimensions_.size() - 1; i >= 0; i--)
        {
            index[i] = offset / offsetFactors_[i];
            offset %= offsetFactors_[i];
        }
        return index;
    }

    template <typename TObjectType>
    inline void hoNDObjectArray<TObjectType>::calculate_index(size_t offset, std::vector<size_t>& index) const
    {
        if (dimensions_.size() == 0)
            throw std::runtime_error("hoNDObjectArray<TObjectType>::calculate_index : array is empty");

        index.resize(dimensions_.size(), 0);
        for (long long i = dimensions_.size() - 1; i >= 0; i--)
        {
            index[i] = offset / offsetFactors_[i];
            offset %= offsetFactors_[i];
        }
    }

    template <typename TObjectType>
    void hoNDObjectArray<TObjectType>::calculate_index(size_t offset, const std::vector<size_t>& offsetFactors, std::vector<size_t>& index)
    {
        index.resize(offsetFactors.size(), 0);

        for (long long i = offsetFactors.size() - 1; i >= 0; i--)
        {
            index[i] = offset / offsetFactors[i];
            offset %= offsetFactors[i];
        }
    }

    template <typename TObjectType>
    void hoNDObjectArray<TObjectType>::clear()
    {
        this->data_.clear();
        this->elements_ = 0;
        this->dimensions_.clear();
        this->offsetFactors_.clear();
    }

    template <typename TObjectType>
    inline TObjectType& hoNDObjectArray<TObjectType>::operator()(const std::vector<size_t>& ind)
    {
        size_t idx = this->calculate_offset(ind);
        GADGET_DEBUG_CHECK_THROW(idx < this->get_number_of_elements());
        return this->data_[idx];
    }

    template <typename TObjectType>
    inline const TObjectType& hoNDObjectArray<TObjectType>::operator()(const std::vector<size_t>& ind) const
    {
        size_t idx = this->calculate_offset(ind);
        GADGET_DEBUG_CHECK_THROW(idx < this->get_number_of_elements());
        return this->data_[idx];
    }

    template <typename TObjectType>
    inline TObjectType& hoNDObjectArray<TObjectType>::operator()(size_t x)
    {
        GADGET_DEBUG_CHECK_THROW(x < this->get_number_of_elements());
        return this->data_[x];
    }

    template <typename TObjectType>
    inline const TObjectType& hoNDObjectArray<TObjectType>::operator()(size_t x) const
    {
        GADGET_DEBUG_CHECK_THROW(x < this->get_number_of_elements());
        return this->data_[x];
    }

    template <typename TObjectType>
    inline TObjectType& hoNDObjectArray<TObjectType>::operator()(size_t x, size_t y)
    {
        size_t idx = this->calculate_offset(x, y);
        GADGET_DEBUG_CHECK_THROW(idx < this->get_number_of_elements());
        return this->data_[idx];
    }

    template <typename TObjectType>
    inline const TObjectType& hoNDObjectArray<TObjectType>::operator()(size_t x, size_t y) const
    {
        size_t idx = this->calculate_offset(x, y);
        GADGET_DEBUG_CHECK_THROW(idx < this->get_number_of_elements());
        return this->data_[idx];
    }

    template <typename TObjectType>
    inline TObjectType& hoNDObjectArray<TObjectType>::operator()(size_t x, size_t y, size_t z)
    {
        size_t idx = this->calculate_offset(x, y, z);
        GADGET_DEBUG_CHECK_THROW(idx < this->get_number_of_elements());
        return this->data_[idx];
    }

    template <typename TObjectType>
    inline const TObjectType& hoNDObjectArray<TObjectType>::operator()(size_t x, size_t y, size_t z) const
    {
        size_t idx = this->calculate_offset(x, y, z);
        GADGET_DEBUG_CHECK_THROW(idx < this->get_number_of_elements());
        return this->data_[idx];
    }

    template <typename TObjectType>
    inline TObjectType& hoNDObjectArray<TObjectType>::operator()(size_t x, size_t y, size_t z, size_t s)
    {
        size_t idx = this->calculate_offset(x, y, z, s);
        GADGET_DEBUG_CHECK_THROW(idx < this->get_number_of_elements());
        return this->data_[idx];
    }

    template <typename TObjectType>
    inline const TObjectType& hoNDObjectArray<TObjectType>::operator()(size_t x, size_t y, size_t z, size_t s) const
    {
        size_t idx = this->calculate_offset(x, y, z, s);
        GADGET_DEBUG_CHECK_THROW(idx < this->get_number_of_elements());
        return this->data_[idx];
    }

    template <typename TObjectType>
    inline TObjectType& hoNDObjectArray<TObjectType>::operator()(size_t x, size_t y, size_t z, size_t s, size_t p)
    {
        size_t idx = this->calculate_offset(x, y, z, s, p);
        GADGET_DEBUG_CHECK_THROW(idx < this->get_number_of_elements());
        return this->data_[idx];
    }

    template <typename TObjectType>
    inline const TObjectType& hoNDObjectArray<TObjectType>::operator()(size_t x, size_t y, size_t z, size_t s, size_t p) const
    {
        size_t idx = this->calculate_offset(x, y, z, s, p);
        GADGET_DEBUG_CHECK_THROW(idx < this->get_number_of_elements());
        return this->data_[idx];
    }

    template <typename TObjectType>
    inline TObjectType& hoNDObjectArray<TObjectType>::operator()(size_t x, size_t y, size_t z, size_t s, size_t p, size_t r)
    {
        size_t idx = this->calculate_offset(x, y, z, s, p, r);
        GADGET_DEBUG_CHECK_THROW(idx < this->get_number_of_elements());
        return this->data_[idx];
    }

    template <typename TObjectType>
    inline const TObjectType& hoNDObjectArray<TObjectType>::operator()(size_t x, size_t y, size_t z, size_t s, size_t p, size_t r) const
    {
        size_t idx = this->calculate_offset(x, y, z, s, p, r);
        GADGET_DEBUG_CHECK_THROW(idx < this->get_number_of_elements());
        return this->data_[idx];
    }

    template <typename TObjectType>
    inline TObjectType& hoNDObjectArray<TObjectType>::operator()(size_t x, size_t y, size_t z, size_t s, size_t p, size_t r, size_t a)
    {
        size_t idx = this->calculate_offset(x, y, z, s, p, r, a);
        GADGET_DEBUG_CHECK_THROW(idx < this->get_number_of_elements());
        return this->data_[idx];
    }

    template <typename TObjectType>
    inline const TObjectType& hoNDObjectArray<TObjectType>::operator()(size_t x, size_t y, size_t z, size_t s, size_t p, size_t r, size_t a) const
    {
        size_t idx = this->calculate_offset(x, y, z, s, p, r, a);
        GADGET_DEBUG_CHECK_THROW(idx < this->get_number_of_elements());
        return this->data_[idx];
    }

    template <typename TObjectType>
    inline TObjectType& hoNDObjectArray<TObjectType>::operator()(size_t x, size_t y, size_t z, size_t s, size_t p, size_t r, size_t a, size_t q)
    {
        size_t idx = this->calculate_offset(x, y, z, s, p, r, a, q);
        GADGET_DEBUG_CHECK_THROW(idx < this->get_number_of_elements());
        return this->data_[idx];
    }

    template <typename TObjectType>
    inline const TObjectType& hoNDObjectArray<TObjectType>::operator()(size_t x, size_t y, size_t z, size_t s, size_t p, size_t r, size_t a, size_t q) const
    {
        size_t idx = this->calculate_offset(x, y, z, s, p, r, a, q);
        GADGET_DEBUG_CHECK_THROW(idx < this->get_number_of_elements());
        return this->data_[idx];
    }

    template <typename TObjectType>
    inline TObjectType& hoNDObjectArray<TObjectType>::operator()(size_t x, size_t y, size_t z, size_t s, size_t p, size_t r, size_t a, size_t q, size_t u)
    {
        size_t idx = this->calculate_offset(x, y, z, s, p, r, a, q, u);
        GADGET_DEBUG_CHECK_THROW(idx < this->get_number_of_elements());
        return this->data_[idx];
    }

    template <typename TObjectType>
    inline const TObjectType& hoNDObjectArray<TObjectType>::operator()(size_t x, size_t y, size_t z, size_t s, size_t p, size_t r, size_t a, size_t q, size_t u) const
    {
        size_t idx = this->calculate_offset(x, y, z, s, p, r, a, q, u);
        GADGET_DEBUG_CHECK_THROW(idx < this->get_number_of_elements());
        return this->data_[idx];
    }

    template <typename TObjectType>
    inline TObjectType& hoNDObjectArray<TObjectType>::at(size_t idx)
    {
        GADGET_DEBUG_CHECK_THROW(idx < this->get_number_of_elements());
        return this->data_[idx];
    }

    template <typename TObjectType>
    inline const TObjectType& hoNDObjectArray<TObjectType>::at(size_t idx) const
    {
        GADGET_DEBUG_CHECK_THROW(idx < this->get_number_of_elements());
        return this->data_[idx];
    }

    template <typename TObjectType>
    inline TObjectType& hoNDObjectArray<TObjectType>::operator[](size_t idx)
    {
        GADGET_DEBUG_CHECK_THROW(idx < this->get_number_of_elements());
        return this->data_[idx];
    }

    template <typename TObjectType>
    inline const TObjectType& hoNDObjectArray<TObjectType>::operator[](size_t idx) const
    {
        GADGET_DEBUG_CHECK_THROW(idx < this->get_number_of_elements());
        return this->data_[idx];
    }

    template <typename TObjectType>
    bool hoNDObjectArray<TObjectType>::serialize(char*& buf, size_t& len) const
    {
        if (buf != NULL) delete[] buf;

        size_t NDim = dimensions_.size();

        // number of dimensions + dimension vector + contents
        len = sizeof(size_t) + sizeof(size_t) * NDim;

        // find len for every element
        size_t N = this->get_number_of_elements();

        std::vector<char*> buf_obj;
        std::vector<size_t> buf_len;
        size_t ii;

        if (N > 0)
        {
            buf_obj.resize(N, NULL);
            buf_len.resize(N, 0);

            for (ii = 0; ii < N; ii++)
            {
                this->data_[ii].serialize(buf_obj[ii], buf_len[ii]);

                len += sizeof(size_t); // buf_len[ii]
                len += buf_len[ii]; // buf[ii]
            }
        }

        buf = new char[len];

        size_t ind = 0;
        memcpy(buf, &NDim, sizeof(size_t));
        ind += sizeof(size_t);

        if (NDim > 0)
        {
            memcpy(buf + ind, &(dimensions_[0]), sizeof(size_t) * NDim);
            ind += sizeof(size_t) * NDim;

            if(N>0)
            {
                for (ii = 0; ii < N; ii++)
                {
                    memcpy(buf + ind, &(buf_len[ii]), sizeof(size_t)); ind += sizeof(size_t);
                    memcpy(buf + ind, buf_obj[ii], buf_len[ii]); ind += buf_len[ii];
                }

                for (ii = 0; ii < N; ii++)
                {
                    delete[] buf_obj[ii];
                    buf_obj[ii] = NULL;
                }
            }
        }

        return true;
    }

    template <typename TObjectType>
    bool hoNDObjectArray<TObjectType>::deserialize(char* buf, size_t& len)
    {
        size_t ind = 0;

        size_t NDim;
        memcpy(&NDim, buf, sizeof(size_t));
        ind += sizeof(size_t);

        if (NDim > 0)
        {
            std::vector<size_t> dimensions(NDim);
            memcpy(&dimensions[0], buf + ind, sizeof(size_t) * NDim);
            ind += sizeof(size_t) * NDim;

            // allocate memory
            this->create(&dimensions);

            // deserialize the content
            size_t N = this->get_number_of_elements();
            size_t ii;
            for (ii = 0; ii < N; ii++)
            {
                TObjectType a;

                size_t len_obj(0);
                memcpy(&len_obj, buf + ind, sizeof(size_t));
                ind += sizeof(size_t);

                a.deserialize(buf + ind, len_obj);
                ind += len_obj;

                this->data_[ii] = a;
            }
        }
        else
        {
            this->clear();
        }

        return true;
    }
}
