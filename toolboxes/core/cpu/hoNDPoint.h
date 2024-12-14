/** \file       hoNDPoint.h
    \brief      N-dimensional point
    \author     Hui Xue
*/

#pragma once

#include <new>
#include <vector>
#include <iostream>
#include <stdexcept>
#include <cmath>
#include <cstring>
#include <limits>

#include "float.h"
#include "log.h"

namespace Gadgetron
{

    template <typename T, unsigned int D>
    class hoNDPoint
    {
    public:

        typedef hoNDPoint<T, D> Self;
        typedef T value_type;

        hoNDPoint();
        hoNDPoint(const Self& p);

        virtual ~hoNDPoint();

        Self& operator=(const Self& p);

        void fill(const T& v);

        T* begin() { return this->data_; }
        const T* begin() const { return this->data_; }

        T& operator[]( size_t idx );
        const T& operator[]( size_t idx ) const;

        T& operator()( size_t idx );
        const T& operator()( size_t idx ) const;

        bool operator==(const Self& p) const;
        bool operator!=(const Self& p) const;

        template<typename T2>
        void copyFrom(const hoNDPoint<T2, D>& aArray)
        {
            unsigned int ii;
            for ( ii=0; ii<D; ii++ )
            {
                this->data_[ii] = static_cast<T>(aArray(ii));
            }
        }

        Self& operator += (const Self& p);
        Self& operator -= (const Self& p);
        Self& operator *= (const Self& p);
        Self& operator /= (const Self& p);

        Self& operator += (const T& p);
        Self& operator -= (const T& p);
        Self& operator *= (const T& p);
        Self& operator /= (const T& p);

        // dot product
        void dot(const Self& p, T& r);

        // the magnitude of point vector
        T abs();

        // normalize the magnitude of point to be 1
        void normalize();

        virtual void print(std::ostream& os) const;

    protected:

        T data_[D];
    };

    template <typename T, unsigned int D>
    hoNDPoint<T, D>::hoNDPoint()
    {
        unsigned int ii;
        for ( ii=0; ii<D; ii++ )
        {
            this->data_[ii] = T(0);
        }
    }

    template <typename T, unsigned int D>
    hoNDPoint<T, D>::hoNDPoint(const Self& p)
    {
        memcpy(this->data_, p.data_, sizeof(T)*D);
    }

    template <typename T, unsigned int D>
    hoNDPoint<T, D>::~hoNDPoint()
    {

    }

    template <typename T, unsigned int D>
    inline hoNDPoint<T, D>& hoNDPoint<T, D>::operator=(const Self& p)
    {
        if ( this == &p ) return *this;
        memcpy(this->data_, p.data_, sizeof(T)*D);
        return *this;
    }

    template <typename T, unsigned int D>
    inline void hoNDPoint<T, D>::fill(const T& v)
    {
        unsigned int ii;
        for ( ii=0; ii<D; ii++ )
        {
            this->data_[ii] = v;
        }
    }

    template <typename T, unsigned int D>
    inline T& hoNDPoint<T, D>::operator[]( size_t idx )
    {
        GADGET_DEBUG_CHECK_THROW(idx < D);
        return this->data_[idx];
    }

    template <typename T, unsigned int D>
    inline const T& hoNDPoint<T, D>::operator[]( size_t idx ) const
    {
        GADGET_DEBUG_CHECK_THROW(idx < D);
        return this->data_[idx];
    }

    template <typename T, unsigned int D>
    inline T& hoNDPoint<T, D>::operator()( size_t idx )
    {
        GADGET_DEBUG_CHECK_THROW(idx < D);
        return this->data_[idx];
    }

    template <typename T, unsigned int D>
    inline const T& hoNDPoint<T, D>::operator()( size_t idx ) const
    {
        GADGET_DEBUG_CHECK_THROW(idx < D);
        return this->data_[idx];
    }

    template <typename T, unsigned int D>
    inline bool hoNDPoint<T, D>::operator==(const Self& p) const
    {
        T minV = std::numeric_limits<T>::epsilon();

        unsigned int ii;
        for ( ii=0; ii<D; ii++ )
        {
            if ( std::abs(this->data_[ii] - p.data_[ii]) > minV ) return false;
        }

        return true;
    }

    template <typename T, unsigned int D>
    inline bool hoNDPoint<T, D>::operator!=(const Self& p) const
    {
        return !(*this==p);
    }

    template <typename T, unsigned int D>
    inline hoNDPoint<T, D>& hoNDPoint<T, D>::operator += (const Self& p)
    {
        unsigned int ii;
        for ( ii=0; ii<D; ii++ )
        {
            this->data_[ii] += p.data_[ii];
        }

        return *this;
    }

    template <typename T, unsigned int D>
    inline hoNDPoint<T, D>& hoNDPoint<T, D>::operator -= (const Self& p)
    {
        unsigned int ii;
        for ( ii=0; ii<D; ii++ )
        {
            this->data_[ii] -= p.data_[ii];
        }

        return *this;
    }

    template <typename T, unsigned int D>
    inline hoNDPoint<T, D>& hoNDPoint<T, D>::operator *= (const Self& p)
    {
        unsigned int ii;
        for ( ii=0; ii<D; ii++ )
        {
            this->data_[ii] *= p.data_[ii];
        }

        return *this;
    }

    template <typename T, unsigned int D>
    inline hoNDPoint<T, D>& hoNDPoint<T, D>::operator /= (const Self& p)
    {
        unsigned int ii;
        for ( ii=0; ii<D; ii++ )
        {
            if ( std::abs(p.data_[ii]) < DBL_EPSILON )
            {
                this->data_[ii] /= (p.data_[ii]+DBL_EPSILON);
            }
            else
            {
                this->data_[ii] /= p.data_[ii];
            }
        }

        return *this;
    }

    template <typename T, unsigned int D>
    inline hoNDPoint<T, D>& hoNDPoint<T, D>::operator += (const T& p)
    {
        unsigned int ii;
        for ( ii=0; ii<D; ii++ )
        {
            this->data_[ii] += p;
        }

        return *this;
    }

    template <typename T, unsigned int D>
    inline hoNDPoint<T, D>& hoNDPoint<T, D>::operator -= (const T& p)
    {
        unsigned int ii;
        for ( ii=0; ii<D; ii++ )
        {
            this->data_[ii] -= p;
        }

        return *this;
    }

    template <typename T, unsigned int D>
    inline hoNDPoint<T, D>& hoNDPoint<T, D>::operator *= (const T& p)
    {
        unsigned int ii;
        for ( ii=0; ii<D; ii++ )
        {
            this->data_[ii] *= p;
        }

        return *this;
    }

    template <typename T, unsigned int D>
    inline hoNDPoint<T, D>& hoNDPoint<T, D>::operator /= (const T& p)
    {
        T pTmp = p;
        if ( std::abs(p) < DBL_EPSILON ) pTmp += DBL_EPSILON;

        unsigned int ii;
        for ( ii=0; ii<D; ii++ )
        {
            this->data_[ii] /= pTmp;
        }

        return *this;
    }

    template <typename T, unsigned int D>
    inline void hoNDPoint<T, D>::dot(const Self& p, T& r)
    {
        r = this->data_[0]*p.data_[0];

        unsigned int ii;
        for ( ii=1; ii<D; ii++ )
        {
            r += (this->data_[ii]*p.data_[ii]);
        }
    }

    template <typename T, unsigned int D>
    inline T hoNDPoint<T, D>::abs()
    {
        T dist = this->data_[0]*this->data_[0];

        unsigned int ii;
        for ( ii=1; ii<D; ii++ )
        {
            dist += (this->data_[ii]*this->data_[ii]);
        }

        dist = std::sqrt(dist);

        return dist;
    }

    template <typename T, unsigned int D>
    inline void hoNDPoint<T, D>::normalize()
    {
        T dist = this->abs();
        if ( std::abs(dist) < DBL_EPSILON ) return;

        unsigned int ii;
        for ( ii=0; ii<D; ii++ )
        {
            this->data_[ii] /= dist;
        }
    }

    template <typename T, unsigned int D>
    void hoNDPoint<T, D>::print(std::ostream& os) const
    {
        using namespace std;

        os.unsetf(std::ios::scientific);
        os.setf(ios::fixed);

        os << "[";
        unsigned int ii;
        for ( ii=0; ii<D-1; ii++ )
        {
            os << this->data_[ii] << ",";
        }
        os << this->data_[D-1] << "]";
    }
}
