/** \file       hoNDBoundaryHandler.h
    \brief      N-dimensional boundary condition handler

                Designed to work with hoNDArray and hoNDImage

    \author     Hui Xue
*/

#pragma once

#include "hoNDArray.h"
#include "hoNDImage.h"

namespace Gadgetron
{
    // define the boundary condition
    enum GT_BOUNDARY_CONDITION
    {
        GT_BOUNDARY_CONDITION_FIXEDVALUE=34568, // a magic number
        GT_BOUNDARY_CONDITION_BORDERVALUE,
        GT_BOUNDARY_CONDITION_PERIODIC,
        GT_BOUNDARY_CONDITION_MIRROR
    };

    inline std::string getBoundaryHandlerName(GT_BOUNDARY_CONDITION bh)
    {
        std::string name;

        switch (bh)
        {
            case GT_BOUNDARY_CONDITION_FIXEDVALUE:
                name = "FixedValue";
                break;

            case GT_BOUNDARY_CONDITION_BORDERVALUE:
                name = "BorderValue";
                break;

            case GT_BOUNDARY_CONDITION_PERIODIC:
                name = "Periodic";
                break;

            case GT_BOUNDARY_CONDITION_MIRROR:
                name = "Mirror";
                break;

            default:
                GERROR_STREAM("Unrecognized boundary handler type : " << bh);
        }

        return name;
    }

    inline GT_BOUNDARY_CONDITION getBoundaryHandlerType(const std::string& bh_name)
    {
        GT_BOUNDARY_CONDITION bh;

        if ( bh_name == "FixedValue" )
        {
            bh = GT_BOUNDARY_CONDITION_FIXEDVALUE;
        }
        else if ( bh_name == "BorderValue" )
        {
            bh = GT_BOUNDARY_CONDITION_BORDERVALUE;
        }
        else if ( bh_name == "Periodic" )
        {
            bh = GT_BOUNDARY_CONDITION_PERIODIC;
        }
        else if ( bh_name == "Mirror" )
        {
            bh = GT_BOUNDARY_CONDITION_MIRROR;
        }
        else
        {
            GERROR_STREAM("Unrecognized boundary handler name : " << bh_name);
        }

        return bh;
    }

    template <typename ArrayType>
    class hoNDBoundaryHandler
    {
    public:

        typedef hoNDBoundaryHandler<ArrayType> Self;
        typedef typename ArrayType::value_type T;

        // enum { D = ArrayType::D }

        hoNDBoundaryHandler() { array_ = NULL; }
        hoNDBoundaryHandler(const ArrayType& a) { array_ = &a; }
        virtual ~hoNDBoundaryHandler() { array_ = NULL ; }

        /// access the pixel value
        virtual T operator()( const std::vector<long long>& ind ) = 0;
        virtual T operator()( long long x ) = 0;
        virtual T operator()( long long x, long long y ) = 0;
        virtual T operator()( long long x, long long y, long long z ) = 0;
        virtual T operator()( long long x, long long y, long long z, long long s ) = 0;
        virtual T operator()( long long x, long long y, long long z, long long s, long long p ) = 0;
        virtual T operator()( long long x, long long y, long long z, long long s, long long p, long long r ) = 0;
        virtual T operator()( long long x, long long y, long long z, long long s, long long p, long long r, long long a ) = 0;
        virtual T operator()( long long x, long long y, long long z, long long s, long long p, long long r, long long a, long long q ) = 0;
        virtual T operator()( long long x, long long y, long long z, long long s, long long p, long long r, long long a, long long q, long long u ) = 0;

        void setArray(const ArrayType& a) { array_ = &a; };

        /// return a%b
        inline long long mod(long long a, long long b)
        {
            a %= b;

            if ( a<0 && b>0 )
            {
                a += b;
            }

            return a;
        }

    protected:

        const ArrayType* array_;
    };

    template <typename ArrayType>
    class hoNDBoundaryHandlerFixedValue : public hoNDBoundaryHandler<ArrayType>
    {
    public:

        typedef hoNDBoundaryHandler<ArrayType> BaseClass;
        typedef hoNDBoundaryHandlerFixedValue<ArrayType> Self;
        typedef typename BaseClass::T T;

        hoNDBoundaryHandlerFixedValue(T v=0) : BaseClass(), value_(v) {}
        hoNDBoundaryHandlerFixedValue(const ArrayType& a, T v=T(0)) : BaseClass(a), value_(v) {}
        virtual ~hoNDBoundaryHandlerFixedValue() {}

        /// access the pixel value
        virtual T operator()( const std::vector<long long>& ind );
        virtual T operator()( long long x );
        virtual T operator()( long long x, long long y );
        virtual T operator()( long long x, long long y, long long z );
        virtual T operator()( long long x, long long y, long long z, long long s );
        virtual T operator()( long long x, long long y, long long z, long long s, long long p );
        virtual T operator()( long long x, long long y, long long z, long long s, long long p, long long r );
        virtual T operator()( long long x, long long y, long long z, long long s, long long p, long long r, long long a );
        virtual T operator()( long long x, long long y, long long z, long long s, long long p, long long r, long long a, long long q );
        virtual T operator()( long long x, long long y, long long z, long long s, long long p, long long r, long long a, long long q, long long u );

    protected:
        using BaseClass::array_;
        T value_;
    };

    template <typename ArrayType>
    class hoNDBoundaryHandlerBorderValue : public hoNDBoundaryHandler<ArrayType>
    {
    public:

        typedef hoNDBoundaryHandler<ArrayType> BaseClass;
        typedef hoNDBoundaryHandlerBorderValue<ArrayType> Self;
        typedef typename BaseClass::T T;

        hoNDBoundaryHandlerBorderValue() : BaseClass() {}
        hoNDBoundaryHandlerBorderValue(const ArrayType& a) : BaseClass(a) {}
        virtual ~hoNDBoundaryHandlerBorderValue() {}

        /// access the pixel value
        virtual T operator()( const std::vector<long long>& ind );
        virtual T operator()( long long x );
        virtual T operator()( long long x, long long y );
        virtual T operator()( long long x, long long y, long long z );
        virtual T operator()( long long x, long long y, long long z, long long s );
        virtual T operator()( long long x, long long y, long long z, long long s, long long p );
        virtual T operator()( long long x, long long y, long long z, long long s, long long p, long long r );
        virtual T operator()( long long x, long long y, long long z, long long s, long long p, long long r, long long a );
        virtual T operator()( long long x, long long y, long long z, long long s, long long p, long long r, long long a, long long q );
        virtual T operator()( long long x, long long y, long long z, long long s, long long p, long long r, long long a, long long q, long long u );

    protected:
        using BaseClass::array_;
    };

    template <typename ArrayType>
    class hoNDBoundaryHandlerPeriodic : public hoNDBoundaryHandler<ArrayType>
    {
    public:

        typedef hoNDBoundaryHandler<ArrayType> BaseClass;
        typedef hoNDBoundaryHandlerPeriodic<ArrayType> Self;
        typedef typename BaseClass::T T;

        hoNDBoundaryHandlerPeriodic() : BaseClass() {}
        hoNDBoundaryHandlerPeriodic(const ArrayType& a) : BaseClass(a) {}
        virtual ~hoNDBoundaryHandlerPeriodic() {}

        /// access the pixel value
        virtual T operator()( const std::vector<long long>& ind );
        virtual T operator()( long long x );
        virtual T operator()( long long x, long long y );
        virtual T operator()( long long x, long long y, long long z );
        virtual T operator()( long long x, long long y, long long z, long long s );
        virtual T operator()( long long x, long long y, long long z, long long s, long long p );
        virtual T operator()( long long x, long long y, long long z, long long s, long long p, long long r );
        virtual T operator()( long long x, long long y, long long z, long long s, long long p, long long r, long long a );
        virtual T operator()( long long x, long long y, long long z, long long s, long long p, long long r, long long a, long long q );
        virtual T operator()( long long x, long long y, long long z, long long s, long long p, long long r, long long a, long long q, long long u );

    protected:
        using BaseClass::array_;
    };

    template <typename ArrayType>
    class hoNDBoundaryHandlerMirror : public hoNDBoundaryHandler<ArrayType>
    {
    public:

        typedef hoNDBoundaryHandler<ArrayType> BaseClass;
        typedef hoNDBoundaryHandlerMirror<ArrayType> Self;
        typedef typename BaseClass::T T;

        hoNDBoundaryHandlerMirror() : BaseClass() {}
        hoNDBoundaryHandlerMirror(const ArrayType& a) : BaseClass(a) {}
        virtual ~hoNDBoundaryHandlerMirror() {}

        /// access the pixel value
        virtual T operator()( const std::vector<long long>& ind );
        virtual T operator()( long long x );
        virtual T operator()( long long x, long long y );
        virtual T operator()( long long x, long long y, long long z );
        virtual T operator()( long long x, long long y, long long z, long long s );
        virtual T operator()( long long x, long long y, long long z, long long s, long long p );
        virtual T operator()( long long x, long long y, long long z, long long s, long long p, long long r );
        virtual T operator()( long long x, long long y, long long z, long long s, long long p, long long r, long long a );
        virtual T operator()( long long x, long long y, long long z, long long s, long long p, long long r, long long a, long long q );
        virtual T operator()( long long x, long long y, long long z, long long s, long long p, long long r, long long a, long long q, long long u );

    protected:
        using BaseClass::array_;
    };

    template <typename ArrayType> 
    hoNDBoundaryHandler<ArrayType>* createBoundaryHandler(GT_BOUNDARY_CONDITION bh)
    {
        hoNDBoundaryHandler<ArrayType>* res=NULL;

        switch (bh)
        {
            case GT_BOUNDARY_CONDITION_FIXEDVALUE:
                res = new hoNDBoundaryHandlerFixedValue<ArrayType>();
                break;

            case GT_BOUNDARY_CONDITION_BORDERVALUE:
                res = new hoNDBoundaryHandlerBorderValue<ArrayType>();
                break;

            case GT_BOUNDARY_CONDITION_PERIODIC:
                res = new hoNDBoundaryHandlerPeriodic<ArrayType>();
                break;

            case GT_BOUNDARY_CONDITION_MIRROR:
                res = new hoNDBoundaryHandlerMirror<ArrayType>();
                break;

            default:
                GERROR_STREAM("Unrecognized boundary handler type : " << bh);
        }

        return res;
    }
}

#include "hoNDBoundaryHandler.hxx"
