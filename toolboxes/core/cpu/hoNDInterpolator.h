/** \file       hoNDInterpolator.h
    \brief      N-dimensional interpolator

                Designed to work with hoNDArray and hoNDImage

    \author     Hui Xue
*/

#pragma once

#include "hoNDArray.h"
#include "hoNDImage.h"
#include "hoNDBoundaryHandler.h"
#include "hoNDBSpline.h"

namespace Gadgetron
{
    // define the image interpolation methods
    enum GT_IMAGE_INTERPOLATOR
    {
        GT_IMAGE_INTERPOLATOR_NEARESTNEIGHBOR=35642, // a magic number
        GT_IMAGE_INTERPOLATOR_LINEAR,
        GT_IMAGE_INTERPOLATOR_BSPLINE
    };

    inline std::string getInterpolatorName(GT_IMAGE_INTERPOLATOR interp)
    {
        std::string name;

        switch (interp)
        {
            case GT_IMAGE_INTERPOLATOR_NEARESTNEIGHBOR:
                name = "NearestNeighbor";
                break;

            case GT_IMAGE_INTERPOLATOR_LINEAR:
                name = "Linear";
                break;

            case GT_IMAGE_INTERPOLATOR_BSPLINE:
                name = "BSpline";
                break;

            default:
                GERROR_STREAM("Unrecognized interpolator type : " << interp);
        }

        return name;
    }

    inline GT_IMAGE_INTERPOLATOR getInterpolatorType(const std::string& interp_name)
    {
        GT_IMAGE_INTERPOLATOR interp;

        if ( interp_name == "NearestNeighbor" )
        {
            interp = GT_IMAGE_INTERPOLATOR_NEARESTNEIGHBOR;
        }
        else if ( interp_name == "Linear" )
        {
            interp = GT_IMAGE_INTERPOLATOR_LINEAR;
        }
        else if ( interp_name == "BSpline" )
        {
            interp = GT_IMAGE_INTERPOLATOR_BSPLINE;
        }
        else
        {
            GERROR_STREAM("Unrecognized interpolator name : " << interp_name);
        }

        return interp;
    }

    /// all interpolation calls must be made thread-safe
    template <typename ArrayType>
    class hoNDInterpolator
    {
    public:

        typedef hoNDInterpolator<ArrayType> Self;
        typedef typename ArrayType::value_type T;
        typedef hoNDBoundaryHandler<ArrayType> BoundHanlderType;
        typedef typename ArrayType::coord_type coord_type;

        hoNDInterpolator() : array_(NULL), data_(NULL), bh_(NULL), sx_(0), sy_(0), sz_(0), st_(0) {}

        hoNDInterpolator(const ArrayType& a, BoundHanlderType& bh)
        {
            array_ = &a;
            data_ = array_->begin();
            bh_ = &bh; bh_->setArray(a);

            sx_ = array_->get_size(0);
            sy_ = array_->get_size(1);
            sz_ = array_->get_size(2);
            st_ = array_->get_size(3);
        }

        virtual ~hoNDInterpolator() { array_ = NULL; bh_ = NULL; }

        virtual void setArray(const ArrayType& a)
        {
            array_ = &a;
            data_ = array_->begin();

            sx_ = array_->get_size(0);
            sy_ = array_->get_size(1);
            sz_ = array_->get_size(2);
            st_ = array_->get_size(3);
        }

        virtual void setBoundaryHandler(BoundHanlderType& bh) { bh_ = &bh; if ( array_!=NULL ) bh_->setArray(*array_); }

        /// access the pixel value
        virtual T operator()( const coord_type* pos ) = 0;
        virtual T operator()( const std::vector<coord_type>& pos ) = 0;
        virtual T operator()( coord_type x ) = 0;
        virtual T operator()( coord_type x, coord_type y ) = 0;
        virtual T operator()( coord_type x, coord_type y, coord_type z ) = 0;
        virtual T operator()( coord_type x, coord_type y, coord_type z, coord_type s ) = 0;
        virtual T operator()( coord_type x, coord_type y, coord_type z, coord_type s, coord_type p ) = 0;
        virtual T operator()( coord_type x, coord_type y, coord_type z, coord_type s, coord_type p, coord_type r ) = 0;
        virtual T operator()( coord_type x, coord_type y, coord_type z, coord_type s, coord_type p, coord_type r, coord_type a ) = 0;
        virtual T operator()( coord_type x, coord_type y, coord_type z, coord_type s, coord_type p, coord_type r, coord_type a, coord_type q ) = 0;
        virtual T operator()( coord_type x, coord_type y, coord_type z, coord_type s, coord_type p, coord_type r, coord_type a, coord_type q, coord_type u ) = 0;

    protected:

        const ArrayType* array_;
        const T* data_;
        BoundHanlderType* bh_;

        size_t sx_;
        size_t sy_;
        size_t sz_;
        size_t st_;
    };

    template <typename ArrayType>
    class hoNDInterpolatorNearestNeighbor : public hoNDInterpolator<ArrayType>
    {
    public:

        typedef hoNDInterpolator<ArrayType> BaseClass;
        typedef hoNDInterpolatorNearestNeighbor<ArrayType> Self;
        typedef typename BaseClass::T T;
        typedef typename BaseClass::coord_type coord_type;
        typedef typename BaseClass::BoundHanlderType BoundHanlderType;

        hoNDInterpolatorNearestNeighbor() : BaseClass() {}
        hoNDInterpolatorNearestNeighbor(const ArrayType& a, BoundHanlderType& bh) : BaseClass(a, bh) {}
        virtual ~hoNDInterpolatorNearestNeighbor() {}

        /// access the pixel value
        virtual T operator()( const coord_type* pos );
        virtual T operator()( const std::vector<coord_type>& pos );
        virtual T operator()( coord_type x );
        virtual T operator()( coord_type x, coord_type y );
        virtual T operator()( coord_type x, coord_type y, coord_type z );
        virtual T operator()( coord_type x, coord_type y, coord_type z, coord_type s );
        virtual T operator()( coord_type x, coord_type y, coord_type z, coord_type s, coord_type p );
        virtual T operator()( coord_type x, coord_type y, coord_type z, coord_type s, coord_type p, coord_type r );
        virtual T operator()( coord_type x, coord_type y, coord_type z, coord_type s, coord_type p, coord_type r, coord_type a );
        virtual T operator()( coord_type x, coord_type y, coord_type z, coord_type s, coord_type p, coord_type r, coord_type a, coord_type q );
        virtual T operator()( coord_type x, coord_type y, coord_type z, coord_type s, coord_type p, coord_type r, coord_type a, coord_type q, coord_type u );

    protected:

        using BaseClass::array_;
        using BaseClass::data_;
        using BaseClass::bh_;

        using BaseClass::sx_;
        using BaseClass::sy_;
        using BaseClass::sz_;
        using BaseClass::st_;
    };

    template <typename ArrayType>
    class hoNDInterpolatorLinear : public hoNDInterpolator<ArrayType>
    {
    public:

        typedef hoNDInterpolator<ArrayType> BaseClass;
        typedef hoNDInterpolatorLinear<ArrayType> Self;
        typedef typename BaseClass::T T;
        typedef typename BaseClass::coord_type coord_type;
        typedef typename BaseClass::BoundHanlderType BoundHanlderType;

        hoNDInterpolatorLinear() : BaseClass() {}

        hoNDInterpolatorLinear(const ArrayType& a, BoundHanlderType& bh) : BaseClass(a, bh)
        {
            number_of_points_ = 1<<a.get_number_of_dimensions();
        }

        virtual ~hoNDInterpolatorLinear() {}

        /// access the pixel value
        virtual T operator()( const coord_type* pos );
        virtual T operator()( const std::vector<coord_type>& pos );
        virtual T operator()( coord_type x );
        virtual T operator()( coord_type x, coord_type y );
        virtual T operator()( coord_type x, coord_type y, coord_type z );
        virtual T operator()( coord_type x, coord_type y, coord_type z, coord_type s );
        virtual T operator()( coord_type x, coord_type y, coord_type z, coord_type s, coord_type p );
        virtual T operator()( coord_type x, coord_type y, coord_type z, coord_type s, coord_type p, coord_type r );
        virtual T operator()( coord_type x, coord_type y, coord_type z, coord_type s, coord_type p, coord_type r, coord_type a );
        virtual T operator()( coord_type x, coord_type y, coord_type z, coord_type s, coord_type p, coord_type r, coord_type a, coord_type q );
        virtual T operator()( coord_type x, coord_type y, coord_type z, coord_type s, coord_type p, coord_type r, coord_type a, coord_type q, coord_type u );

    protected:

        using BaseClass::array_;
        using BaseClass::data_;
        using BaseClass::bh_;

        using BaseClass::sx_;
        using BaseClass::sy_;
        using BaseClass::sz_;
        using BaseClass::st_;

        // number of points involved in interpolation
        unsigned int number_of_points_;
    };

    template <typename ArrayType, unsigned int D>
    class hoNDInterpolatorBSpline : public hoNDInterpolator<ArrayType>
    {
    public:

        typedef hoNDInterpolator<ArrayType> BaseClass;
        typedef hoNDInterpolatorBSpline<ArrayType, D> Self;
        typedef typename BaseClass::T T;
        typedef typename BaseClass::coord_type coord_type;
        typedef typename BaseClass::BoundHanlderType BoundHanlderType;

        hoNDInterpolatorBSpline(unsigned int order=5) : BaseClass(), order_(order) { derivative_.resize(D, 0); }
        hoNDInterpolatorBSpline(const ArrayType& a, BoundHanlderType& bh, unsigned int order=5);
        virtual ~hoNDInterpolatorBSpline();

        void setArray(const ArrayType& a) override ;

        void setDerivative(const std::vector<unsigned int>& derivative) { GADGET_CHECK_THROW(derivative.size()>=D); derivative_ = derivative; }

        /// access the pixel value
        virtual T operator() ( const coord_type* pos ) override;
        virtual T operator() ( const std::vector<coord_type>& pos ) override;
        virtual T operator() ( coord_type x ) override;
        virtual T operator() ( coord_type x, coord_type y ) override;
        virtual T operator() ( coord_type x, coord_type y, coord_type z ) override;
        virtual T operator() ( coord_type x, coord_type y, coord_type z, coord_type s ) override;
        virtual T operator() ( coord_type x, coord_type y, coord_type z, coord_type s, coord_type p ) override;
        virtual T operator() ( coord_type x, coord_type y, coord_type z, coord_type s, coord_type p, coord_type r ) override;
        virtual T operator() ( coord_type x, coord_type y, coord_type z, coord_type s, coord_type p, coord_type r, coord_type a ) override;
        virtual T operator() ( coord_type x, coord_type y, coord_type z, coord_type s, coord_type p, coord_type r, coord_type a, coord_type q ) override;
        virtual T operator() ( coord_type x, coord_type y, coord_type z, coord_type s, coord_type p, coord_type r, coord_type a, coord_type q, coord_type u ) override;

     protected:

        using BaseClass::array_;
        using BaseClass::data_;
        using BaseClass::bh_;

        using BaseClass::sx_;
        using BaseClass::sy_;
        using BaseClass::sz_;
        using BaseClass::st_;

        hoNDBSpline<T, D,coord_type> bspline_;
        std::vector<size_t> dimension_;
        std::vector<unsigned int> derivative_;
        unsigned int order_;
        hoNDArray<T> coeff_;
    };

    template <typename ArrayType, unsigned int D>
    inline hoNDInterpolator<ArrayType>* createInterpolator(GT_IMAGE_INTERPOLATOR interp)
    {
        hoNDInterpolator<ArrayType>* res = NULL;

        switch (interp)
        {
            case GT_IMAGE_INTERPOLATOR_NEARESTNEIGHBOR:
                res = new hoNDInterpolatorNearestNeighbor<ArrayType>();
                break;

            case GT_IMAGE_INTERPOLATOR_LINEAR:
                res = new hoNDInterpolatorLinear<ArrayType>();
                break;

            case GT_IMAGE_INTERPOLATOR_BSPLINE:
                res = new hoNDInterpolatorBSpline<ArrayType, D>();
                break;

            default:
                GERROR_STREAM("Unrecognized interpolator type : " << interp);
        }

        return res;
    }
}

#include "hoNDInterpolatorNearestNeighbor.hxx"
#include "hoNDInterpolatorLinear.hxx"
#include "hoNDInterpolatorBSpline.hxx"
