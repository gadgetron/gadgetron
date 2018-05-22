/** \file   hoImageRegTransformation.h
    \brief  Define the base class for the geometric transformation in gadgetron registration
    \author Hui Xue
*/

#ifndef hoImageRegTransformation_H_
#define hoImageRegTransformation_H_

#pragma once

#include "hoNDArray.h"
#include "hoNDImage.h"
#include "hoMatrix.h"
#include "hoNDInterpolator.h"
#include "hoNDBoundaryHandler.h"
#include "hoMatrix.h"
#include "hoNDArray_utils.h"
#include "hoNDArray_elemwise.h"
#include "hoNDImage_util.h"
#include "GadgetronTimer.h"
#include "ImageIOAnalyze.h"

#ifdef USE_OMP
    #include <omp.h>
#endif // USE_OMP

namespace Gadgetron {

    enum GT_IMAGE_REG_TRANSFORMATION
    {
        GT_IMAGE_REG_TRANSFORMATION_RIGID,
        GT_IMAGE_REG_TRANSFORMATION_AFFINE,
        GT_IMAGE_REG_TRANSFORMATION_DEFORMATION_FIELD,
        GT_IMAGE_REG_TRANSFORMATION_DEFORMATION_FIELD_BIDIRECTIONAL
    };

    inline std::string getImageRegTransformationName(GT_IMAGE_REG_TRANSFORMATION v)
    {
        std::string name;

        switch (v)
        {
            case GT_IMAGE_REG_TRANSFORMATION_RIGID:
                name = "Rigid";
                break;

            case GT_IMAGE_REG_TRANSFORMATION_AFFINE:
                name = "Affine";
                break;

            case GT_IMAGE_REG_TRANSFORMATION_DEFORMATION_FIELD:
                name = "DeformationField";
                break;

            case GT_IMAGE_REG_TRANSFORMATION_DEFORMATION_FIELD_BIDIRECTIONAL:
                name = "DeformationFieldBidirectional";
                break;

            default:
                GERROR_STREAM("Unrecognized image registration transformation type : " << v);
        }

        return name;
    }

    inline GT_IMAGE_REG_TRANSFORMATION getImageRegTransformationType(const std::string& name)
    {
        GT_IMAGE_REG_TRANSFORMATION v;

        if ( name == "Rigid" )
        {
            v = GT_IMAGE_REG_TRANSFORMATION_RIGID;
        }
        else if ( name == "Affine" )
        {
            v = GT_IMAGE_REG_TRANSFORMATION_AFFINE;
        }
        else if ( name == "DeformationField" )
        {
            v = GT_IMAGE_REG_TRANSFORMATION_DEFORMATION_FIELD;
        }
        else if ( name == "DeformationFieldBidirectional" )
        {
            v = GT_IMAGE_REG_TRANSFORMATION_DEFORMATION_FIELD_BIDIRECTIONAL;
        }
        else
        {
            GERROR_STREAM("Unrecognized image registration transformation name : " << name);
        }

        return v;
    }

    /// transform a spatial position to another spatial position
    /// input and output can have different dimensions
    /// input has DIn dimension and output has DOut dimension
    /// a transformation is defined as a vector function M*1
    /// [T1; T2; T3; ...; TDOut] = T( [x1; x2; x3; ...; xDIn], [a1, a2, a3, ..., ak])
    /// transforms from n dimension to m dimension with k parameters
    /// therefore, the jacobian matrix to the parameters (Jac_parameter) is a DOut*k matrix
    /// the jacobian matrix to the spatial position (Jac_position) is a DOut*DIn matrix
    template<typename ValueType, unsigned int DIn, unsigned int DOut> 
    class hoImageRegTransformation
    {
    public:

        typedef hoImageRegTransformation<ValueType, DIn, DOut> Self;

        typedef ValueType T;
        typedef ValueType element_type;
        typedef ValueType value_type;

        typedef hoNDPoint<T, DIn> input_point_type;
        typedef hoNDPoint<T, DOut> output_point_type;

        /// there are two types of jacobian for transformations
        /// one is the jacobian to the transformation paramerters
        /// Jacobian matrix to paramters DOut*k matrix
        typedef hoMatrix<T> jacobian_parameter_type;

        /// Jacobian matrix to spatial position DOut*DIn matrix
        typedef hoMatrix<T> jacobian_position_type;

        hoImageRegTransformation() : performTiming_(false) 
        {
            gt_timer1_.set_timing_in_destruction(false);
            gt_timer2_.set_timing_in_destruction(false);
            gt_timer3_.set_timing_in_destruction(false); 
        }

        virtual ~hoImageRegTransformation() {}

        /// invert the transformation, after calling this, the transformation is replace by its inverse transformation
        virtual bool invertTransformation() = 0;

        /// set the transformation to be identical transformation
        virtual bool setIdentity() = 0;

        /// transform a point
        /// pt_in, pt_out stores a point as an array
        virtual bool transform(const T* pt_in, T* pt_out) const = 0;
        /// transform a point
        virtual bool transform( const input_point_type& in, output_point_type& out ) const;
        /// transform a group of points
        virtual bool transform( input_point_type* in, size_t N, output_point_type* out ) const;
        /// hoNDArray stores input and output points
        /// pt_in: [DIn N]; pt_out: [DOut N]
        virtual bool transform(const hoNDArray<T>& pt_in, hoNDArray<T>& pt_out) const;
        /// pt_in, pt_out stores the points as an array
        virtual bool transform(const T* pt_in, size_t N, T* pt_out) const;
        /// for the DIn==DOut
        virtual bool transform(T* pt_inout, size_t N) const;

        /// for 2D - 2D transformation
        virtual bool transform(const T& xi, const T& yi, T& xo, T& yo) const = 0;
        virtual bool transform(const T* xi, const T* yi, size_t N, T* xo, T* yo) const;
        virtual bool transform(T* x_inout, T* y_inout, size_t N) const;

        /// for 3D - 3D transformation
        virtual bool transform(const T& xi, const T& yi, const T& zi, T& xo, T& yo, T& zo) const = 0;
        virtual bool transform(const T* xi, const T* yi, const T* zi, size_t N, T* xo, T* yo, T* zo) const;
        virtual bool transform(T* x_inout, T* y_inout, T* z_inout, size_t N) const;

        /// transform a point
        /// the point is in the integer image pixel indexes
        /// image interpolator is not used
        /// pt_in, pt_out stores a point as an array
        virtual bool transform(const size_t* pt_in, T* pt_out) const = 0;
        virtual bool transform(const size_t* pt_in, size_t N, T* pt_out) const = 0;

        /// for 2D - 2D transformation
        virtual bool transform(const size_t& xi, const size_t& yi, T& xo, T& yo) const = 0;
        virtual bool transform(const size_t* xi, const size_t* yi, size_t N, T* xo, T* yo) const = 0;

        /// for 3D - 3D transformation
        virtual bool transform(const size_t& xi, const size_t& yi, const size_t& zi, T& xo, T& yo, T& zo) const = 0;
        virtual bool transform(const size_t* xi, const size_t* yi, const size_t* zi, size_t N, T* xo, T* yo, T* zo) const = 0;

        /// serialize/deserialize the transformation
        virtual bool serialize(char*& buf, size_t& len) const = 0;
        virtual bool deserialize(char* buf, size_t& len) = 0;

        virtual void print(std::ostream& os) const
        {
            using namespace std;
            os << "--------------Gagdgetron geometric transformation -------------" << endl;
            os << "Input dimension is : " << DIn << endl;
            os << "Output dimension is : " << DOut << endl;

            std::string elemTypeName = std::string(typeid(T).name());
            os << "Transformation data type is : " << elemTypeName << std::endl;
        }

        virtual std::string transformationName() const
        {
            return std::string("hoImageRegTransformation"); 
        }

        // ----------------------------------
        // debug and timing
        // ----------------------------------
        // clock for timing
        Gadgetron::GadgetronTimer gt_timer1_;
        Gadgetron::GadgetronTimer gt_timer2_;
        Gadgetron::GadgetronTimer gt_timer3_;

        bool performTiming_;

        // exporter
        Gadgetron::ImageIOAnalyze gt_exporter_;

        // debug folder
        std::string debugFolder_;
    };

    template<typename ValueType, unsigned int DIn, unsigned int DOut> 
    inline bool hoImageRegTransformation<ValueType, DIn, DOut>::
    transform( const input_point_type& in, output_point_type& out ) const
    {
        return this->transform(in.begin(), out.begin());
    }

    template<typename ValueType, unsigned int DIn, unsigned int DOut> 
    inline bool hoImageRegTransformation<ValueType, DIn, DOut>::
    transform( input_point_type* in, size_t N, output_point_type* out ) const
    {
        try
        {
            long long ii;

            #pragma omp parallel for default(none) private(ii) shared(in, out, N)
            for ( ii=0; ii<(long long)N; ii++ )
            {
                this->transform(in[ii].begin(), out[ii].begin());
            }
        }
        catch(...)
        {
            GERROR_STREAM("Errors happen in hoImageRegTransformation<ValueType, DIn, DOut>::transform( input_point_type* in, size_t N, output_point_type* out ) ... ");
            return false;
        }

        return true;
    }

    template<typename ValueType, unsigned int DIn, unsigned int DOut> 
    inline bool hoImageRegTransformation<ValueType, DIn, DOut>::
    transform(const hoNDArray<T>& pt_in, hoNDArray<T>& pt_out) const
    {
        const T* pIn = pt_in.begin();
        T* pOut = pt_out.begin();
        size_t N = pt_in.get_size(1);

        return this->transform(pIn, N, pOut);
    }

    template<typename ValueType, unsigned int DIn, unsigned int DOut> 
    inline bool hoImageRegTransformation<ValueType, DIn, DOut>::
    transform(const T* pt_in, size_t N, T* pt_out) const
    {
        try
        {
            long long ii;

            #pragma omp parallel for default(none) private(ii) shared(pt_in, N, pt_out)
            for ( ii=0; ii<(long long)N; ii++ )
            {
                this->transform(pt_in+ii*DIn, pt_out+ii*DOut);
            }
        }
        catch(...)
        {
            GERROR_STREAM("Errors happen in hoImageRegTransformation<ValueType, DIn, DOut>::transform(T* pt_in, size_t N, T* pt_out) ... ");
            return false;
        }

        return true;
    }

    template<typename ValueType, unsigned int DIn, unsigned int DOut> 
    inline bool hoImageRegTransformation<ValueType, DIn, DOut>::
    transform(T* pt_inout, size_t N) const
    {
        try
        {
            GADGET_CHECK_RETURN_FALSE(DIn>=DOut);

            long long ii;

            #pragma omp parallel default(none) private(ii) shared(pt_inout, N)
            {
                T pt_out[DOut];

                #pragma omp for 
                for ( ii=0; ii<(long long)N; ii++ )
                {
                    this->transform(pt_inout+ii*DIn, pt_out);
                    memcpy(pt_inout+ii*DIn, pt_out, sizeof(T)*DOut);
                }
            }
        }
        catch(...)
        {
            GERROR_STREAM("Errors happen in hoImageRegTransformation<ValueType, DIn, DOut>::transform(T* pt_inout, size_t N) ... ");
            return false;
        }

        return true;
    }

    template<typename ValueType, unsigned int DIn, unsigned int DOut> 
    inline bool hoImageRegTransformation<ValueType, DIn, DOut>::
    transform(const T* xi, const T* yi, size_t N, T* xo, T* yo) const
    {
        try
        {
            long long ii;

            #pragma omp parallel for default(none) private(ii) shared(xi, yi, xo, yo, N)
            for ( ii=0; ii<(long long)N; ii++ )
            {
                this->transform(xi[ii], yi[ii], xo[ii], yo[ii]);
            }
        }
        catch(...)
        {
            GERROR_STREAM("Errors happen in hoImageRegTransformation<ValueType, DIn, DOut>::transform(T* xi, T* yi, size_t N, T* xo, T* yo) ... ");
            return false;
        }

        return true;
    }

    template<typename ValueType, unsigned int DIn, unsigned int DOut> 
    inline bool hoImageRegTransformation<ValueType, DIn, DOut>::
    transform(T* x_inout, T* y_inout, size_t N) const
    {
        try
        {
            long long ii;

            T xo, yo;

            #pragma omp parallel for default(none) private(ii, xo, yo) shared(x_inout, y_inout, N)
            for ( ii=0; ii<(long long)N; ii++ )
            {
                this->transform(x_inout[ii], y_inout[ii], xo, yo);
                x_inout[ii] = xo;
                y_inout[ii] = yo;
            }
        }
        catch(...)
        {
            GERROR_STREAM("Errors happen in hoImageRegTransformation<ValueType, DIn, DOut>::transform(T* x_inout, T* y_inout, size_t N) ... ");
            return false;
        }

        return true;
    }

    template<typename ValueType, unsigned int DIn, unsigned int DOut> 
    inline bool hoImageRegTransformation<ValueType, DIn, DOut>::
    transform(const T* xi, const T* yi, const T* zi, size_t N, T* xo, T* yo, T* zo) const
    {
        try
        {
            long long ii;

            #pragma omp parallel for default(none) private(ii) shared(xi, yi, zi, xo, yo, zo, N)
            for ( ii=0; ii<(long long)N; ii++ )
            {
                this->transform(xi[ii], yi[ii], zi[ii], xo[ii], yo[ii], zo[ii]);
            }
        }
        catch(...)
        {
            GERROR_STREAM("Errors happen in hoImageRegTransformation<ValueType, DIn, DOut>::transform(T* xi, T* yi, T* zi, size_t N, T* xo, T* yo, T* zo) ... ");
            return false;
        }

        return true;
    }

    template<typename ValueType, unsigned int DIn, unsigned int DOut> 
    inline bool hoImageRegTransformation<ValueType, DIn, DOut>::
    transform(T* x_inout, T* y_inout, T* z_inout, size_t N) const
    {
        try
        {
            long long ii;

            T xo, yo, zo;

            #pragma omp parallel for default(none) private(ii, xo, yo, zo) shared(x_inout, y_inout, z_inout, N)
            for ( ii=0; ii<(long long)N; ii++ )
            {
                this->transform(x_inout[ii], y_inout[ii], z_inout[ii], xo, yo, zo);
                x_inout[ii] = xo;
                y_inout[ii] = yo;
                z_inout[ii] = zo;
            }
        }
        catch(...)
        {
            GERROR_STREAM("Errors happen in hoImageRegTransformation<ValueType, DIn, DOut>::transform(T* x_inout, T* y_inout, T* z_inout, size_t N) ... ");
            return false;
        }

        return true;
    }
}
#endif // hoImageRegTransformation_H_
