/** \file       BSplineFFD2D.h
    \brief      Implement 2D BSpline FreeFormDeformation
    \author     Hui Xue
*/

#pragma once

#include "BSplineFFD.h"

namespace Gadgetron { 

template <typename T, typename CoordType, unsigned int DOut>
class BSplineFFD2D : public BSplineFFD<T, CoordType, 2, DOut>
{
public:

    typedef BSplineFFD<T, CoordType, 2, DOut> BaseClass;
    typedef BSplineFFD2D<T, CoordType, DOut> Self;

    typedef typename BaseClass::real_value_type real_value_type;
    typedef real_value_type bspline_float_type;

    typedef typename BaseClass::coord_type coord_type;

    using BaseClass::D;
    using BaseClass::BSPLINELUTSIZE;
    using BaseClass::BSPLINEPADDINGSIZE;

    typedef typename BaseClass::LUTType             LUTType;
    typedef typename BaseClass::CoordArrayType      CoordArrayType;
    typedef typename BaseClass::ValueArrayType      ValueArrayType;
    typedef typename BaseClass::ArrayType           ArrayType;
    typedef typename BaseClass::FFDCtrlPtGridType   FFDCtrlPtGridType;
    typedef typename BaseClass::PointType           PointType;
    typedef typename BaseClass::ImageType           ImageType;

    /// constructors
    BSplineFFD2D();
    /// define the FFD over a region with specific control point spacing
    BSplineFFD2D(const PointType& start, const PointType& end, CoordType dx, CoordType dy);
    /// define the FFD over the image region with specific control point spacing
    BSplineFFD2D(const ImageType& im, CoordType dx, CoordType dy);
    /// define the FFD over the image region with specific number of control points
    BSplineFFD2D(const ImageType& im, size_t sx, size_t sy);
    /// define the FFD over an array region with specific number of control points
    BSplineFFD2D(const ArrayType& a, size_t sx, size_t sy);
    /// copy constructor
    BSplineFFD2D(const Self& bffd);

    virtual ~BSplineFFD2D();

    /// evaluate the FFD at a grid location
    virtual bool evaluateFFD(const CoordType pt[D], T r[DOut]) const;
    virtual bool evaluateFFD(CoordType px, CoordType py, T r[DOut]) const;

    virtual bool evaluateFFDDX(const CoordType pt[D], T dx[DOut]) const;
    virtual bool evaluateFFDDY(const CoordType pt[D], T dy[DOut]) const;

    virtual bool evaluateWorldDX(const CoordType pt[D], T dx[DOut]) const;
    virtual bool evaluateWorldDY(const CoordType pt[D], T dy[DOut]) const;

    /// evaluate the 1st order derivative of FFD at a grid location
    virtual bool evaluateFFDDerivative(const CoordType pt[D], T deriv[D][DOut]) const;
    virtual bool evaluateFFDDerivative(CoordType px, CoordType py, T deriv[D][DOut]) const;

    /// evaluate the 2nd order derivative of FFD at a grid location
    /// dderiv : D*D vector, stores dxx dxy dxz ...; dyx dyy dyz ...; dzx dzy dzz ...
    virtual bool evaluateFFDSecondOrderDerivative(const CoordType pt[D], T dderiv[D*D][DOut]) const;
    virtual bool evaluateFFDSecondOrderDerivative(CoordType px, CoordType py, T dderiv[D*D][DOut]) const;

    /// compute the FFD approximation once
    /// pos : the position of input points, 2 by N
    /// value : the value on input points, DOut by N
    /// residual : the approximation residual after computing FFD, DOut by N
    virtual bool ffdApprox(const CoordArrayType& pos, ValueArrayType& value, ValueArrayType& residual, real_value_type& totalResidual, size_t N);

    /// As suggested in ref [2], the BSpline FFD can be refined to achieve better approximation
    virtual bool refine();

    /// general print function
    virtual void print(std::ostream& os) const;

    using BaseClass::performTiming_;
    using BaseClass::debugFolder_;

protected:

    using BaseClass::ctrl_pt_;
    //using BaseClass::gt_timer1_;
    //using BaseClass::gt_timer2_;
    //using BaseClass::gt_timer3_;
    //using BaseClass::gt_exporter_;
    //using BaseClass::gtPlus_util_;
    //using BaseClass::gtPlus_util_complex_;
    using BaseClass::LUT_;
    using BaseClass::LUT1_;
    using BaseClass::LUT2_;

    /// evaluate the FFD
    /// px and py are at FFD grid
    /// ordx, ordy indicates the order of derivative; 0/1/2 for 0/1st/2nd derivative
    virtual bool evaluateFFD2D(CoordType px, CoordType py, size_t ordx, size_t ordy, T r[DOut]) const;
};

template <typename T, typename CoordType, unsigned int DOut>
BSplineFFD2D<T, CoordType, DOut>::BSplineFFD2D() : BaseClass()
{
}

template <typename T, typename CoordType, unsigned int DOut>
BSplineFFD2D<T, CoordType, DOut>::BSplineFFD2D(const PointType& start, const PointType& end, CoordType dx, CoordType dy) : BaseClass()
{
    GADGET_CHECK_THROW(this->initializeBFFD(start, end, dx, dy));
}

template <typename T, typename CoordType, unsigned int DOut>
BSplineFFD2D<T, CoordType, DOut>::BSplineFFD2D(const ImageType& im, CoordType dx, CoordType dy) : BaseClass()
{
    typename ImageType::coord_type x, y;

    PointType start, end;

    im.image_to_world( (size_t)0, (size_t)0, x, y);
    start(0) = x;
    start(1) = y;

    im.image_to_world(im.get_size(0)-1, im.get_size(1)-1, x, y);
    end(0) = x;
    end(1) = y;

    GADGET_CHECK_THROW(this->initializeBFFD(im, start, end, dx, dy));
}

template <typename T, typename CoordType, unsigned int DOut>
BSplineFFD2D<T, CoordType, DOut>::BSplineFFD2D(const ImageType& im, size_t sx, size_t sy) : BaseClass()
{
    PointType start, end;

    typename ImageType::coord_type x, y;

    im.image_to_world( (size_t)0, (size_t)0, x, y);
    start(0) = x;
    start(1) = y;

    im.image_to_world(im.get_size(0)-1, im.get_size(1)-1, x, y);
    end(0) = x;
    end(1) = y;

    GADGET_CHECK_THROW(this->initializeBFFD(im, start, end, sx, sy));
}

template <typename T, typename CoordType, unsigned int DOut>
BSplineFFD2D<T, CoordType, DOut>::BSplineFFD2D(const ArrayType& a, size_t sx, size_t sy) : BaseClass()
{
    PointType start, end;

    start(0) = 0;
    start(1) = 0;

    end(0) = (CoordType)(a.get_size(0)-1);
    end(1) = (CoordType)(a.get_size(1)-1);

    GADGET_CHECK_THROW(this->initializeBFFD(start, end, sx, sy));
}

template <typename T, typename CoordType, unsigned int DOut>
BSplineFFD2D<T, CoordType, DOut>::BSplineFFD2D(const Self& bffd) : BaseClass()
{
    unsigned int d;
    for ( d=0; d<DOut; d++ )
    {
        this->ctrl_pt_[d].copyFrom( bffd.get_ctrl_pt(d) );
    }
}

template <typename T, typename CoordType, unsigned int DOut>
BSplineFFD2D<T, CoordType, DOut>::~BSplineFFD2D()
{
}

template <typename T, typename CoordType, unsigned int DOut>
bool BSplineFFD2D<T, CoordType, DOut>::evaluateFFD2D(CoordType px, CoordType py, size_t ordx, size_t ordy, T r[DOut]) const
{
    try
    {
        GADGET_DEBUG_CHECK_RETURN_FALSE( (px>=-2) && (px<=this->get_size(0)+1) );
        GADGET_DEBUG_CHECK_RETURN_FALSE( (py>=-2) && (py<=this->get_size(1)+1) );

        GADGET_DEBUG_CHECK_RETURN_FALSE(ordx>=0 && ordx<=2);
        GADGET_DEBUG_CHECK_RETURN_FALSE(ordy>=0 && ordy<=2);
        GADGET_DEBUG_CHECK_RETURN_FALSE(ordx+ordy<=2);

        long long ix = (long long)std::floor(px);
        CoordType deltaX = px-(CoordType)ix;
        long long lx = FFD_MKINT(BSPLINELUTSIZE*deltaX);

        long long iy = (long long)std::floor(py);
        CoordType deltaY = py-(CoordType)iy;
        long long ly = FFD_MKINT(BSPLINELUTSIZE*deltaY);

        unsigned int d, jj;
        size_t offset[4];
        offset[0] = this->calculate_offset(ix-1, iy-1);
        offset[1] = this->calculate_offset(ix-1, iy);
        offset[2] = this->calculate_offset(ix-1, iy+1);
        offset[3] = this->calculate_offset(ix-1, iy+2);

        const LUTType* p_xLUT= &this->LUT_;
        const LUTType* p_yLUT= &this->LUT_;

        if ( ordx == 1 )
        {
            p_xLUT= &this->LUT1_;
        }
        else if ( ordx == 2 )
        {
            p_xLUT= &this->LUT2_;
        }

        if ( ordy == 1 )
        {
            p_yLUT= &this->LUT1_;
        }
        else if ( ordy == 2 )
        {
            p_yLUT= &this->LUT2_;
        }

        const LUTType& xLUT= *p_xLUT;
        const LUTType& yLUT= *p_yLUT;

        for ( d=0; d<DOut; d++ )
        {
            r[d] = 0;

            T v(0);
            for (jj=0; jj<4; jj++)
            {
                v =   ( this->ctrl_pt_[d](offset[jj]  ) * xLUT[lx][0] )
                    + ( this->ctrl_pt_[d](offset[jj]+1) * xLUT[lx][1] )
                    + ( this->ctrl_pt_[d](offset[jj]+2) * xLUT[lx][2] )
                    + ( this->ctrl_pt_[d](offset[jj]+3) * xLUT[lx][3] );

                r[d] += v * yLUT[ly][jj];
            }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Error happened in evaluateFFD2D(CoordType px, CoordType py, size_t ordx, size_t ordy, T r[DOut]) const ... ");
        return false;
    }

    return true;
}

template <typename T, typename CoordType, unsigned int DOut>
inline bool BSplineFFD2D<T, CoordType, DOut>::evaluateFFD(const CoordType pt[D], T r[DOut]) const
{
    return this->evaluateFFD2D(pt[0], pt[1], 0, 0, r);
}

template <typename T, typename CoordType, unsigned int DOut>
inline bool BSplineFFD2D<T, CoordType, DOut>::evaluateFFD(CoordType px, CoordType py, T r[DOut]) const
{
    return this->evaluateFFD2D(px, py, 0, 0, r);
}

template <typename T, typename CoordType, unsigned int DOut>
inline bool BSplineFFD2D<T, CoordType, DOut>::evaluateFFDDX(const CoordType pt[D], T dx[DOut]) const
{
    return this->evaluateFFD2D(pt[0], pt[1], 1, 0, dx);
}

template <typename T, typename CoordType, unsigned int DOut>
inline bool BSplineFFD2D<T, CoordType, DOut>::evaluateFFDDY(const CoordType pt[D], T dy[DOut]) const
{
    return this->evaluateFFD2D(pt[0], pt[1], 0, 1, dy);
}

template <typename T, typename CoordType, unsigned int DOut>
inline bool BSplineFFD2D<T, CoordType, DOut>::evaluateWorldDX(const CoordType pt[D], T dx[DOut]) const
{
    GADGET_CHECK_RETURN_FALSE(this->evaluateFFD2D(pt[0], pt[1], 1, 0, dx));
    coord_type sx = coord_type(1.0)/this->get_spacing(0);
    unsigned int d;
    for ( d=0; d<DOut; d++ )
    {
        dx[d] *= sx;
    }
    return true;
}

template <typename T, typename CoordType, unsigned int DOut>
inline bool BSplineFFD2D<T, CoordType, DOut>::evaluateWorldDY(const CoordType pt[D], T dy[DOut]) const
{
    GADGET_CHECK_RETURN_FALSE(this->evaluateFFD2D(pt[0], pt[1], 0, 1, dy));
    coord_type sy = coord_type(1.0)/this->get_spacing(1);
    unsigned int d;
    for ( d=0; d<DOut; d++ )
    {
        dy[d] *= sy;
    }
    return true;
}

template <typename T, typename CoordType, unsigned int DOut>
inline bool BSplineFFD2D<T, CoordType, DOut>::evaluateFFDDerivative(const CoordType pt[D], T deriv[D][DOut]) const
{
    GADGET_CHECK_RETURN_FALSE(this->evaluateFFD2D(pt[0], pt[1], 1, 0, deriv[0]));
    GADGET_CHECK_RETURN_FALSE(this->evaluateFFD2D(pt[0], pt[1], 0, 1, deriv[1]));
    return true;
}

template <typename T, typename CoordType, unsigned int DOut>
inline bool BSplineFFD2D<T, CoordType, DOut>::evaluateFFDDerivative(CoordType px, CoordType py, T deriv[D][DOut]) const
{
    GADGET_CHECK_RETURN_FALSE(this->evaluateFFD2D(px, py, 1, 0, deriv[0]));
    GADGET_CHECK_RETURN_FALSE(this->evaluateFFD2D(px, py, 0, 1, deriv[1]));
    return true;
}

template <typename T, typename CoordType, unsigned int DOut>
inline bool BSplineFFD2D<T, CoordType, DOut>::evaluateFFDSecondOrderDerivative(const CoordType pt[D], T dderiv[D*D][DOut]) const
{
    GADGET_CHECK_RETURN_FALSE(this->evaluateFFD2D(pt[0], pt[1], 2, 0, dderiv[0]));
    GADGET_CHECK_RETURN_FALSE(this->evaluateFFD2D(pt[0], pt[1], 1, 1, dderiv[1]));
    memcpy(dderiv[2], dderiv[1], DOut*sizeof(T));
    GADGET_CHECK_RETURN_FALSE(this->evaluateFFD2D(pt[0], pt[1], 0, 2, dderiv[3]));
    return true;
}

template <typename T, typename CoordType, unsigned int DOut>
inline bool BSplineFFD2D<T, CoordType, DOut>::evaluateFFDSecondOrderDerivative(CoordType px, CoordType py, T dderiv[D*D][DOut]) const
{
    GADGET_CHECK_RETURN_FALSE(this->evaluateFFD2D(px, py, 2, 0, dderiv[0]));
    GADGET_CHECK_RETURN_FALSE(this->evaluateFFD2D(px, py, 1, 1, dderiv[1]));
    memcpy(dderiv[2], dderiv[1], DOut*sizeof(T));
    GADGET_CHECK_RETURN_FALSE(this->evaluateFFD2D(px, py, 0, 2, dderiv[3]));
    return true;
}

template <typename T, typename CoordType, unsigned int DOut>
bool BSplineFFD2D<T, CoordType, DOut>::ffdApprox(const CoordArrayType& pos, ValueArrayType& value, ValueArrayType& residual, real_value_type& totalResidual, size_t N)
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(pos.get_size(0)==2);
        GADGET_CHECK_RETURN_FALSE(pos.get_size(1)==N);

        GADGET_CHECK_RETURN_FALSE(value.get_size(0)==DOut);
        GADGET_CHECK_RETURN_FALSE(value.get_size(1)==N);

        if ( !residual.dimensions_equal(&value) )
        {
            residual.create(value.dimensions());
            Gadgetron::clear(residual);
        }

        size_t sx = this->get_size(0);
        size_t sy = this->get_size(1);

        /// following the definition of ref[2]
        ho3DArray<T> dx(sx, sy, DOut), ds(sx, sy, DOut);
        Gadgetron::clear(dx);
        Gadgetron::clear(ds);

        /// compute the current approximation values
        ValueArrayType approxValue;
        approxValue = value;

        /// compute current residual
        GADGET_CHECK_RETURN_FALSE(this->evaluateFFDArray(pos, approxValue));
        GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::subtract(value, approxValue, residual));

        /// compute the update of control points
        unsigned int d;

        long long n;
        for (n=0; n<(long long)N; n++)
        {
            coord_type px = pos(0, n);
            coord_type py = pos(1, n);

            if ( px<-2 || px>sx+2 || py<-2 || py>sy+2 )
            {
                continue;
            }

            long long ix = (long long)std::floor(px);
            CoordType deltaX = px-(CoordType)ix;

            long long iy = (long long)std::floor(py);
            CoordType deltaY = py-(CoordType)iy;

            long long i, j, I, J;

            T dist=0, v, vv, vvv;
            for (j=0; j<4; j++)
            {
                for (i=0; i<4; i++)
                {
                    v = this->BSpline(i, deltaX) * this->BSpline(j, deltaY);
                    dist += v*v;
                }
            }

            for (j=0; j<4; j++)
            {
                J = j + iy - 1;
                if ( (J>=0) && (J<(long long)sy) )
                {
                    for (i=0; i<4; i++)
                    {
                        I = i + ix - 1;
                        if ( (I>=0) && (I<(long long)sx) )
                        {
                            v = this->BSpline(i, deltaX) * this->BSpline(j, deltaY);
                            vv = v*v;
                            vvv = vv*v;

                            for ( d=0; d<DOut; d++ )
                            {
                                dx(I, J, d) += vvv*residual(d, n)/dist;
                                ds(I, J, d) += vv;
                            }
                        }
                    }
                }
            }
        }

        /// update the control point values
        GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::addEpsilon(ds));
        GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::divide(dx, ds, dx));

        std::vector<size_t> startND(2, BSPLINEPADDINGSIZE), size(2);
        size[0] = sx;
        size[1] = sy;

        hoNDArray<T> ctrlPtWithoutPadding(sx, sy);

        for ( d=0; d<DOut; d++ )
        {
            hoNDArray<T> dx2D(sx, sy, dx.begin()+d*sx*sy*sizeof(T));

            auto dim = this->ctrl_pt_[d].dimensions();
            hoNDArray<T> tmpCtrlPt(dim, this->ctrl_pt_[d].begin(), false);

            vector_td<size_t, 2> crop_offset;
            crop_offset[0] = startND[0];
            crop_offset[1] = startND[1];

            vector_td<size_t, 2> crop_size;
            crop_size[0] = size[0];
            crop_size[1] = size[1];

            Gadgetron::crop(crop_offset, crop_size, tmpCtrlPt, ctrlPtWithoutPadding);
            Gadgetron::add(ctrlPtWithoutPadding, dx2D, ctrlPtWithoutPadding);
            Gadgetron::fill(ctrlPtWithoutPadding, crop_offset, tmpCtrlPt);
        }

        /*for (j=0; j<sy; j++)
        {
            for (i=0; i<sx; i++)
            {
                for ( d=0; d<DOut; d++ )
                {
                    if ( ds(i, j, d) > 0)
                    {
                        this->ctrl_pt_[d](i, j) += dx(i, j, d)/ds(i, j, d);
                    }
                }
            }
        }*/

        /// calculate residual error
        totalResidual = 0;
        GADGET_CHECK_RETURN_FALSE(this->evaluateFFDArray(pos, approxValue));
        GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::subtract(value, approxValue, residual));
        totalResidual = Gadgetron::nrm2(residual);
        totalResidual = totalResidual / (real_value_type)N;
    }
    catch(...)
    {
        GERROR_STREAM("Error happened in ffdApprox(const CoordArrayType& pos, ValueArrayType& value, ValueArrayType& residual, real_value_type& totalResidual, size_t N) ... ");
        return false;
    }

    return true;
}

template <typename T, typename CoordType, unsigned int DOut>
bool BSplineFFD2D<T, CoordType, DOut>::refine()
{
    try
    {
        size_t sx = this->get_size(0);
        size_t sy = this->get_size(1);

        /// the refined control point grid definition

        std::vector<size_t> dim(2);
        dim[0] = 2*sx-1 + 2*BSPLINEPADDINGSIZE;
        dim[1] = 2*sy-1 + 2*BSPLINEPADDINGSIZE;

        std::vector<coord_type> spacing;
        this->get_spacing(spacing);
        spacing[0] /= 2;
        spacing[1] /= 2;

        std::vector<coord_type> oldOrigin;
        this->ctrl_pt_[0].get_origin(oldOrigin);

        std::vector<coord_type> gridOrigin(2);
        this->ctrl_pt_[0].image_to_world( (CoordType)(BSPLINEPADDINGSIZE), (CoordType)(BSPLINEPADDINGSIZE), gridOrigin[0], gridOrigin[1]);

        std::vector<coord_type> origin(2);
        origin[0] = (oldOrigin[0] + gridOrigin[0])/2;
        origin[1] = (oldOrigin[1] + gridOrigin[1])/2;

        typename ImageType::axis_type axis;
        this->ctrl_pt_[0].get_axis(axis);

        /// allocate new control points
        FFDCtrlPtGridType new_ctrl_pt[DOut];

        unsigned int d;
        for( d=0; d<DOut; d++ )
        {
            new_ctrl_pt[d].create(dim, spacing, origin, axis);
            Gadgetron::clear(new_ctrl_pt[d]);
        }

        /// refinement weights, see ref[2]
        T w[2][3];

        w[0][0] = T(0.125); w[0][1] = T(0.75);  w[0][2] = T(0.125);
        w[1][0] = 0;        w[1][1] = T(0.5);   w[1][2] = T(0.5);

        /// compute refined control point values
        int x, y, i_new, j_new, i_old, j_old;
        for (y=0; y<sy; y++)
        {
            for (x=0; x<sx; x++)
            {
                for (j_new=0; j_new<2; j_new++)
                {
                    for (i_new=0; i_new<2; i_new++)
                    {
                        size_t offsetNew = new_ctrl_pt[0].calculate_offset(2*x+i_new+BSPLINEPADDINGSIZE, 2*y+j_new+BSPLINEPADDINGSIZE);

                        for (j_old=0; j_old<3; j_old++)
                        {
                            for (i_old=0; i_old<3; i_old++)
                            {
                                size_t offsetOld = this->calculate_offset(x+i_old-1, y+j_old-1);
                                for ( d=0; d<DOut; d++ )
                                {
                                    new_ctrl_pt[d](offsetNew) += w[i_new][i_old]*w[j_new][j_old] * this->ctrl_pt_[d](offsetOld);
                                }
                            }
                        }
                    }
                }
            }
        }

        for ( d=0; d<DOut; d++ )
        {
            this->ctrl_pt_[d].create(dim, spacing, origin, axis, new_ctrl_pt[d].begin(), true);
            new_ctrl_pt[d].delete_data_on_destruct(false);
        }
    }
    catch(...)
    {
        GERROR_STREAM("Error happened in refine() ... ");
        return false;
    }

    return true;
}

template <typename T, typename CoordType, unsigned int DOut>
void BSplineFFD2D<T, CoordType, DOut>::print(std::ostream& os) const
{
    using namespace std;

    os << "---------------------- BSpline 2D Free Form Deformation ------------------" << endl;
    os << "Implement 2D BSpline Free Form Deformation (BFFD) " << endl;

    std::string elemTypeName = std::string(typeid(T).name());
    os << "FFD value type is : " << elemTypeName << endl;

    elemTypeName = std::string(typeid(CoordType).name());
    os << "FFD coord type is : " << elemTypeName << endl;

    os << "Output dimension is : " << DOut << endl;
    os << "---------------------------------------------------" << endl;
    os << "BFFD grid information : " << endl;
    this->ctrl_pt_[0].printContent(os);
    os << "------------------------------------------------------------------------------" << endl;
}

}
