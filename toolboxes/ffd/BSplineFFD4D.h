/** \file       BSplineFFD4D.h
    \brief      Implement 4D BSpline FreeFormDeformation
    \author     Hui Xue
*/

#pragma once

#include "FFDBase.h"

namespace Gadgetron { 

template <typename T, typename CoordType, unsigned int DOut>
class BSplineFFD4D : public BSplineFFD<T, CoordType, 4, DOut>
{
public:

    typedef BSplineFFD<T, CoordType, 4, DOut> BaseClass;
    typedef BSplineFFD4D<T, CoordType, DOut> Self;

    typedef typename BaseClass::real_value_type real_value_type;
    typedef real_value_type bspline_float_type;

    typedef typename BaseClass::coord_type coord_type;

    enum { D = 4 };
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
    BSplineFFD4D();
    /// define the FFD over a region with specific control point spacing
    BSplineFFD4D(const PointType& start, const PointType& end, CoordType dx, CoordType dy, CoordType dz, CoordType ds);
    /// define the FFD over the image region with specific control point spacing
    BSplineFFD4D(const ImageType& im, CoordType dx, CoordType dy, CoordType dz, CoordType ds);
    /// define the FFD over the image region with specific number of control points
    BSplineFFD4D(const ImageType& im, size_t sx, size_t sy, size_t sz, size_t ss);
    /// define the FFD over an array region with specific number of control points
    BSplineFFD4D(const ArrayType& a, size_t sx, size_t sy, size_t sz, size_t ss);
    /// copy constructor
    BSplineFFD4D(const Self& bffd);

    virtual ~BSplineFFD4D();

    /// evaluate the FFD at a grid location
    virtual bool evaluateFFD(const CoordType pt[D], T r[DOut]) const;
    virtual bool evaluateFFD(CoordType px, CoordType py, CoordType pz, CoordType ps, T r[DOut]) const;

    virtual bool evaluateFFDDX(const CoordType pt[D], T dx[DOut]) const;
    virtual bool evaluateFFDDY(const CoordType pt[D], T dy[DOut]) const;
    virtual bool evaluateFFDDZ(const CoordType pt[D], T dz[DOut]) const;
    virtual bool evaluateFFDDS(const CoordType pt[D], T ds[DOut]) const;

    virtual bool evaluateWorldDX(const CoordType pt[D], T dx[DOut]) const;
    virtual bool evaluateWorldDY(const CoordType pt[D], T dy[DOut]) const;
    virtual bool evaluateWorldDZ(const CoordType pt[D], T dz[DOut]) const;
    virtual bool evaluateWorldDS(const CoordType pt[D], T ds[DOut]) const;

    /// evaluate the 1st order derivative of FFD at a grid location
    virtual bool evaluateFFDDerivative(const CoordType pt[D], T deriv[D][DOut]) const;
    virtual bool evaluateFFDDerivative(CoordType px, CoordType py, CoordType pz, CoordType ps, T deriv[D][DOut]) const;

    /// evaluate the 2nd order derivative of FFD at a grid location
    /// dderiv : D*D vector, stores dxx dxy dxz ...; dyx dyy dyz ...; dzx dzy dzz ...
    virtual bool evaluateFFDSecondOrderDerivative(const CoordType pt[D], T dderiv[D*D][DOut]) const;
    virtual bool evaluateFFDSecondOrderDerivative(CoordType px, CoordType py, CoordType pz, CoordType ps, T dderiv[D*D][DOut]) const;

    /// compute the FFD approximation once
    /// pos : the position of input points, D by N
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
    virtual bool evaluateFFD4D(CoordType px, CoordType py, CoordType pz, CoordType ps, size_t ordx, size_t ordy, size_t ordz, size_t ords, T r[DOut]) const;
};

template <typename T, typename CoordType, unsigned int DOut>
BSplineFFD4D<T, CoordType, DOut>::BSplineFFD4D() : BaseClass()
{
}

template <typename T, typename CoordType, unsigned int DOut>
BSplineFFD4D<T, CoordType, DOut>::BSplineFFD4D(const PointType& start, const PointType& end, CoordType dx, CoordType dy, CoordType dz, CoordType ds) : BaseClass()
{
    GADGET_CHECK_THROW(this->initializeBFFD(start, end, dx, dy, dz, ds));
}

template <typename T, typename CoordType, unsigned int DOut>
BSplineFFD4D<T, CoordType, DOut>::BSplineFFD4D(const ImageType& im, CoordType dx, CoordType dy, CoordType dz, CoordType ds) : BaseClass()
{
    PointType start, end;

    typename ImageType::coord_type x, y, z, s;

    im.image_to_world(0, 0, 0, 0, x, y, z, s);
    start(0) = x;
    start(1) = y;
    start(2) = z;
    start(3) = s;

    im.image_to_world(im.get_size(0)-1, im.get_size(1)-1, im.get_size(2)-1, im.get_size(3)-1, x, y, z, s);
    end(0) = x;
    end(1) = y;
    end(2) = z;
    end(3) = s;

    GADGET_CHECK_THROW(this->initializeBFFD(im, start, end, dx, dy, dz, ds));
}

template <typename T, typename CoordType, unsigned int DOut>
BSplineFFD4D<T, CoordType, DOut>::BSplineFFD4D(const ImageType& im, size_t sx, size_t sy, size_t sz, size_t ss) : BaseClass()
{
    PointType start, end;

    typename ImageType::coord_type x, y, z, s;

    im.image_to_world( (size_t)0, (size_t)0, (size_t)0, (size_t)0, x, y, z, s);
    start(0) = x;
    start(1) = y;
    start(2) = z;
    start(3) = s;

    im.image_to_world(im.get_size(0)-1, im.get_size(1)-1, im.get_size(2)-1, im.get_size(3)-1, x, y, z, s);
    end(0) = x;
    end(1) = y;
    end(2) = z;
    end(3) = s;

    GADGET_CHECK_THROW(this->initializeBFFD(im, start, end, sx, sy, sz, ss));
}

template <typename T, typename CoordType, unsigned int DOut>
BSplineFFD4D<T, CoordType, DOut>::BSplineFFD4D(const ArrayType& a, size_t sx, size_t sy, size_t sz, size_t ss) : BaseClass()
{
    PointType start, end;

    start(0) = 0;
    start(1) = 0;
    start(2) = 0;
    start(3) = 0;

    end(0) = a.get_size(0)-1;
    end(1) = a.get_size(1)-1;
    end(2) = a.get_size(2)-1;
    end(3) = a.get_size(3)-1;

    GADGET_CHECK_THROW(this->initializeBFFD(start, end, sx, sy, sz, ss));
}

template <typename T, typename CoordType, unsigned int DOut>
BSplineFFD4D<T, CoordType, DOut>::BSplineFFD4D(const Self& bffd) : BaseClass()
{
    unsigned int d;
    for ( d=0; d<DOut; d++ )
    {
        this->ctrl_pt_[d].copyFrom( bffd.get_ctrl_pt(d) );
    }
}

template <typename T, typename CoordType, unsigned int DOut>
BSplineFFD4D<T, CoordType, DOut>::~BSplineFFD4D()
{
}

template <typename T, typename CoordType, unsigned int DOut>
bool BSplineFFD4D<T, CoordType, DOut>::evaluateFFD4D(CoordType px, CoordType py, CoordType pz, CoordType ps, size_t ordx, size_t ordy, size_t ordz, size_t ords, T r[DOut]) const
{
    try
    {
        GADGET_DEBUG_CHECK_RETURN_FALSE( (px>=-2) && (px<=this->get_size(0)+1) );
        GADGET_DEBUG_CHECK_RETURN_FALSE( (py>=-2) && (py<=this->get_size(1)+1) );
        GADGET_DEBUG_CHECK_RETURN_FALSE( (pz>=-2) && (pz<=this->get_size(2)+1) );
        GADGET_DEBUG_CHECK_RETURN_FALSE( (ps>=-2) && (ps<=this->get_size(3)+1) );

        GADGET_DEBUG_CHECK_RETURN_FALSE(ordx>=0 && ordx<=2);
        GADGET_DEBUG_CHECK_RETURN_FALSE(ordy>=0 && ordy<=2);
        GADGET_DEBUG_CHECK_RETURN_FALSE(ordz>=0 && ordz<=2);
        GADGET_DEBUG_CHECK_RETURN_FALSE(ords>=0 && ords<=2);
        GADGET_DEBUG_CHECK_RETURN_FALSE(ordx+ordy+ordz+ords<=2);

        long long ix = (long long)std::floor(px);
        CoordType deltaX = px-(CoordType)ix;
        long long lx = FFD_MKINT(BSPLINELUTSIZE*deltaX);

        long long iy = (long long)std::floor(py);
        CoordType deltaY = py-(CoordType)iy;
        long long ly = FFD_MKINT(BSPLINELUTSIZE*deltaY);

        long long iz = (long long)std::floor(pz);
        CoordType deltaZ = pz-(CoordType)iz;
        long long lz = FFD_MKINT(BSPLINELUTSIZE*deltaZ);

        long long is = (long long)std::floor(ps);
        CoordType deltaS = ps-(CoordType)is;
        long long ls = FFD_MKINT(BSPLINELUTSIZE*deltaS);

        unsigned int d, jj, kk, ss;
        size_t offset[4][4][4]; // s, z, y

        for ( ss=0; ss<4; ss++ )
        {
            offset[ss][0][0] = this->calculate_offset(ix-1, iy-1, iz-1, is+ss-1);
            offset[ss][0][1] = this->calculate_offset(ix-1, iy  , iz-1, is+ss-1);
            offset[ss][0][2] = this->calculate_offset(ix-1, iy+1, iz-1, is+ss-1);
            offset[ss][0][3] = this->calculate_offset(ix-1, iy+2, iz-1, is+ss-1);

            offset[ss][1][0] = this->calculate_offset(ix-1, iy-1, iz, is+ss-1);
            offset[ss][1][1] = this->calculate_offset(ix-1, iy  , iz, is+ss-1);
            offset[ss][1][2] = this->calculate_offset(ix-1, iy+1, iz, is+ss-1);
            offset[ss][1][3] = this->calculate_offset(ix-1, iy+2, iz, is+ss-1);

            offset[ss][2][0] = this->calculate_offset(ix-1, iy-1, iz+1, is+ss-1);
            offset[ss][2][1] = this->calculate_offset(ix-1, iy  , iz+1, is+ss-1);
            offset[ss][2][2] = this->calculate_offset(ix-1, iy+1, iz+1, is+ss-1);
            offset[ss][2][3] = this->calculate_offset(ix-1, iy+2, iz+1, is+ss-1);

            offset[ss][3][0] = this->calculate_offset(ix-1, iy-1, iz+2, is+ss-1);
            offset[ss][3][1] = this->calculate_offset(ix-1, iy  , iz+2, is+ss-1);
            offset[ss][3][2] = this->calculate_offset(ix-1, iy+1, iz+2, is+ss-1);
            offset[ss][3][3] = this->calculate_offset(ix-1, iy+2, iz+2, is+ss-1);
        }

        const LUTType* p_xLUT= &this->LUT_;
        const LUTType* p_yLUT= &this->LUT_;
        const LUTType* p_zLUT= &this->LUT_;
        const LUTType* p_sLUT= &this->LUT_;

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

        if ( ordz == 1 )
        {
            p_zLUT= &this->LUT1_;
        }
        else if ( ordz == 2 )
        {
            p_zLUT= &this->LUT2_;
        }

        if ( ords == 1 )
        {
            p_sLUT= &this->LUT1_;
        }
        else if ( ords == 2 )
        {
            p_sLUT= &this->LUT2_;
        }

        const LUTType& xLUT= *p_xLUT;
        const LUTType& yLUT= *p_yLUT;
        const LUTType& zLUT= *p_zLUT;
        const LUTType& sLUT= *p_sLUT;

        for ( d=0; d<DOut; d++ )
        {
            r[d] = 0;
            for (ss=0; ss<4; ss++)
            {
                T rs=0;
                for (kk=0; kk<4; kk++)
                {
                    T rv = 0;
                    for (jj=0; jj<4; jj++)
                    {
                        T v  =  ( this->ctrl_pt_[d](offset[ss][kk][jj])   * xLUT[lx][0] )
                            + ( this->ctrl_pt_[d](offset[ss][kk][jj]+1) * xLUT[lx][1] )
                            + ( this->ctrl_pt_[d](offset[ss][kk][jj]+2) * xLUT[lx][2] )
                            + ( this->ctrl_pt_[d](offset[ss][kk][jj]+3) * xLUT[lx][3] );

                        rv += v * yLUT[ly][jj];
                    }

                    rs += rv * zLUT[lz][kk];
                }

                r[d] += rs * sLUT[ls][ss];
            }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Error happened in evaluateFFD4D(CoordType px, CoordType py, CoordType pz, CoordType ps, size_t ordx, size_t ordy, size_t ordz, size_t ords, T r[DOut]) const ... ");
        return false;
    }

    return true;
}

template <typename T, typename CoordType, unsigned int DOut>
inline bool BSplineFFD4D<T, CoordType, DOut>::evaluateFFD(const CoordType pt[D], T r[DOut]) const
{
    return this->evaluateFFD4D(pt[0], pt[1], pt[2], pt[3], 0, 0, 0, 0, r);
}

template <typename T, typename CoordType, unsigned int DOut>
inline bool BSplineFFD4D<T, CoordType, DOut>::evaluateFFD(CoordType px, CoordType py, CoordType pz, CoordType ps, T r[DOut]) const
{
    return this->evaluateFFD4D(px, py, pz, ps, 0, 0, 0, 0, r);
}

template <typename T, typename CoordType, unsigned int DOut>
inline bool BSplineFFD4D<T, CoordType, DOut>::evaluateFFDDX(const CoordType pt[D], T dx[DOut]) const
{
    return this->evaluateFFD4D(pt[0], pt[1], pt[2], pt[3], 1, 0, 0, 0, dx);
}

template <typename T, typename CoordType, unsigned int DOut>
inline bool BSplineFFD4D<T, CoordType, DOut>::evaluateFFDDY(const CoordType pt[D], T dy[DOut]) const
{
    return this->evaluateFFD4D(pt[0], pt[1], pt[2], pt[3], 0, 1, 0, 0, dy);
}

template <typename T, typename CoordType, unsigned int DOut>
inline bool BSplineFFD4D<T, CoordType, DOut>::evaluateFFDDZ(const CoordType pt[D], T dz[DOut]) const
{
    return this->evaluateFFD4D(pt[0], pt[1], pt[2], pt[3], 0, 0, 1, 0, dz);
}

template <typename T, typename CoordType, unsigned int DOut>
inline bool BSplineFFD4D<T, CoordType, DOut>::evaluateFFDDS(const CoordType pt[D], T ds[DOut]) const
{
    return this->evaluateFFD4D(pt[0], pt[1], pt[2], pt[3], 0, 0, 0, 1, ds);
}

template <typename T, typename CoordType, unsigned int DOut>
inline bool BSplineFFD4D<T, CoordType, DOut>::evaluateWorldDX(const CoordType pt[D], T dx[DOut]) const
{
    GADGET_CHECK_RETURN_FALSE(this->evaluateFFD4D(pt[0], pt[1], pt[2], pt[3], 1, 0, 0, 0, dx));
    coord_type sx = coord_type(1.0)/this->get_spacing(0);
    unsigned int d;
    for ( d=0; d<DOut; d++ )
    {
        dx[d] *= sx;
    }
    return true;
}

template <typename T, typename CoordType, unsigned int DOut>
inline bool BSplineFFD4D<T, CoordType, DOut>::evaluateWorldDY(const CoordType pt[D], T dy[DOut]) const
{
    GADGET_CHECK_RETURN_FALSE(this->evaluateFFD4D(pt[0], pt[1], pt[2], pt[3], 0, 1, 0, 0, dy));
    coord_type sy = coord_type(1.0)/this->get_spacing(1);
    unsigned int d;
    for ( d=0; d<DOut; d++ )
    {
        dy[d] *= sy;
    }
    return true;
}

template <typename T, typename CoordType, unsigned int DOut>
inline bool BSplineFFD4D<T, CoordType, DOut>::evaluateWorldDZ(const CoordType pt[D], T dz[DOut]) const
{
    GADGET_CHECK_RETURN_FALSE(this->evaluateFFD4D(pt[0], pt[1], pt[2], pt[3], 0, 0, 1, 0, dz));
    coord_type sz = coord_type(1.0)/this->get_spacing(2);
    unsigned int d;
    for ( d=0; d<DOut; d++ )
    {
        dz[d] *= sz;
    }
    return true;
}

template <typename T, typename CoordType, unsigned int DOut>
inline bool BSplineFFD4D<T, CoordType, DOut>::evaluateWorldDS(const CoordType pt[D], T ds[DOut]) const
{
    GADGET_CHECK_RETURN_FALSE(this->evaluateFFD4D(pt[0], pt[1], pt[2], pt[3], 0, 0, 0, 1, ds));
    coord_type ss = coord_type(1.0)/this->get_spacing(3);
    unsigned int d;
    for ( d=0; d<DOut; d++ )
    {
        ds[d] *= ss;
    }
    return true;
}

template <typename T, typename CoordType, unsigned int DOut>
inline bool BSplineFFD4D<T, CoordType, DOut>::evaluateFFDDerivative(const CoordType pt[D], T deriv[D][DOut]) const
{
    GADGET_CHECK_RETURN_FALSE(this->evaluateFFD4D(pt[0], pt[1], pt[2], pt[3], 1, 0, 0, 0, deriv[0]));
    GADGET_CHECK_RETURN_FALSE(this->evaluateFFD4D(pt[0], pt[1], pt[2], pt[3], 0, 1, 0, 0, deriv[1]));
    GADGET_CHECK_RETURN_FALSE(this->evaluateFFD4D(pt[0], pt[1], pt[2], pt[3], 0, 0, 1, 0, deriv[2]));
    GADGET_CHECK_RETURN_FALSE(this->evaluateFFD4D(pt[0], pt[1], pt[2], pt[3], 0, 0, 0, 1, deriv[3]));
    return true;
}

template <typename T, typename CoordType, unsigned int DOut>
inline bool BSplineFFD4D<T, CoordType, DOut>::evaluateFFDDerivative(CoordType px, CoordType py, CoordType pz, CoordType ps, T deriv[D][DOut]) const
{
    GADGET_CHECK_RETURN_FALSE(this->evaluateFFD4D(px, py, pz, ps, 1, 0, 0, 0, deriv[0]));
    GADGET_CHECK_RETURN_FALSE(this->evaluateFFD4D(px, py, pz, ps, 0, 1, 0, 0, deriv[1]));
    GADGET_CHECK_RETURN_FALSE(this->evaluateFFD4D(px, py, pz, ps, 0, 0, 1, 0, deriv[2]));
    GADGET_CHECK_RETURN_FALSE(this->evaluateFFD4D(px, py, pz, ps, 0, 0, 0, 1, deriv[3]));
    return true;
}

template <typename T, typename CoordType, unsigned int DOut>
inline bool BSplineFFD4D<T, CoordType, DOut>::evaluateFFDSecondOrderDerivative(const CoordType pt[D], T dderiv[D*D][DOut]) const
{
    // dxx
    GADGET_CHECK_RETURN_FALSE(this->evaluateFFD4D(pt[0], pt[1], pt[2], pt[3], 2, 0, 0, 0, dderiv[0]));
    // dxy
    GADGET_CHECK_RETURN_FALSE(this->evaluateFFD4D(pt[0], pt[1], pt[2], pt[3], 1, 1, 0, 0, dderiv[1]));
    // dxz
    GADGET_CHECK_RETURN_FALSE(this->evaluateFFD4D(pt[0], pt[1], pt[2], pt[3], 1, 0, 1, 0, dderiv[2]));
    // dxs
    GADGET_CHECK_RETURN_FALSE(this->evaluateFFD4D(pt[0], pt[1], pt[2], pt[3], 1, 0, 0, 1, dderiv[3]));

    // dyx
    memcpy(dderiv[4], dderiv[1], DOut*sizeof(T));
    // dyy
    GADGET_CHECK_RETURN_FALSE(this->evaluateFFD4D(pt[0], pt[1], pt[2], pt[3], 0, 2, 0, 0, dderiv[5]));
    // dyz
    GADGET_CHECK_RETURN_FALSE(this->evaluateFFD4D(pt[0], pt[1], pt[2], pt[3], 0, 1, 1, 0, dderiv[6]));
    // dys
    GADGET_CHECK_RETURN_FALSE(this->evaluateFFD4D(pt[0], pt[1], pt[2], pt[3], 0, 1, 0, 1, dderiv[7]));

    // dzx
    memcpy(dderiv[8], dderiv[2], DOut*sizeof(T));
    // dzy
    memcpy(dderiv[9], dderiv[6], DOut*sizeof(T));
    // dzz
    GADGET_CHECK_RETURN_FALSE(this->evaluateFFD4D(pt[0], pt[1], pt[2], pt[3], 0, 0, 2, 0, dderiv[10]));
    // dzs
    GADGET_CHECK_RETURN_FALSE(this->evaluateFFD4D(pt[0], pt[1], pt[2], pt[3], 0, 0, 1, 1, dderiv[11]));

    // dsx
    memcpy(dderiv[12], dderiv[3], DOut*sizeof(T));
    // dsy
    memcpy(dderiv[13], dderiv[7], DOut*sizeof(T));
    // dsz
    memcpy(dderiv[14], dderiv[11], DOut*sizeof(T));
    // dss
    GADGET_CHECK_RETURN_FALSE(this->evaluateFFD4D(pt[0], pt[1], pt[2], pt[3], 0, 0, 0, 2, dderiv[15]));

    return true;
}

template <typename T, typename CoordType, unsigned int DOut>
inline bool BSplineFFD4D<T, CoordType, DOut>::evaluateFFDSecondOrderDerivative(CoordType px, CoordType py, CoordType pz, CoordType ps, T dderiv[D*D][DOut]) const
{
    // dxx
    GADGET_CHECK_RETURN_FALSE(this->evaluateFFD4D(px, py, pz, ps, 2, 0, 0, 0, dderiv[0]));
    // dxy
    GADGET_CHECK_RETURN_FALSE(this->evaluateFFD4D(px, py, pz, ps, 1, 1, 0, 0, dderiv[1]));
    // dxz
    GADGET_CHECK_RETURN_FALSE(this->evaluateFFD4D(px, py, pz, ps, 1, 0, 1, 0, dderiv[2]));
    // dxs
    GADGET_CHECK_RETURN_FALSE(this->evaluateFFD4D(px, py, pz, ps, 1, 0, 0, 1, dderiv[3]));

    // dyx
    memcpy(dderiv[4], dderiv[1], DOut*sizeof(T));
    // dyy
    GADGET_CHECK_RETURN_FALSE(this->evaluateFFD4D(px, py, pz, ps, 0, 2, 0, 0, dderiv[5]));
    // dyz
    GADGET_CHECK_RETURN_FALSE(this->evaluateFFD4D(px, py, pz, ps, 0, 1, 1, 0, dderiv[6]));
    // dys
    GADGET_CHECK_RETURN_FALSE(this->evaluateFFD4D(px, py, pz, ps, 0, 1, 0, 1, dderiv[7]));

    // dzx
    memcpy(dderiv[8], dderiv[2], DOut*sizeof(T));
    // dzy
    memcpy(dderiv[9], dderiv[6], DOut*sizeof(T));
    // dzz
    GADGET_CHECK_RETURN_FALSE(this->evaluateFFD4D(px, py, pz, ps, 0, 0, 2, 0, dderiv[10]));
    // dzs
    GADGET_CHECK_RETURN_FALSE(this->evaluateFFD4D(px, py, pz, ps, 0, 0, 1, 1, dderiv[11]));

    // dsx
    memcpy(dderiv[12], dderiv[3], DOut*sizeof(T));
    // dsy
    memcpy(dderiv[13], dderiv[7], DOut*sizeof(T));
    // dsz
    memcpy(dderiv[14], dderiv[11], DOut*sizeof(T));
    // dss
    GADGET_CHECK_RETURN_FALSE(this->evaluateFFD4D(px, py, pz, ps, 0, 0, 0, 2, dderiv[15]));

    return true;
}

template <typename T, typename CoordType, unsigned int DOut>
bool BSplineFFD4D<T, CoordType, DOut>::ffdApprox(const CoordArrayType& pos, ValueArrayType& value, ValueArrayType& residual, real_value_type& totalResidual, size_t N)
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(pos.get_size(0)==D);
        GADGET_CHECK_RETURN_FALSE(pos.get_size(1)==N);

        GADGET_CHECK_RETURN_FALSE(value.get_size(0)==DOut);
        GADGET_CHECK_RETURN_FALSE(value.get_size(1)==N);

        std::vector<size_t> dim;
        value.get_dimensions(dim);
        if ( !residual.dimensions_equal(&dim) )
        {
            residual.create(value.get_dimensions());
            Gadgetron::clear(residual);
        }

        size_t sx = this->get_size(0);
        size_t sy = this->get_size(1);
        size_t sz = this->get_size(2);
        size_t ss = this->get_size(3);

        /// following the definition of ref[2]
        ho5DArray<T> dx(sx, sy, sz, ss, DOut), ds(sx, sy, sz, ss, DOut);
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
        for (n=0; n<N; n++)
        {
            coord_type px = pos(0, n);
            coord_type py = pos(1, n);
            coord_type pz = pos(2, n);
            coord_type ps = pos(3, n);

            if ( px<-2 || px>sx+2
                || py<-2 || py>sy+2
                || pz<-2 || pz>sz+2
                || ps<-2 || ps>ss+2 )
            {
                continue;
            }

            long long ix = (long long)std::floor(px);
            CoordType deltaX = px-(CoordType)ix;

            long long iy = (long long)std::floor(py);
            CoordType deltaY = py-(CoordType)iy;

            long long iz = (long long)std::floor(pz);
            CoordType deltaZ = pz-(CoordType)iz;

            long long is = (long long)std::floor(ps);
            CoordType deltaS = ps-(CoordType)is;

            long long i, j, k, s, I, J, K, S;

            T dist=0, v, vv, vvv;
            for (s=0; s<4; s++)
            {
                for (k=0; k<4; k++)
                {
                    for (j=0; j<4; j++)
                    {
                        for (i=0; i<4; i++)
                        {
                            v = (this->BSpline(i, deltaX) * this->BSpline(j, deltaY)) * (this->BSpline(k, deltaZ) * this->BSpline(s, deltaS));
                            dist += v*v;
                        }
                    }
                }
            }

            for (s=0; s<4; s++)
            {
                S = s + is - 1;
                if ( (S>=0) && (S<(long long)ss) )
                {
                    for (k=0; k<4; k++)
                    {
                        K = k + iz - 1;
                        if ( (K>=0) && (K<(long long)sz) )
                        {
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
                                            v = this->BSpline(i, deltaX) * this->BSpline(j, deltaY) * this->BSpline(k, deltaZ) * this->BSpline(s, deltaS);
                                            vv = v*v;
                                            vvv = vv*v;

                                            for ( d=0; d<DOut; d++ )
                                            {
                                                dx(I, J, K, S, d) += vvv*residual(d, n)/dist;
                                                ds(I, J, K, S, d) += vv;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        /// update the control point values
        GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::addEpsilon(ds));
        GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::divide(dx, ds, dx));

        std::vector<size_t> startND(4, BSPLINEPADDINGSIZE), size(4);
        size[0] = sx;
        size[1] = sy;
        size[2] = sz;
        size[3] = ss;

        hoNDArray<T> ctrlPtWithoutPadding(sx, sy, sz, ss);

        for ( d=0; d<DOut; d++ )
        {
            hoNDArray<T> dx4D(sx, sy, sz, ss, dx.begin()+d*sx*sy*sz*ss*sizeof(T));

            std::vector<size_t> dim;
            this->ctrl_pt_[d].get_dimensions(dim);
            hoNDArray<T> tmpCtrlPt(dim, this->ctrl_pt_[d].begin(), false);

            vector_td<size_t, 4> crop_offset;
            crop_offset[0] = startND[0];
            crop_offset[1] = startND[1];
            crop_offset[2] = startND[2];
            crop_offset[3] = startND[3];

            vector_td<size_t, 4> crop_size;
            crop_size[0] = size[0];
            crop_size[1] = size[1];
            crop_size[2] = size[2];
            crop_size[3] = size[3];

            Gadgetron::crop(crop_offset, crop_size, tmpCtrlPt, ctrlPtWithoutPadding);
            Gadgetron::add(ctrlPtWithoutPadding, dx4D, ctrlPtWithoutPadding);
            Gadgetron::fill(ctrlPtWithoutPadding, crop_offset, tmpCtrlPt);
        }

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
bool BSplineFFD4D<T, CoordType, DOut>::refine()
{
    try
    {
        size_t sx = this->get_size(0);
        size_t sy = this->get_size(1);
        size_t sz = this->get_size(2);
        size_t ss = this->get_size(3);

        /// the refined control point grid definition

        std::vector<size_t> dim(4);
        dim[0] = 2*sx-1 + 2*BSPLINEPADDINGSIZE;
        dim[1] = 2*sy-1 + 2*BSPLINEPADDINGSIZE;
        dim[2] = 2*sz-1 + 2*BSPLINEPADDINGSIZE;
        dim[3] = 2*ss-1 + 2*BSPLINEPADDINGSIZE;

        std::vector<coord_type> spacing;
        this->get_spacing(spacing);
        spacing[0] /= 2;
        spacing[1] /= 2;
        if ( sz > 1 ) spacing[2] /= 2;
        if ( ss > 1 ) spacing[3] /= 2;

        std::vector<coord_type> oldOrigin;
        this->ctrl_pt_[0].get_origin(oldOrigin);

        std::vector<coord_type> gridOrigin(4);
        this->ctrl_pt_[0].image_to_world( (CoordType)(BSPLINEPADDINGSIZE),
                                          (CoordType)(BSPLINEPADDINGSIZE),
                                          (CoordType)(BSPLINEPADDINGSIZE),
                                          (CoordType)(BSPLINEPADDINGSIZE),
                                          gridOrigin[0], gridOrigin[1],
                                          gridOrigin[2], gridOrigin[3]);

        std::vector<coord_type> origin(4);
        origin[0] = (oldOrigin[0] + gridOrigin[0])/2;
        origin[1] = (oldOrigin[1] + gridOrigin[1])/2;
        origin[2] = (oldOrigin[2] + gridOrigin[2])/2;
        origin[3] = (oldOrigin[3] + gridOrigin[3])/2;

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
        int x, y, z, s, i_new, j_new, k_new, s_new, i_old, j_old, k_old, s_old;

        if ( ss>1 && sz>1 )
        {
            for (s=0; s<ss; s++)
            {
                for (z=0; z<sz; z++)
                {
                    for (y=0; y<sy; y++)
                    {
                        for (x=0; x<sx; x++)
                        {
                            for (s_new=0; s_new<2; s_new++)
                            {
                                for (k_new=0; k_new<2; k_new++)
                                {
                                    for (j_new=0; j_new<2; j_new++)
                                    {
                                        for (i_new=0; i_new<2; i_new++)
                                        {
                                            size_t offsetNew = new_ctrl_pt[0].calculate_offset(2*x+i_new+BSPLINEPADDINGSIZE, 2*y+j_new+BSPLINEPADDINGSIZE, 2*z+k_new+BSPLINEPADDINGSIZE, 2*s+s_new+BSPLINEPADDINGSIZE);

                                            for (s_old=0; s_old<3; s_old++)
                                            {
                                                for (k_old=0; k_old<3; k_old++)
                                                {
                                                    for (j_old=0; j_old<3; j_old++)
                                                    {
                                                        for (i_old=0; i_old<3; i_old++)
                                                        {
                                                            size_t offsetOld = this->calculate_offset(x+i_old-1, y+j_old-1, z+k_old-1, s+s_old-1);
                                                            for ( d=0; d<DOut; d++ )
                                                            {
                                                                new_ctrl_pt[d](offsetNew) += w[i_new][i_old]*w[j_new][j_old]*w[k_new][k_old]*w[s_new][s_old] * this->ctrl_pt_[d](offsetOld);
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        else if ( ss==1 && sz>1 )
        {
            for (z=0; z<sz; z++)
            {
                for (y=0; y<sy; y++)
                {
                    for (x=0; x<sx; x++)
                    {
                        for (k_new=0; k_new<2; k_new++)
                        {
                            for (j_new=0; j_new<2; j_new++)
                            {
                                for (i_new=0; i_new<2; i_new++)
                                {
                                    size_t offsetNew = new_ctrl_pt[0].calculate_offset(2*x+i_new+BSPLINEPADDINGSIZE, 2*y+j_new+BSPLINEPADDINGSIZE, 2*z+k_new+BSPLINEPADDINGSIZE, BSPLINEPADDINGSIZE);

                                    for (k_old=0; k_old<3; k_old++)
                                    {
                                        for (j_old=0; j_old<3; j_old++)
                                        {
                                            for (i_old=0; i_old<3; i_old++)
                                            {
                                                size_t offsetOld = this->calculate_offset(x+i_old-1, y+j_old-1, z+k_old-1, 0);
                                                for ( d=0; d<DOut; d++ )
                                                {
                                                    new_ctrl_pt[d](offsetNew) += w[i_new][i_old]*w[j_new][j_old]*w[k_new][k_old] * this->ctrl_pt_[d](offsetOld);
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        else if ( ss==1 && sz==1 )
        {
            for (y=0; y<sy; y++)
            {
                for (x=0; x<sx; x++)
                {
                    for (j_new=0; j_new<2; j_new++)
                    {
                        for (i_new=0; i_new<2; i_new++)
                        {
                            size_t offsetNew = new_ctrl_pt[0].calculate_offset(2*x+i_new+BSPLINEPADDINGSIZE, 2*y+j_new+BSPLINEPADDINGSIZE, BSPLINEPADDINGSIZE, BSPLINEPADDINGSIZE);

                            for (j_old=0; j_old<3; j_old++)
                            {
                                for (i_old=0; i_old<3; i_old++)
                                {
                                    size_t offsetOld = this->calculate_offset(x+i_old-1, y+j_old-1, 0, 0);
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
void BSplineFFD4D<T, CoordType, DOut>::print(std::ostream& os) const
{
    using namespace std;

    os << "---------------------- BSpline 4D Free Form Deformation ------------------" << endl;
    os << "Implement 4D BSpline Free Form Deformation (BFFD) " << endl;

    std::string elemTypeName = std::string( typeid(T).name() );
    os << "FFD value type is : " << elemTypeName << endl;

    elemTypeName = std::string( typeid(CoordType).name() );
    os << "FFD coord type is : " << elemTypeName << endl;

    os << "Output dimension is : " << DOut << endl;
    os << "---------------------------------------------------" << endl;
    os << "BFFD grid information : " << endl;
    this->ctrl_pt_[0].printContent(os);
    os << "---------------------------------------------------------------------------------" << endl;
}

}
