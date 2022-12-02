/** \file       BSplineFFD.h
    \brief      Class for BSpline FreeFormDeformation
    \author     Hui Xue
*/

#pragma once

#include "FFDBase.h"

namespace Gadgetron { 

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut>
class BSplineFFD : public FFDBase<T, CoordType, DIn, DOut>
{
public:

    typedef Gadgetron::FFDBase<T, CoordType, DIn, DOut> BaseClass;
    typedef BSplineFFD<T, CoordType, DIn, DOut> Self;

    typedef typename BaseClass::real_value_type real_value_type;
    typedef real_value_type bspline_float_type;

    typedef typename BaseClass::coord_type coord_type;

    using BaseClass::D;
    enum { BSPLINELUTSIZE = 1000 };
    enum { BSPLINEPADDINGSIZE = 4 };

    typedef real_value_type LUTType[BSPLINELUTSIZE][BSPLINEPADDINGSIZE];

    typedef typename BaseClass::CoordArrayType      CoordArrayType;
    typedef typename BaseClass::ValueArrayType      ValueArrayType;
    typedef typename BaseClass::ArrayType           ArrayType;
    typedef typename BaseClass::FFDCtrlPtGridType   FFDCtrlPtGridType;
    typedef typename BaseClass::PointType           PointType;
    typedef typename BaseClass::ImageType           ImageType;
    typedef typename BaseClass::MaskArrayType       MaskArrayType;

    BSplineFFD();
    virtual ~BSplineFFD();
    /// although BSpline grid has the padding, every index is defined on the unpadded grid

    /// get the size of control point arrays
    virtual size_t get_size(size_t dimension) const override { return ctrl_pt_[0].get_size(dimension)-2*BSPLINEPADDINGSIZE; }
    virtual std::vector<size_t> get_dimensions() const
    {
        std::vector<size_t> dim;
        ctrl_pt_[0].get_dimensions(dim);

        unsigned int d;
        for ( d=0; d<DIn; d++ )
        {
            dim[d] -= 2*BSPLINEPADDINGSIZE;
        }

        return dim;
    }

    /// get the spacing of of control point arrays
    coord_type get_spacing(size_t dimension) const override { return ctrl_pt_[0].get_pixel_size(dimension); }
    void get_spacing(std::vector<coord_type>& spacing) const override { ctrl_pt_[0].get_pixel_size(spacing); }

    /// get/set a control point value
    T get(size_t x, size_t y, size_t d) const override { return ctrl_pt_[d](x+BSPLINEPADDINGSIZE, y+BSPLINEPADDINGSIZE); }
    void set(size_t x, size_t y, size_t d, T v) override { ctrl_pt_[d](x+BSPLINEPADDINGSIZE, y+BSPLINEPADDINGSIZE) = v; }

    T get(size_t x, size_t y, size_t z, size_t d) const override { return ctrl_pt_[d](x+BSPLINEPADDINGSIZE, y+BSPLINEPADDINGSIZE, z+BSPLINEPADDINGSIZE); }
    void set(size_t x, size_t y, size_t z, size_t d, T v) override { ctrl_pt_[d](x+BSPLINEPADDINGSIZE, y+BSPLINEPADDINGSIZE, z+BSPLINEPADDINGSIZE) = v; }

    T get(size_t x, size_t y, size_t z, size_t s, size_t d) const override { return ctrl_pt_[d](x+BSPLINEPADDINGSIZE, y+BSPLINEPADDINGSIZE, z+BSPLINEPADDINGSIZE, s+BSPLINEPADDINGSIZE); }
    void set(size_t x, size_t y, size_t z, size_t s, size_t d, T v) override { ctrl_pt_[d](x+BSPLINEPADDINGSIZE, y+BSPLINEPADDINGSIZE, z+BSPLINEPADDINGSIZE, s+BSPLINEPADDINGSIZE) = v; }

    /// offset to/from indexes for control points
    size_t calculate_offset(size_t x, size_t y) const override  { return ctrl_pt_[0].calculate_offset(x+BSPLINEPADDINGSIZE, y+BSPLINEPADDINGSIZE); }

    void calculate_index( size_t offset, size_t& x, size_t& y ) const override
    {
        ctrl_pt_[0].calculate_index(offset, x, y);
        x -= BSPLINEPADDINGSIZE;
        y -= BSPLINEPADDINGSIZE;
    }

    size_t calculate_offset(size_t x, size_t y, size_t z) const override { return ctrl_pt_[0].calculate_offset(x+BSPLINEPADDINGSIZE, y+BSPLINEPADDINGSIZE, z+BSPLINEPADDINGSIZE); }
    void calculate_index( size_t offset, size_t& x, size_t& y, size_t& z ) const override
    {
        ctrl_pt_[0].calculate_index(offset, x, y, z);
        x -= BSPLINEPADDINGSIZE;
        y -= BSPLINEPADDINGSIZE;
        z -= BSPLINEPADDINGSIZE;
    }

    size_t calculate_offset(size_t x, size_t y, size_t z, size_t s) const override { return ctrl_pt_[0].calculate_offset(x+BSPLINEPADDINGSIZE, y+BSPLINEPADDINGSIZE, z+BSPLINEPADDINGSIZE, s+BSPLINEPADDINGSIZE); }
    void calculate_index( size_t offset, size_t& x, size_t& y, size_t& z, size_t& s ) const override
    {
        ctrl_pt_[0].calculate_index(offset, x, y, z, s);
        x -= BSPLINEPADDINGSIZE;
        y -= BSPLINEPADDINGSIZE;
        z -= BSPLINEPADDINGSIZE;
        s -= BSPLINEPADDINGSIZE;
    }

    /// compute the control point location in world coordinates
    void get_location(size_t x, size_t y, CoordType& sx, CoordType& sy) const override { ctrl_pt_[0].image_to_world(x+BSPLINEPADDINGSIZE, y+BSPLINEPADDINGSIZE, sx, sy); }
    void get_location(size_t x, size_t y, size_t z, CoordType& sx, CoordType& sy, CoordType& sz) const override { ctrl_pt_[0].image_to_world(x+BSPLINEPADDINGSIZE, y+BSPLINEPADDINGSIZE, z+BSPLINEPADDINGSIZE, sx, sy, sz); }
    void get_location(size_t x, size_t y, size_t z, size_t s, CoordType& sx, CoordType& sy, CoordType& sz, CoordType& ss) const override { ctrl_pt_[0].image_to_world(x+BSPLINEPADDINGSIZE, y+BSPLINEPADDINGSIZE, z+BSPLINEPADDINGSIZE, s+BSPLINEPADDINGSIZE, sx, sy, sz, ss); }

    /// convert a world coordinate point to FFD grid location
    bool world_to_grid(const CoordType pt_w[D], CoordType pt_g[D]) const override;
    bool world_to_grid(CoordType px_w, CoordType py_w, CoordType& px_g, CoordType& py_g) const override;
    bool world_to_grid(CoordType px_w, CoordType py_w, CoordType pz_w, CoordType& px_g, CoordType& py_g, CoordType& pz_g) const override;
    bool world_to_grid(CoordType px_w, CoordType py_w, CoordType pz_w, CoordType ps_w, CoordType& px_g, CoordType& py_g, CoordType& pz_g, CoordType& ps_g) const override;

    bool grid_to_world(const CoordType pt_g[D], CoordType pt_w[D]) const override;
    bool grid_to_world(CoordType px_g, CoordType py_g, CoordType& px_w, CoordType& py_w) const override;
    bool grid_to_world(CoordType px_g, CoordType py_g, CoordType pz_g, CoordType& px_w, CoordType& py_w, CoordType& pz_w) const override;
    virtual bool grid_to_world(CoordType px_g, CoordType py_g, CoordType pz_g, CoordType ps_g, CoordType& px_w, CoordType& py_w, CoordType& pz_w, CoordType& ps_w) const;

    /// print info
    void print(std::ostream& os) const override;

    /// compute four BSpline basis functions
    static bspline_float_type BSpline0(bspline_float_type t)
    {
        return (1-t)*(1-t)*(1-t)/(bspline_float_type)6.0;
    }

    static bspline_float_type BSpline1(bspline_float_type t)
    {
        return (3*t*t*t - 6*t*t + 4)/(bspline_float_type)6.0;
    }

    static bspline_float_type BSpline2(bspline_float_type t)
    {
        return (-3*t*t*t + 3*t*t + 3*t + 1)/(bspline_float_type)6.0;
    }

    static bspline_float_type BSpline3(bspline_float_type t)
    {
        return (t*t*t)/(bspline_float_type)6.0;
    }

    static bspline_float_type BSpline(size_t ind, bspline_float_type t)
    {
        switch (ind)
        {
        case 0:
            return BSpline0(t);
        case 1:
            return BSpline1(t);
        case 2:
            return BSpline2(t);
        case 3:
            return BSpline3(t);
        default:
            throw std::invalid_argument("Index must be smaller than 3");
        }

    }

    /// compute 1st order derivatives of four BSpline basis functions
    static bspline_float_type BSpline0FirstOrderDeriv(bspline_float_type t)
    {
        return -(1-t)*(1-t)/(bspline_float_type)2.0;
    }

    static bspline_float_type BSpline1FirstOrderDeriv(bspline_float_type t)
    {
        return (9*t*t - 12*t)/(bspline_float_type)6.0;
    }

    static bspline_float_type BSpline2FirstOrderDeriv(bspline_float_type t)
    {
        return (-9*t*t + 6*t + 3)/(bspline_float_type)6.0;
    }

    static bspline_float_type BSpline3FirstOrderDeriv(bspline_float_type t)
    {
        return (t*t)/(bspline_float_type)2.0;
    }

    static bspline_float_type BSplineFirstOrderDeriv(size_t ind, bspline_float_type t)
    {
        switch (ind)
        {
        case 0:
            return BSpline0FirstOrderDeriv(t);
        case 1:
            return BSpline1FirstOrderDeriv(t);
        case 2:
            return BSpline2FirstOrderDeriv(t);
        case 3:
            return BSpline3FirstOrderDeriv(t);
        }

        return 0;
    }

    /// compute 2nd order derivatives of four BSpline basis functions
    static bspline_float_type BSpline0SecondOrderDeriv(bspline_float_type t)
    {
        return 1 - t;
    }

    static bspline_float_type BSpline1SecondOrderDeriv(bspline_float_type t)
    {
        return 3*t - 2;
    }

    static bspline_float_type BSpline2SecondOrderDeriv(bspline_float_type t)
    {
        return -3*t + 1;
    }

    static bspline_float_type BSpline3SecondOrderDeriv(bspline_float_type t)
    {
        return t;
    }

    static bspline_float_type BSplineSecondOrderDeriv(size_t ind, bspline_float_type t)
    {
        switch (ind)
        {
        case 0:
            return BSpline0SecondOrderDeriv(t);
        case 1:
            return BSpline1SecondOrderDeriv(t);
        case 2:
            return BSpline2SecondOrderDeriv(t);
        case 3:
            return BSpline3SecondOrderDeriv(t);
        }

        return 0;
    }

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

    /// load the look up table for BSpline functions
    virtual bool loadLookUpTable();

    /// initialize the FFD
    /// define the FFD over a region
    bool initializeBFFD(const PointType& start, const PointType& end, CoordType gridCtrlPtSpacing[DIn]);
    bool initializeBFFD(const PointType& start, const PointType& end, size_t gridCtrlPtNum[DIn]);
    /// define the FFD over the region covered by an image
    bool initializeBFFD(const ImageType& im, const PointType& start, const PointType& end, CoordType gridCtrlPtSpacing[DIn]);
    bool initializeBFFD(const ImageType& im, const PointType& start, const PointType& end, size_t gridCtrlPtNum[DIn]);

    /// 2D
    bool initializeBFFD(const PointType& start, const PointType& end, CoordType dx, CoordType dy);
    bool initializeBFFD(const PointType& start, const PointType& end, size_t sx, size_t sy);
    bool initializeBFFD(const ImageType& im, const PointType& start, const PointType& end, CoordType dx, CoordType dy);
    bool initializeBFFD(const ImageType& im, const PointType& start, const PointType& end, size_t sx, size_t sy);

    /// 3D
    bool initializeBFFD(const PointType& start, const PointType& end, CoordType dx, CoordType dy, CoordType dz);
    bool initializeBFFD(const PointType& start, const PointType& end, size_t sx, size_t sy, size_t sz);
    bool initializeBFFD(const ImageType& im, const PointType& start, const PointType& end, CoordType dx, CoordType dy, CoordType dz);
    bool initializeBFFD(const ImageType& im, const PointType& start, const PointType& end, size_t sx, size_t sy, size_t sz);

    /// 4D
    bool initializeBFFD(const PointType& start, const PointType& end, CoordType dx, CoordType dy, CoordType dz, CoordType ds);
    bool initializeBFFD(const PointType& start, const PointType& end, size_t sx, size_t sy, size_t sz, size_t ss);
    bool initializeBFFD(const ImageType& im, const PointType& start, const PointType& end, CoordType dx, CoordType dy, CoordType dz, CoordType ds);
    bool initializeBFFD(const ImageType& im, const PointType& start, const PointType& end, size_t sx, size_t sy, size_t sz, size_t ss);

    /// look up table for BSpline and its first and second order derivatives
    LUTType LUT_;
    LUTType LUT1_;
    LUTType LUT2_;
};

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
BSplineFFD<T, CoordType, DIn, DOut>::BSplineFFD()
{
    this->loadLookUpTable();
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
BSplineFFD<T, CoordType, DIn, DOut>::~BSplineFFD()
{
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
bool BSplineFFD<T, CoordType, DIn, DOut>::loadLookUpTable()
{
    try
    {
        long long ii;
        double gapInLUT = (double)(BSPLINELUTSIZE-1);

        #pragma omp parallel for default(none) private(ii) shared(gapInLUT)
        for (ii=0; ii<(long long)BSPLINELUTSIZE; ii++)
        {
            bspline_float_type g = (bspline_float_type)(ii/gapInLUT);

            LUT_[ii][0]   = BSpline0(g);
            LUT_[ii][1]   = BSpline1(g);
            LUT_[ii][2]   = BSpline2(g);
            LUT_[ii][3]   = BSpline3(g);

            LUT1_[ii][0]  = BSpline0FirstOrderDeriv(g);
            LUT1_[ii][1]  = BSpline1FirstOrderDeriv(g);
            LUT1_[ii][2]  = BSpline2FirstOrderDeriv(g);
            LUT1_[ii][3]  = BSpline3FirstOrderDeriv(g);

            LUT2_[ii][0]  = BSpline0SecondOrderDeriv(g);
            LUT2_[ii][1]  = BSpline1SecondOrderDeriv(g);
            LUT2_[ii][2]  = BSpline2SecondOrderDeriv(g);
            LUT2_[ii][3]  = BSpline3SecondOrderDeriv(g);
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors happened in BSplineFFD<T, CoordType, DIn, DOut>::loadLookUpTable() ...");
        return false;
    }

    return true;
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool BSplineFFD<T, CoordType, DIn, DOut>::world_to_grid(const CoordType pt_w[D], CoordType pt_g[D]) const
{
    try
    {
        this->ctrl_pt_[0].world_to_image(pt_w, pt_g);
        unsigned int d;
        for ( d=0; d<D; d++ )
        {
            pt_g[d] -= BSPLINEPADDINGSIZE;
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors happened in world_to_grid(const CoordType pt_w[D], CoordType pt_g[D]) const ... ");
        return false;
    }

    return true;
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool BSplineFFD<T, CoordType, DIn, DOut>::world_to_grid(CoordType px_w, CoordType py_w, CoordType& px_g, CoordType& py_g) const
{
    GADGET_CHECK_RETURN_FALSE(DIn==2);

    try
    {
        this->ctrl_pt_[0].world_to_image(px_w, py_w, px_g, py_g);
        px_g -= BSPLINEPADDINGSIZE;
        py_g -= BSPLINEPADDINGSIZE;
    }
    catch(...)
    {
        GERROR_STREAM("Errors happened in world_to_grid(CoordType px_w, CoordType py_w, CoordType& px_g, CoordType& py_g) const ... ");
        return false;
    }

    return true;
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool BSplineFFD<T, CoordType, DIn, DOut>::world_to_grid(CoordType px_w, CoordType py_w, CoordType pz_w, CoordType& px_g, CoordType& py_g, CoordType& pz_g) const
{
    GADGET_CHECK_RETURN_FALSE(DIn==3);

    try
    {
        this->ctrl_pt_[0].world_to_image(px_w, py_w, pz_w, px_g, py_g, pz_g);
        px_g -= BSPLINEPADDINGSIZE;
        py_g -= BSPLINEPADDINGSIZE;
        pz_g -= BSPLINEPADDINGSIZE;
    }
    catch(...)
    {
        GERROR_STREAM("Errors happened in world_to_grid(CoordType px_w, CoordType py_w, CoordType pz_w, CoordType& px_g, CoordType& py_g, CoordType& pz_g) const ... ");
        return false;
    }

    return true;
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool BSplineFFD<T, CoordType, DIn, DOut>::world_to_grid(CoordType px_w, CoordType py_w, CoordType pz_w, CoordType ps_w, CoordType& px_g, CoordType& py_g, CoordType& pz_g, CoordType& ps_g) const
{
    GADGET_CHECK_RETURN_FALSE(DIn==4);

    try
    {
        this->ctrl_pt_[0].world_to_image(px_w, py_w, pz_w, ps_w, px_g, py_g, pz_g, ps_g);
        px_g -= BSPLINEPADDINGSIZE;
        py_g -= BSPLINEPADDINGSIZE;
        pz_g -= BSPLINEPADDINGSIZE;
        ps_g -= BSPLINEPADDINGSIZE;
    }
    catch(...)
    {
        GERROR_STREAM("Errors happened in world_to_grid(CoordType px_w, CoordType py_w, CoordType pz_w, CoordType ps_w, CoordType& px_g, CoordType& py_g, CoordType& pz_g, CoordType& ps_g) const ... ");
        return false;
    }

    return true;
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool BSplineFFD<T, CoordType, DIn, DOut>::grid_to_world(const CoordType pt_g[D], CoordType pt_w[D]) const
{
    try
    {
        CoordType pt_g_padded[D];
        unsigned int d;
        for ( d=0; d<D; d++ )
        {
            pt_g_padded[d] = pt_g[d] + BSPLINEPADDINGSIZE;
        }

        this->ctrl_pt_[0].image_to_world(pt_g_padded, pt_w);
    }
    catch(...)
    {
        GERROR_STREAM("Errors happened in grid_to_world(const CoordType pt_g[D], CoordType pt_w[D]) const ... ");
        return false;
    }

    return true;
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool BSplineFFD<T, CoordType, DIn, DOut>::grid_to_world(CoordType px_g, CoordType py_g, CoordType& px_w, CoordType& py_w) const
{
    GADGET_CHECK_RETURN_FALSE(DIn==2);

    try
    {
        px_g += BSPLINEPADDINGSIZE;
        py_g += BSPLINEPADDINGSIZE;
        this->ctrl_pt_[0].image_to_world(px_g, py_g, px_w, py_w);
    }
    catch(...)
    {
        GERROR_STREAM("Errors happened in grid_to_world(CoordType px_g, CoordType py_g, CoordType& px_w, CoordType& py_w) const ... ");
        return false;
    }

    return true;
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool BSplineFFD<T, CoordType, DIn, DOut>::grid_to_world(CoordType px_g, CoordType py_g, CoordType pz_g, CoordType& px_w, CoordType& py_w, CoordType& pz_w) const
{
    GADGET_CHECK_RETURN_FALSE(DIn==3);

    try
    {
        px_g += BSPLINEPADDINGSIZE;
        py_g += BSPLINEPADDINGSIZE;
        pz_g += BSPLINEPADDINGSIZE;
        this->ctrl_pt_[0].image_to_world(px_g, py_g, pz_g, px_w, py_w, pz_w);
    }
    catch(...)
    {
        GERROR_STREAM("Errors happened in grid_to_world(CoordType px_g, CoordType py_g, CoordType pz_g, CoordType& px_w, CoordType& py_w, CoordType& pz_w) const ... ");
        return false;
    }

    return true;
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool BSplineFFD<T, CoordType, DIn, DOut>::grid_to_world(CoordType px_g, CoordType py_g, CoordType pz_g, CoordType ps_g, CoordType& px_w, CoordType& py_w, CoordType& pz_w, CoordType& ps_w) const
{
    GADGET_CHECK_RETURN_FALSE(DIn==4);

    try
    {
        px_g += BSPLINEPADDINGSIZE;
        py_g += BSPLINEPADDINGSIZE;
        pz_g += BSPLINEPADDINGSIZE;
        ps_g += BSPLINEPADDINGSIZE;
        this->ctrl_pt_[0].image_to_world(px_g, py_g, pz_g, ps_g, px_w, py_w, pz_w, ps_w);
    }
    catch(...)
    {
        GERROR_STREAM("Errors happened in grid_to_world(CoordType px_g, CoordType py_g, CoordType pz_g, CoordType ps_g, CoordType& px_w, CoordType& py_w, CoordType& pz_w, CoordType& ps_w) const ... ");
        return false;
    }

    return true;
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool BSplineFFD<T, CoordType, DIn, DOut>::initializeBFFD(const PointType& start, const PointType& end, CoordType gridCtrlPtSpacing[DIn])
{
    try
    {
        unsigned int d;
        for ( d=0; d<DIn; d++ )
        {
            GADGET_CHECK_RETURN_FALSE(end(d) > start(d));
        }

        std::vector<size_t> dim(DIn, 2);
        std::vector<coord_type> pixelSize(DIn, 1);
        std::vector<coord_type> origin(DIn, 0);

        for ( d=0; d<DIn; d++ )
        {
            dim[d] = FFD_MKINT( (end(d)-start(d))/gridCtrlPtSpacing[d] ) + 1;
            pixelSize[d] = (end(d)-start(d))/(dim[d]-1);

            /// add the padding
            dim[d] += 2*BSPLINEPADDINGSIZE;

            origin[d] = -pixelSize[d]*BSPLINEPADDINGSIZE;
        }

        for ( d=0; d<DOut; d++ )
        {
            this->ctrl_pt_[d].create(dim, pixelSize, origin);
            Gadgetron::clear(this->ctrl_pt_[d]);
        }
    }
    catch(...)
    {
        GERROR_STREAM("Error happened in initializeBFFD(const PointType& start, const PointType& end, CoordType gridCtrlPtSpacing[DIn]) ... ");
        return false;
    }

    return true;
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool BSplineFFD<T, CoordType, DIn, DOut>::initializeBFFD(const PointType& start, const PointType& end, size_t gridCtrlPtNum[DIn])
{
    try
    {
        unsigned int d;
        for ( d=0; d<DIn; d++ )
        {
            GADGET_CHECK_RETURN_FALSE(end(d) > start(d));
        }

        std::vector<size_t> dim(DIn, 2);
        std::vector<coord_type> pixelSize(DIn, 1);
        std::vector<coord_type> origin(DIn, 0);

        for ( d=0; d<DIn; d++ )
        {
            dim[d] = gridCtrlPtNum[d];
            if ( dim[d] < 3 ) dim[d] = 3;

            pixelSize[d] = (end(d)-start(d))/(dim[d]-1);

            /// add the padding
            dim[d] += 2*BSPLINEPADDINGSIZE;

            origin[d] = -pixelSize[d]*BSPLINEPADDINGSIZE;
        }

        for ( d=0; d<DOut; d++ )
        {
            this->ctrl_pt_[d].create(dim, pixelSize, origin);
            Gadgetron::clear(this->ctrl_pt_[d]);
        }
    }
    catch(...)
    {
        GERROR_STREAM("Error happened in initializeBFFD(const PointType& start, const PointType& end, CoordType gridCtrlPtNum[DIn]) ... ");
        return false;
    }

    return true;
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool BSplineFFD<T, CoordType, DIn, DOut>::initializeBFFD(const ImageType& im, const PointType& start, const PointType& end, CoordType gridCtrlPtSpacing[DIn])
{
    try
    {
        unsigned int d;
        for ( d=0; d<DIn; d++ )
        {
            GADGET_CHECK_RETURN_FALSE(end(d) > start(d));
        }

        std::vector<size_t> dim(DIn, 2);
        std::vector<coord_type> pixelSize(DIn, 1);
        std::vector<coord_type> origin(DIn, 0);

        std::vector<coord_type> firstCtrlPt(DIn);

        for ( d=0; d<DIn; d++ )
        {
            dim[d] = FFD_MKINT( (end(d)-start(d))/gridCtrlPtSpacing[d] ) + 1;
            pixelSize[d] = (end(d)-start(d))/(dim[d]-1);

            /// add the padding
            dim[d] += 2*BSPLINEPADDINGSIZE;

            firstCtrlPt[d] = -pixelSize[d]*BSPLINEPADDINGSIZE/im.get_pixel_size(d);
        }
        im.image_to_world( firstCtrlPt, origin);

        for ( d=0; d<DOut; d++ )
        {
            this->ctrl_pt_[d].create(dim, pixelSize, origin);
            Gadgetron::clear(this->ctrl_pt_[d]);
        }
    }
    catch(...)
    {
        GERROR_STREAM("Error happened in initializeBFFD(const ImageType& im, const PointType& start, const PointType& end, CoordType gridCtrlPtSpacing[DIn]) ... ");
        return false;
    }

    return true;
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool BSplineFFD<T, CoordType, DIn, DOut>::initializeBFFD(const ImageType& im, const PointType& start, const PointType& end, size_t gridCtrlPtNum[DIn])
{
    try
    {
        unsigned int d;
        for ( d=0; d<DIn; d++ )
        {
            GADGET_CHECK_RETURN_FALSE(end(d) > start(d));
        }

        std::vector<size_t> dim(DIn, 2);
        std::vector<coord_type> pixelSize(DIn, 1);
        std::vector<coord_type> origin(DIn, 0);

        std::vector<coord_type> firstCtrlPt(DIn);

        for ( d=0; d<DIn; d++ )
        {
            dim[d] = gridCtrlPtNum[d];
            if ( dim[d] < 3 ) dim[d] = 3;

            pixelSize[d] = (end(d)-start(d))/(dim[d]-1);

            /// add the padding
            dim[d] += 2*BSPLINEPADDINGSIZE;

            firstCtrlPt[d] = -pixelSize[d]*BSPLINEPADDINGSIZE/im.get_pixel_size(d);
        }
        im.image_to_world( firstCtrlPt, origin);

        for ( d=0; d<DOut; d++ )
        {
            this->ctrl_pt_[d].create(dim, pixelSize, origin);
            Gadgetron::clear(this->ctrl_pt_[d]);
        }
    }
    catch(...)
    {
        GERROR_STREAM("Error happened in initializeBFFD(const ImageType& im, const PointType& start, const PointType& end, CoordType gridCtrlPtNum[DIn]) ... ");
        return false;
    }

    return true;
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool BSplineFFD<T, CoordType, DIn, DOut>::initializeBFFD(const PointType& start, const PointType& end, CoordType dx, CoordType dy)
{
    CoordType gridCtrlPtSpacing[2];
    gridCtrlPtSpacing[0] = dx;
    gridCtrlPtSpacing[1] = dy;
    return this->initializeBFFD(start, end, gridCtrlPtSpacing);
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool BSplineFFD<T, CoordType, DIn, DOut>::initializeBFFD(const PointType& start, const PointType& end, size_t sx, size_t sy)
{
    size_t gridCtrlPtNum[2];
    gridCtrlPtNum[0] = sx;
    gridCtrlPtNum[1] = sy;
    return this->initializeBFFD(start, end, gridCtrlPtNum);
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool BSplineFFD<T, CoordType, DIn, DOut>::initializeBFFD(const ImageType& im, const PointType& start, const PointType& end, CoordType dx, CoordType dy)
{
    CoordType gridCtrlPtSpacing[2];
    gridCtrlPtSpacing[0] = dx;
    gridCtrlPtSpacing[1] = dy;
    return this->initializeBFFD(im, start, end, gridCtrlPtSpacing);
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool BSplineFFD<T, CoordType, DIn, DOut>::initializeBFFD(const ImageType& im, const PointType& start, const PointType& end, size_t sx, size_t sy)
{
    size_t gridCtrlPtNum[2];
    gridCtrlPtNum[0] = sx;
    gridCtrlPtNum[1] = sy;
    return this->initializeBFFD(im, start, end, gridCtrlPtNum);
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool BSplineFFD<T, CoordType, DIn, DOut>::initializeBFFD(const PointType& start, const PointType& end, CoordType dx, CoordType dy, CoordType dz)
{
    CoordType gridCtrlPtSpacing[3];
    gridCtrlPtSpacing[0] = dx;
    gridCtrlPtSpacing[1] = dy;
    gridCtrlPtSpacing[2] = dz;
    return this->initializeBFFD(start, end, gridCtrlPtSpacing);
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool BSplineFFD<T, CoordType, DIn, DOut>::initializeBFFD(const PointType& start, const PointType& end, size_t sx, size_t sy, size_t sz)
{
    size_t gridCtrlPtNum[3];
    gridCtrlPtNum[0] = sx;
    gridCtrlPtNum[1] = sy;
    gridCtrlPtNum[2] = sz;
    return this->initializeBFFD(start, end, gridCtrlPtNum);
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool BSplineFFD<T, CoordType, DIn, DOut>::initializeBFFD(const ImageType& im, const PointType& start, const PointType& end, CoordType dx, CoordType dy, CoordType dz)
{
    CoordType gridCtrlPtSpacing[3];
    gridCtrlPtSpacing[0] = dx;
    gridCtrlPtSpacing[1] = dy;
    gridCtrlPtSpacing[2] = dz;
    return this->initializeBFFD(im, start, end, gridCtrlPtSpacing);
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool BSplineFFD<T, CoordType, DIn, DOut>::initializeBFFD(const ImageType& im, const PointType& start, const PointType& end, size_t sx, size_t sy, size_t sz)
{
    size_t gridCtrlPtNum[3];
    gridCtrlPtNum[0] = sx;
    gridCtrlPtNum[1] = sy;
    gridCtrlPtNum[2] = sz;
    return this->initializeBFFD(im, start, end, gridCtrlPtNum);
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool BSplineFFD<T, CoordType, DIn, DOut>::initializeBFFD(const PointType& start, const PointType& end, CoordType dx, CoordType dy, CoordType dz, CoordType ds)
{
    CoordType gridCtrlPtSpacing[4];
    gridCtrlPtSpacing[0] = dx;
    gridCtrlPtSpacing[1] = dy;
    gridCtrlPtSpacing[2] = dz;
    gridCtrlPtSpacing[3] = ds;
    return this->initializeBFFD(start, end, gridCtrlPtSpacing);
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool BSplineFFD<T, CoordType, DIn, DOut>::initializeBFFD(const PointType& start, const PointType& end, size_t sx, size_t sy, size_t sz, size_t ss)
{
    size_t gridCtrlPtNum[4];
    gridCtrlPtNum[0] = sx;
    gridCtrlPtNum[1] = sy;
    gridCtrlPtNum[2] = sz;
    gridCtrlPtNum[3] = ss;
    return this->initializeBFFD(start, end, gridCtrlPtNum);
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool BSplineFFD<T, CoordType, DIn, DOut>::initializeBFFD(const ImageType& im, const PointType& start, const PointType& end, CoordType dx, CoordType dy, CoordType dz, CoordType ds)
{
    CoordType gridCtrlPtSpacing[4];
    gridCtrlPtSpacing[0] = dx;
    gridCtrlPtSpacing[1] = dy;
    gridCtrlPtSpacing[2] = dz;
    gridCtrlPtSpacing[3] = ds;
    return this->initializeBFFD(im, start, end, gridCtrlPtSpacing);
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool BSplineFFD<T, CoordType, DIn, DOut>::initializeBFFD(const ImageType& im, const PointType& start, const PointType& end, size_t sx, size_t sy, size_t sz, size_t ss)
{
    size_t gridCtrlPtNum[4];
    gridCtrlPtNum[0] = sx;
    gridCtrlPtNum[1] = sy;
    gridCtrlPtNum[2] = sz;
    gridCtrlPtNum[3] = ss;
    return this->initializeBFFD(im, start, end, gridCtrlPtNum);
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
void BSplineFFD<T, CoordType, DIn, DOut>::print(std::ostream& os) const
{
    using namespace std;

    os << "---------------------- BSpline Free Form Deformation ------------------" << endl;
    os << "Define the interface for BSpline Free Form Deformation (BFFD) " << endl;
    os << "------------------------------------------------------------------------------" << endl;
}

}
