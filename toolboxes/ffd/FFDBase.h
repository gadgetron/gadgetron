/** \file       FFDBase.h
    \brief      Base class for FreeFormDeformation package

                FreeFormDeformation (FFD) is a general purpose scatter interpolation algorithm. It is widely used in numerical applications, 
                such as image registration, data inteprolation and geometric modelling etc.

                [1] http://en.wikipedia.org/wiki/Free-form_deformation

                [2] Seungyong Lee ; Dept. of Comput. Sci., Pohang Inst. of Sci. & Technol., South Korea ; Wolberg, G. ; Sung Yong Shin. Scattered data interpolation with multilevel B-splines. IEEE 
                    Transactions on Visualization and Computer Graphics, Volume 3, Issue 3, 1997.

                [3] D Rueckert, LI Sonoda, C Hayes, DLG Hill, MO Leach, DJ Hawkes. Nonrigid registration using free-form deformations: application to breast MR images. IEEE 
                    Transactions on Medical Imaging, Volume 18, Issue 8, 1999.

    \author     Hui Xue
*/

#pragma once

#include <typeinfo>
#include <cmath>
#include <complex>
#include "hoNDArray.h"
#include "hoNDImage.h"
#include "GadgetronTimer.h"

// #include "gtPlusIOAnalyze.h"

#ifdef USE_OMP
    #include "omp.h"
#endif // USE_OMP

#define FFD_MKINT(a) (((a)>=0)?((int)((a)+0.5)):((int)((a)-0.5)))

namespace Gadgetron { 

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut>
class FFDBase
{
public:

    typedef FFDBase<T, CoordType, DIn, DOut> Self;

    typedef typename realType<T>::Type real_value_type;

    typedef CoordType coord_type;

    enum { D = DIn };

    /// array to store the coordinates of spatial points
    /// has the dimension of DIn by N for N points
    typedef hoNDArray<CoordType> CoordArrayType;

    /// array to store the point value
    /// for N points, the dimension of array is DOut by N
    /// DOut is equal or larger than 1; if larger than 1, the 
    /// vectorized FFD is computed
    typedef hoNDArray<T> ValueArrayType;
    typedef ValueArrayType ArrayType;

    typedef hoNDArray<float> MaskArrayType;

    /// control point grip type
    typedef hoNDImage<T, DIn> FFDCtrlPtGridType;

    /// point type
    typedef hoNDPoint<CoordType, DIn> PointType;

    /// image type
    typedef hoNDImage<T, DIn> ImageType;

    FFDBase();
    virtual ~FFDBase();

    /// evaluate the FFD at a grid location
    /// the input points are in the FFD grid
    virtual bool evaluateFFD(const CoordType pt[D], T r[DOut]) const = 0;
    virtual bool evaluateFFD(const CoordType* pt[D], T* r[DOut], size_t N) const;
    virtual bool evaluateFFD(const PointType& pt, T r[DOut]) const;
    virtual bool evaluateFFDArray(const CoordArrayType& pts, ValueArrayType& r) const;

    /// evaluate the 1st order derivative of FFD at a grid location
    /// deriv: derivative for all D dimensions and all DOut values
    virtual bool evaluateFFDDerivative(const CoordType pt[D], T deriv[D][DOut]) const = 0;
    virtual bool evaluateFFDDerivative(const PointType& pt, T deriv[D][DOut]) const;

    virtual bool evaluateFFDDX(const CoordType pt[D], T dx[DOut]) const;
    virtual bool evaluateFFDDY(const CoordType pt[D], T dy[DOut]) const;
    virtual bool evaluateFFDDZ(const CoordType pt[D], T dz[DOut]) const;
    virtual bool evaluateFFDDS(const CoordType pt[D], T ds[DOut]) const;

    /// calculate the 1st order derivative of FFD at a world coordinate location with the world coordinate unit
    virtual bool evaluateWorldDerivative(const CoordType pt[D], T deriv[D][DOut]) const;

    virtual bool evaluateWorldDX(const CoordType pt[D], T dx[DOut]) const;
    virtual bool evaluateWorldDY(const CoordType pt[D], T dy[DOut]) const;
    virtual bool evaluateWorldDZ(const CoordType pt[D], T dz[DOut]) const;
    virtual bool evaluateWorldDS(const CoordType pt[D], T ds[DOut]) const;

    /// evaluate the 2nd order derivative of FFD at a grid location
    /// dderiv : D*D vector, stores dxx dxy dxz ...; dyx dyy dyz ...; dzx dzy dzz ...
    virtual bool evaluateFFDSecondOrderDerivative(const CoordType pt[D], T dderiv[D*D][DOut]) const = 0;
    virtual bool evaluateFFDSecondOrderDerivative(const PointType& pt, T dderiv[D*D][DOut]) const;

    /// evaluate the FFD at a world location
    virtual bool evaluateFFDW(const CoordType pt[D], T r[DOut]) const;
    virtual bool evaluateFFDW(CoordType px, CoordType py, T r[DOut]) const;
    virtual bool evaluateFFDW(CoordType px, CoordType py, CoordType pz, T r[DOut]) const;
    virtual bool evaluateFFDW(CoordType px, CoordType py, CoordType pz, CoordType ps, T r[DOut]) const;

    virtual bool evaluateFFDDerivativeW(const CoordType pt[D], T deriv[D][DOut]) const;
    virtual bool evaluateFFDDerivativeW(CoordType px, CoordType py, T deriv[D][DOut]) const;
    virtual bool evaluateFFDDerivativeW(CoordType px, CoordType py, CoordType pz, T deriv[D][DOut]) const;
    virtual bool evaluateFFDDerivativeW(CoordType px, CoordType py, CoordType pz, CoordType ps, T deriv[D][DOut]) const;

    virtual bool evaluateFFDSecondOrderDerivativeW(const CoordType pt[D], T dderiv[D*D][DOut]) const;
    virtual bool evaluateFFDSecondOrderDerivativeW(CoordType px, CoordType py, T dderiv[D*D][DOut]) const;
    virtual bool evaluateFFDSecondOrderDerivativeW(CoordType px, CoordType py, CoordType pz, T dderiv[D*D][DOut]) const;
    virtual bool evaluateFFDSecondOrderDerivativeW(CoordType px, CoordType py, CoordType pz, CoordType ps, T dderiv[D*D][DOut]) const;

    /// compute the FFD approximation once
    /// pos : the position of input points, DIn by N
    /// value : the value on input points, DOut by N
    /// residual : the approximation residual after computing FFD, DOut by N
    /// N : the number of points
    virtual bool ffdApprox(const CoordArrayType& pos, ValueArrayType& value, ValueArrayType& residual, real_value_type& totalResidual, size_t N) = 0;

    /// compute the FFD approximation with refinement, see ref [2]
    /// numOfRefinement : number of grid refinement
    virtual bool ffdApprox(const CoordArrayType& pos, ValueArrayType& value, ValueArrayType& residual, real_value_type& totalResidual, size_t N, size_t numOfRefinement);

    /// keep refine the FFD until either the maximal refinement level is reached or total residual is less than a threshold
    virtual bool ffdApprox(const CoordArrayType& pos, ValueArrayType& value, ValueArrayType& residual, real_value_type& totalResidual, size_t N, size_t& numOfRefinement, real_value_type thresResidual, size_t maxNumOfRefinement);

    //// fft approximation with input in the world coordinates
    virtual bool ffdApproxW(const CoordArrayType& pos, ValueArrayType& value, ValueArrayType& residual, real_value_type& totalResidual, size_t N, size_t numOfRefinement);
    virtual bool ffdApproxW(const CoordArrayType& pos, ValueArrayType& value, ValueArrayType& residual, real_value_type& totalResidual, size_t N, size_t& numOfRefinement, real_value_type thresResidual, size_t maxNumOfRefinement);

    /// easy-to-use function calls for image and array

    /// convert every pixel in the image to FFD point inputs with world coordiantes
    virtual bool imageToFFDInputsW(ImageType target[DOut], CoordArrayType& pos, ValueArrayType& value);
    /// mask == 0 means this point is excluded from approximation
    virtual bool imageToFFDInputsW(ImageType target[DOut], const MaskArrayType& mask, CoordArrayType& pos, ValueArrayType& value);

    /// convert every pixel in the array to FFD point inputs with world coordiantes
    virtual bool arrayToFFDInputsW(ArrayType target[DOut], CoordArrayType& pos, ValueArrayType& value);
    /// mask == 0 means this point is excluded from approximation
    virtual bool arrayToFFDInputsW(ArrayType target[DOut], const MaskArrayType& mask, CoordArrayType& pos, ValueArrayType& value);

    /// for Image type
    virtual bool ffdApproxImage(ImageType target[DOut], real_value_type& totalResidual, size_t numOfRefinement);
    virtual bool ffdApproxImage(ImageType target[DOut], real_value_type& totalResidual, size_t& numOfRefinement, real_value_type thresResidual, size_t maxNumOfRefinement);
    virtual bool ffdApproxImage(ImageType& target, real_value_type& totalResidual, size_t numOfRefinement);
    virtual bool ffdApproxImage(ImageType& target, real_value_type& totalResidual, size_t& numOfRefinement, real_value_type thresResidual, size_t maxNumOfRefinement);

    virtual bool ffdApproxImage(ImageType target[DOut], const MaskArrayType& mask, real_value_type& totalResidual, size_t numOfRefinement);
    virtual bool ffdApproxImage(ImageType target[DOut], const MaskArrayType& mask, real_value_type& totalResidual, size_t& numOfRefinement, real_value_type thresResidual, size_t maxNumOfRefinement);

    /// for Array type
    virtual bool ffdApproxArray(ArrayType target[DOut], real_value_type& totalResidual, size_t numOfRefinement);
    virtual bool ffdApproxArray(ArrayType target[DOut], real_value_type& totalResidual, size_t& numOfRefinement, real_value_type thresResidual, size_t maxNumOfRefinement);
    virtual bool ffdApproxArray(ArrayType& target, real_value_type& totalResidual, size_t numOfRefinement);
    virtual bool ffdApproxArray(ArrayType& target, real_value_type& totalResidual, size_t& numOfRefinement, real_value_type thresResidual, size_t maxNumOfRefinement);

    virtual bool ffdApproxArray(ArrayType target[DOut], const MaskArrayType& mask, real_value_type& totalResidual, size_t numOfRefinement);
    virtual bool ffdApproxArray(ArrayType target[DOut], const MaskArrayType& mask, real_value_type& totalResidual, size_t& numOfRefinement, real_value_type thresResidual, size_t maxNumOfRefinement);

    /// As suggested in ref [2], the BSpline FFD can be refined to achieve better approximation
    virtual bool refine() = 0;

    /// utility functions for easy-to-use

    /// get control points
    FFDCtrlPtGridType& get_ctrl_pt(unsigned int d) { return this->ctrl_pt_[d]; }
    const FFDCtrlPtGridType& get_ctrl_pt(unsigned int d) const { return this->ctrl_pt_[d]; }

    /// get the size of control point arrays
    virtual size_t get_size(size_t dimension) const { return ctrl_pt_[0].get_size(dimension); }
    virtual std::vector<size_t> dimensions() const { return ctrl_pt_[0].dimensions(); }

    /// get the spacing of of control point arrays
    virtual coord_type get_spacing(size_t dimension) const { return ctrl_pt_[0].get_pixel_size(dimension); }
    virtual void get_spacing(std::vector<coord_type>& spacing) const { ctrl_pt_[0].get_pixel_size(spacing); }

    /// get/set a control point value
    virtual T get(size_t x, size_t y, size_t d) const { return ctrl_pt_[d](x, y); }
    virtual void set(size_t x, size_t y, size_t d, T v) { ctrl_pt_[d](x, y) = v; }

    virtual T get(size_t x, size_t y, size_t z, size_t d) const { return ctrl_pt_[d](x, y, z); }
    virtual void set(size_t x, size_t y, size_t z, size_t d, T v) { ctrl_pt_[d](x, y, z) = v; }

    virtual T get(size_t x, size_t y, size_t z, size_t s, size_t d) const { return ctrl_pt_[d](x, y, z, s); }
    virtual void set(size_t x, size_t y, size_t z, size_t s, size_t d, T v) { ctrl_pt_[d](x, y, z, s) = v; }

    /// offset to/from indexes for control points
    virtual size_t calculate_offset(size_t x, size_t y) const { return ctrl_pt_[0].calculate_offset(x, y); }
    virtual void calculate_index( size_t offset, size_t& x, size_t& y ) const { ctrl_pt_[0].calculate_index(offset, x, y); }

    virtual size_t calculate_offset(size_t x, size_t y, size_t z) const { return ctrl_pt_[0].calculate_offset(x, y, z); }
    virtual void calculate_index( size_t offset, size_t& x, size_t& y, size_t& z ) const { ctrl_pt_[0].calculate_index(offset, x, y, z); }

    virtual size_t calculate_offset(size_t x, size_t y, size_t z, size_t s) const { return ctrl_pt_[0].calculate_offset(x, y, z, s); }
    virtual void calculate_index( size_t offset, size_t& x, size_t& y, size_t& z, size_t& s ) const { ctrl_pt_[0].calculate_index(offset, x, y, z, s); }

    /// compute the control point location in world coordinates
    virtual void get_location(size_t x, size_t y, CoordType& sx, CoordType& sy) const { ctrl_pt_[0].image_to_world(x, y, sx, sy); }
    virtual void get_location(size_t x, size_t y, size_t z, CoordType& sx, CoordType& sy, CoordType& sz) const { ctrl_pt_[0].image_to_world(x, y, z, sx, sy, sz); }
    virtual void get_location(size_t x, size_t y, size_t z, size_t s, CoordType& sx, CoordType& sy, CoordType& sz, CoordType& ss) const { ctrl_pt_[0].image_to_world(x, y, z, s, sx, sy, sz, ss); }

    /// convert a world coordinate point to FFD grid location
    virtual bool world_to_grid(const CoordArrayType& pt_w, CoordArrayType& pt_g) const;
    virtual bool world_to_grid(const CoordType pt_w[D], CoordType pt_g[D]) const;
    virtual bool world_to_grid(CoordType px_w, CoordType py_w, CoordType& px_g, CoordType& py_g) const;
    virtual bool world_to_grid(CoordType px_w, CoordType py_w, CoordType pz_w, CoordType& px_g, CoordType& py_g, CoordType& pz_g) const;
    virtual bool world_to_grid(CoordType px_w, CoordType py_w, CoordType pz_w, CoordType ps_w, CoordType& px_g, CoordType& py_g, CoordType& pz_g, CoordType& ps_g) const;

    virtual bool grid_to_world(const CoordArrayType& pt_g, CoordArrayType& pt_w) const;
    virtual bool grid_to_world(const CoordType pt_g[D], CoordType pt_w[D]) const;
    virtual bool grid_to_world(CoordType px_g, CoordType py_g, CoordType& px_w, CoordType& py_w) const;
    virtual bool grid_to_world(CoordType px_g, CoordType py_g, CoordType pz_g, CoordType& px_w, CoordType& py_w, CoordType& pz_w) const;

    /// evaluate FFD for every pixel in the target image
    /// the image pixel will first be converted to world-coordinate
    /// and then converted to FFD grid location
    virtual bool evaluateFFDOnImage(ImageType& target) const;
    virtual bool evaluateFFDOnImage(ImageType target[DOut]) const;

    /// evaluate FFD for every elements in an array
    /// the point indexes will be taken as the FFD grid location
    virtual bool evaluateFFDOnArray(hoNDArray<T>& target) const;
    virtual bool evaluateFFDOnArray(hoNDArray<T> target[DOut]) const;

    /// clear the control points
    virtual bool clear(T v=0);

    /// print info
    virtual void print(std::ostream& os) const;

    /// whether to perform timing and print out messages
    bool performTiming_;

    /// debug folder
    std::string debugFolder_;

protected:

    /// control points
    FFDCtrlPtGridType ctrl_pt_[DOut];

    /// clock for timing
    //Gadgetron::GadgetronTimer gt_timer1_;
    //Gadgetron::GadgetronTimer gt_timer2_;
    //Gadgetron::GadgetronTimer gt_timer3_;

    /// exporter
    // Gadgetron::gtPlus::gtPlusIOAnalyze gt_exporter_;
};

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
FFDBase<T, CoordType, DIn, DOut>::FFDBase()
{
    //gt_timer1_.set_timing_in_destruction(false);
    //gt_timer2_.set_timing_in_destruction(false);
    //gt_timer3_.set_timing_in_destruction(false);
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
FFDBase<T, CoordType, DIn, DOut>::~FFDBase()
{
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool FFDBase<T, CoordType, DIn, DOut>::evaluateFFD(const CoordType* pt[D], T* r[DOut], size_t N) const
{
    try
    {
        long long n;
#pragma omp parallel for private(n) shared(N, pt, r)
        for ( n=0; n<(long long)N; n++ )
        {
            this->evaluateFFD(pt[n], r[n]);
        }
    }
    catch(...)
    {
        GERROR_STREAM("Error happened in evaluateFFD(const CoordType* pt[D], T* r[DOut], size_t N) const ... ");
        return false;
    }

    return true;
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool FFDBase<T, CoordType, DIn, DOut>::evaluateFFD(const PointType& pt, T r[DOut]) const
{
    GADGET_CHECK_RETURN_FALSE(this->evaluateFFD(pt.begin(), r));
    return true;
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool FFDBase<T, CoordType, DIn, DOut>::evaluateFFDArray(const CoordArrayType& pts, ValueArrayType& r) const
{
    try
    {
        size_t N = pts.get_size(1);
        GADGET_CHECK_RETURN_FALSE(pts.get_size(0)==DIn);

        if ( r.get_size(1)!=N || r.get_size(0)!=DOut )
        {
            r.create(DOut, N);
        }

        const CoordType* pPts = pts.begin();
        T* pR = r.begin();

        long long n;
#pragma omp parallel for private(n) shared(N, pPts, pR, r)
        for ( n=0; n<(long long)N; n++ )
        {
            this->evaluateFFD(pPts+n*DIn, pR+n*DOut);
        }
    }
    catch(...)
    {
        GERROR_STREAM("Error happened in evaluateFFDArray(const CoordArrayType& pts, ValueArrayType& r) const ... ");
        return false;
    }

    return true;
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool FFDBase<T, CoordType, DIn, DOut>::evaluateFFDDerivative(const PointType& pt, T deriv[D][DOut]) const
{
    GADGET_CHECK_RETURN_FALSE(this->evaluateFFDDerivative(pt.begin(), deriv));
    return true;
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool FFDBase<T, CoordType, DIn, DOut>::evaluateFFDDX(const CoordType pt[D], T dx[DOut]) const
{
    T deriv[D][DOut];
    GADGET_CHECK_RETURN_FALSE(this->evaluateFFDDerivative(pt, deriv));
    memcpy(dx, deriv, sizeof(T)*DOut);
    return true;
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool FFDBase<T, CoordType, DIn, DOut>::evaluateFFDDY(const CoordType pt[D], T dy[DOut]) const
{
    T deriv[D][DOut];
    GADGET_CHECK_RETURN_FALSE(this->evaluateFFDDerivative(pt, deriv));
    memcpy(dy, deriv+sizeof(T)*DOut, sizeof(T)*DOut);
    return true;
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool FFDBase<T, CoordType, DIn, DOut>::evaluateFFDDZ(const CoordType pt[D], T dz[DOut]) const
{
    T deriv[D][DOut];
    GADGET_CHECK_RETURN_FALSE(this->evaluateFFDDerivative(pt, deriv));
    memcpy(dz, deriv+2*sizeof(T)*DOut, sizeof(T)*DOut);
    return true;
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool FFDBase<T, CoordType, DIn, DOut>::evaluateFFDDS(const CoordType pt[D], T ds[DOut]) const
{
    T deriv[D][DOut];
    GADGET_CHECK_RETURN_FALSE(this->evaluateFFDDerivative(pt, deriv));
    memcpy(ds, deriv+3*sizeof(T)*DOut, sizeof(T)*DOut);
    return true;
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool FFDBase<T, CoordType, DIn, DOut>::evaluateWorldDerivative(const CoordType pt[D], T deriv[D][DOut]) const
{
    CoordType pt_g[D];
    this->world_to_grid(pt, pt_g);
    GADGET_CHECK_RETURN_FALSE(this->evaluateFFDDerivative(pt_g, deriv));

    std::vector<coord_type> spacing;
    this->get_spacing(spacing);

    unsigned int d, d2;
    for ( d=0; d<DIn; d++ )
    {
        for ( d2=0; d2<DOut; d2++ )
        {
            deriv[d][d2] /= spacing[d];
        }
    }

    return true;
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool FFDBase<T, CoordType, DIn, DOut>::evaluateWorldDX(const CoordType pt[D], T dx[DOut]) const
{
    CoordType pt_g[D];
    this->world_to_grid(pt, pt_g);
    GADGET_CHECK_RETURN_FALSE(this->evaluateFFDDX(pt_g, dx));

    coord_type sx = coord_type(1.0)/this->get_spacing(0);

    unsigned int d;
    for ( d=0; d<DOut; d++ )
    {
        dx[d] *= sx;
    }

    return true;
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool FFDBase<T, CoordType, DIn, DOut>::evaluateWorldDY(const CoordType pt[D], T dy[DOut]) const
{
    CoordType pt_g[D];
    this->world_to_grid(pt, pt_g);
    GADGET_CHECK_RETURN_FALSE(this->evaluateFFDDY(pt_g, dy));

    coord_type sy = coord_type(1.0)/this->get_spacing(1);

    unsigned int d;
    for ( d=0; d<DOut; d++ )
    {
        dy[d] *= sy;
    }

    return true;
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool FFDBase<T, CoordType, DIn, DOut>::evaluateWorldDZ(const CoordType pt[D], T dz[DOut]) const
{
    CoordType pt_g[D];
    this->world_to_grid(pt, pt_g);
    GADGET_CHECK_RETURN_FALSE(this->evaluateFFDDZ(pt_g, dz));

    coord_type sz = coord_type(1.0)/this->get_spacing(2);

    unsigned int d;
    for ( d=0; d<DOut; d++ )
    {
        dz[d] *= sz;
    }

    return true;
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool FFDBase<T, CoordType, DIn, DOut>::evaluateWorldDS(const CoordType pt[D], T ds[DOut]) const
{
    CoordType pt_g[D];
    this->world_to_grid(pt, pt_g);
    GADGET_CHECK_RETURN_FALSE(this->evaluateFFDDS(pt_g, ds));

    coord_type ss = coord_type(1.0)/this->get_spacing(3);

    unsigned int d;
    for ( d=0; d<DOut; d++ )
    {
        ds[d] *= ss;
    }

    return true;
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool FFDBase<T, CoordType, DIn, DOut>::evaluateFFDSecondOrderDerivative(const PointType& pt, T dderiv[D*D][DOut]) const
{
    GADGET_CHECK_RETURN_FALSE(this->evaluateFFDDerivative(pt.begin(), dderiv));
    return true;
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool FFDBase<T, CoordType, DIn, DOut>::evaluateFFDW(const CoordType pt[D], T r[DOut]) const
{
    CoordType pg[D];
    this->world_to_grid(pt, pg);
    return this->evaluateFFD(pg, r);
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool FFDBase<T, CoordType, DIn, DOut>::evaluateFFDW(CoordType px, CoordType py, T r[DOut]) const
{
    CoordType pg[2];
    this->world_to_grid(px, py, pg[0], pg[1]);
    return this->evaluateFFD(pg, r);
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool FFDBase<T, CoordType, DIn, DOut>::evaluateFFDW(CoordType px, CoordType py, CoordType pz, T r[DOut]) const
{
    CoordType pg[3];
    this->world_to_grid(px, py, pz, pg[0], pg[1], pg[2]);
    return this->evaluateFFD(pg, r);
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool FFDBase<T, CoordType, DIn, DOut>::evaluateFFDW(CoordType px, CoordType py, CoordType pz, CoordType ps, T r[DOut]) const
{
    CoordType pg[4];
    this->world_to_grid(px, py, pz, ps, pg[0], pg[1], pg[2], pg[3]);
    return this->evaluateFFD(pg, r);
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool FFDBase<T, CoordType, DIn, DOut>::evaluateFFDDerivativeW(const CoordType pt[D], T deriv[D][DOut]) const
{
    CoordType pg[D];
    this->world_to_grid(pt, pg);
    return this->evaluateFFDDerivative(pg, deriv);
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool FFDBase<T, CoordType, DIn, DOut>::evaluateFFDDerivativeW(CoordType px, CoordType py, T deriv[D][DOut]) const
{
    CoordType pg[2];
    this->world_to_grid(px, py, pg[0], pg[1]);
    return this->evaluateFFDDerivative(pg, deriv);
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool FFDBase<T, CoordType, DIn, DOut>::evaluateFFDDerivativeW(CoordType px, CoordType py, CoordType pz, T deriv[D][DOut]) const
{
    CoordType pg[3];
    this->world_to_grid(px, py, pz, pg[0], pg[1], pg[2]);
    return this->evaluateFFDDerivative(pg, deriv);
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool FFDBase<T, CoordType, DIn, DOut>::evaluateFFDDerivativeW(CoordType px, CoordType py, CoordType pz, CoordType ps, T deriv[D][DOut]) const
{
    CoordType pg[4];
    this->world_to_grid(px, py, pz, ps, pg[0], pg[1], pg[2], pg[3]);
    return this->evaluateFFDDerivative(pg, deriv);
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool FFDBase<T, CoordType, DIn, DOut>::evaluateFFDSecondOrderDerivativeW(const CoordType pt[D], T dderiv[D*D][DOut]) const
{
    CoordType pg[D];
    this->world_to_grid(pt, pg);
    return this->evaluateFFDDerivative(pg, dderiv);
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool FFDBase<T, CoordType, DIn, DOut>::evaluateFFDSecondOrderDerivativeW(CoordType px, CoordType py, T dderiv[D*D][DOut]) const
{
    CoordType pg[2];
    this->world_to_grid(px, py, pg[0], pg[1]);
    return this->evaluateFFDSecondOrderDerivative(pg, dderiv);
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool FFDBase<T, CoordType, DIn, DOut>::evaluateFFDSecondOrderDerivativeW(CoordType px, CoordType py, CoordType pz, T dderiv[D*D][DOut]) const
{
    CoordType pg[3];
    this->world_to_grid(px, py, pz, pg[0], pg[1], pg[2]);
    return this->evaluateFFDSecondOrderDerivative(pg, dderiv);
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool FFDBase<T, CoordType, DIn, DOut>::evaluateFFDSecondOrderDerivativeW(CoordType px, CoordType py, CoordType pz, CoordType ps, T dderiv[D*D][DOut]) const
{
    CoordType pg[4];
    this->world_to_grid(px, py, pz, ps, pg[0], pg[1], pg[2], pg[3]);
    return this->evaluateFFDSecondOrderDerivative(pg, dderiv);
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool FFDBase<T, CoordType, DIn, DOut>::world_to_grid(const CoordArrayType& pt_w, CoordArrayType& pt_g) const
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(pt_w.get_size(0)==DIn);

        if ( pt_g.dimensions_equal(&pt_w) )
        {
            pt_g = pt_w;
        }

        const CoordType* pW = pt_w.begin();
        CoordType* pG = pt_g.begin();

        size_t N = pt_w.get_size(1);

        long long n;

#pragma omp parallel for default(none) private(n) shared(N, pW, pG)
        for ( n=0; n<(long long)N; n++ )
        {
            this->world_to_grid(pW+n*DIn, pG+n*DIn);
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors happened in world_to_grid(const CoordArrayType& pt_w, CoordArrayType& pt_g) ... ");
        return false;
    }

    return true;
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool FFDBase<T, CoordType, DIn, DOut>::world_to_grid(const CoordType pt_w[D], CoordType pt_g[D]) const
{
    try
    {
        this->ctrl_pt_[0].world_to_image(pt_w, pt_g);
    }
    catch(...)
    {
        GERROR_STREAM("Errors happened in world_to_grid(const CoordType pt_w[D], CoordType pt_g[D]) ... ");
        return false;
    }

    return true;
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool FFDBase<T, CoordType, DIn, DOut>::world_to_grid(CoordType px_w, CoordType py_w, CoordType& px_g, CoordType& py_g) const
{
    GADGET_CHECK_RETURN_FALSE(DIn==2);

    try
    {
        this->ctrl_pt_[0].world_to_image(px_w, py_w, px_g, py_g);
    }
    catch(...)
    {
        GERROR_STREAM("Errors happened in world_to_grid(CoordType px_w, CoordType py_w, CoordType& px_g, CoordType& py_g) ... ");
        return false;
    }

    return true;
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool FFDBase<T, CoordType, DIn, DOut>::world_to_grid(CoordType px_w, CoordType py_w, CoordType pz_w, CoordType& px_g, CoordType& py_g, CoordType& pz_g) const
{
    GADGET_CHECK_RETURN_FALSE(DIn==3);

    try
    {
        this->ctrl_pt_[0].world_to_image(px_w, py_w, pz_w, px_g, py_g, pz_g);
    }
    catch(...)
    {
        GERROR_STREAM("Errors happened in world_to_grid(CoordType px_w, CoordType py_w, CoordType pz_w, CoordType& px_g, CoordType& py_g, CoordType& pz_g) ... ");
        return false;
    }

    return true;
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool FFDBase<T, CoordType, DIn, DOut>::world_to_grid(CoordType px_w, CoordType py_w, CoordType pz_w, CoordType ps_w, CoordType& px_g, CoordType& py_g, CoordType& pz_g, CoordType& ps_g) const
{
    GADGET_CHECK_RETURN_FALSE(DIn==4);

    try
    {
        this->ctrl_pt_[0].world_to_image(px_w, py_w, pz_w, ps_w, px_g, py_g, pz_g, ps_g);
    }
    catch(...)
    {
        GERROR_STREAM("Errors happened in world_to_grid(CoordType px_w, CoordType py_w, CoordType pz_w, CoordType ps_w, CoordType& px_g, CoordType& py_g, CoordType& pz_g, CoordType& ps_g) ... ");
        return false;
    }

    return true;
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool FFDBase<T, CoordType, DIn, DOut>::grid_to_world(const CoordArrayType& pt_g, CoordArrayType& pt_w) const
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(pt_g.get_size(0)==DIn);

        if ( pt_w.dimensions_equal(&pt_g) )
        {
            pt_w = pt_g;
        }

        const CoordType* pG = pt_g.begin();
        CoordType* pW = pt_w.begin();

        size_t N = pt_g.get_size(1);

        long long n;

#pragma omp parallel for default(none) private(n) shared(N, pG, pW)
        for ( n=0; n<(long long)N; n++ )
        {
            this->grid_to_world(pG+n*DIn, pW+n*DIn);
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors happened in grid_to_world(const CoordArrayType& pt_g, CoordArrayType& pt_w) ... ");
        return false;
    }

    return true;
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool FFDBase<T, CoordType, DIn, DOut>::grid_to_world(const CoordType pt_g[D], CoordType pt_w[D]) const
{
    try
    {
        this->ctrl_pt_[0].image_to_world(pt_g, pt_w);
    }
    catch(...)
    {
        GERROR_STREAM("Errors happened in grid_to_world(const CoordType pt_g[D], CoordType pt_w[D]) ... ");
        return false;
    }

    return true;
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool FFDBase<T, CoordType, DIn, DOut>::grid_to_world(CoordType px_g, CoordType py_g, CoordType& px_w, CoordType& py_w) const
{
    GADGET_CHECK_RETURN_FALSE(DIn==2);

    try
    {
        this->ctrl_pt_[0].image_to_world(px_g, py_g, px_w, py_w);
    }
    catch(...)
    {
        GERROR_STREAM("Errors happened in grid_to_world(CoordType px_g, CoordType py_g, CoordType& px_w, CoordType& py_w) ... ");
        return false;
    }

    return true;
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool FFDBase<T, CoordType, DIn, DOut>::grid_to_world(CoordType px_g, CoordType py_g, CoordType pz_g, CoordType& px_w, CoordType& py_w, CoordType& pz_w) const
{
    GADGET_CHECK_RETURN_FALSE(DIn==3);

    try
    {
        this->ctrl_pt_[0].image_to_world(px_g, py_g, pz_g, px_w, py_w, pz_w);
    }
    catch(...)
    {
        GERROR_STREAM("Errors happened in grid_to_world(CoordType px_g, CoordType py_g, CoordType pz_g, CoordType& px_w, CoordType& py_w, CoordType& pz_w) ... ");
        return false;
    }

    return true;
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool FFDBase<T, CoordType, DIn, DOut>::ffdApprox(const CoordArrayType& pos, ValueArrayType& value, ValueArrayType& residual, real_value_type& totalResidual, size_t N, size_t numOfRefinement)
{
    size_t num;
    return this->ffdApprox(pos, value, residual, totalResidual, N, num, FLT_EPSILON, numOfRefinement);
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool FFDBase<T, CoordType, DIn, DOut>::ffdApprox(const CoordArrayType& pos, ValueArrayType& value, ValueArrayType& residual, real_value_type& totalResidual, size_t N, size_t& numOfRefinement, real_value_type thresResidual, size_t maxNumOfRefinement)
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(pos.get_size(0)==N);
        GADGET_CHECK_RETURN_FALSE(pos.get_size(1)==DIn);

        GADGET_CHECK_RETURN_FALSE(value.get_size(0)==N);
        GADGET_CHECK_RETURN_FALSE(value.get_size(1)==DOut);

        totalResidual = 0;

        if ( !residual.dimensions_equal(&value) )
        {
            residual.create(value.dimensions());
            Gadgetron::clear(residual);
        }

        CoordArrayType posg(pos);
        CoordArrayType posw(pos);
        GADGET_CHECK_RETURN_FALSE(this->grid_to_world(posg, posw));

        size_t num;
        for ( num=0; num<maxNumOfRefinement; num++ )
        {
            GADGET_CHECK_RETURN_FALSE(this->ffdApprox(posg, value, residual, totalResidual, N));

            GDEBUG_CONDITION_STREAM(performTiming_, "BSpline FFD refinement " << num << " has residual of " << totalResidual);

            if ( totalResidual < thresResidual )
            {
                GDEBUG_STREAM("BSpline FFD residual is too small : " << totalResidual);
                GDEBUG_STREAM("No further refinement will be computed ... ");
                break;
            }

            if ( num<maxNumOfRefinement-1 )
            {
                GADGET_CHECK_RETURN_FALSE(this->refine());
                GADGET_CHECK_RETURN_FALSE(this->world_to_grid(posw, posg));
            }
        }

        numOfRefinement = num;
    }
    catch(...)
    {
        GERROR_STREAM("Error happened in ffdApprox(const CoordArrayType& pos, ValueArrayType& value, ValueArrayType& residual, real_value_type& totalResidual, size_t N, size_t& numOfRefinement, real_value_type thresResidual, size_t maxNumOfRefinement) ... ");
        return false;
    }

    return true;
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool FFDBase<T, CoordType, DIn, DOut>::ffdApproxW(const CoordArrayType& pos, ValueArrayType& value, ValueArrayType& residual, real_value_type& totalResidual, size_t N, size_t numOfRefinement)
{
    size_t num;
    return this->ffdApproxW(pos, value, residual, totalResidual, N, num, FLT_EPSILON, numOfRefinement);
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool FFDBase<T, CoordType, DIn, DOut>::ffdApproxW(const CoordArrayType& pos, ValueArrayType& value, ValueArrayType& residual, real_value_type& totalResidual, size_t N, size_t& numOfRefinement, real_value_type thresResidual, size_t maxNumOfRefinement)
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(pos.get_size(0)==DIn);
        GADGET_CHECK_RETURN_FALSE(pos.get_size(1)==N);

        GADGET_CHECK_RETURN_FALSE(value.get_size(0)==DOut);
        GADGET_CHECK_RETURN_FALSE(value.get_size(1)==N);

        totalResidual = 0;

        if ( !residual.dimensions_equal(&value) )
        {
            residual.create(value.dimensions());
            Gadgetron::clear(residual);
        }

        CoordArrayType posg(pos);
        GADGET_CHECK_RETURN_FALSE(this->world_to_grid(pos, posg));

        size_t num;
        for ( num=0; num<maxNumOfRefinement; num++ )
        {
            GADGET_CHECK_RETURN_FALSE(this->ffdApprox(posg, value, residual, totalResidual, N));

            GDEBUG_CONDITION_STREAM(performTiming_, "BSpline FFD refinement " << num << " has residual of " << totalResidual);

            if ( totalResidual < thresResidual )
            {
                GDEBUG_STREAM("BSpline FFD residual is too small : " << totalResidual);
                GDEBUG_STREAM("No further refinement will be computed ... ");
                break;
            }

            if ( num<maxNumOfRefinement-1 )
            {
                GADGET_CHECK_RETURN_FALSE(this->refine());
                GADGET_CHECK_RETURN_FALSE(this->world_to_grid(pos, posg));
            }
        }

        numOfRefinement = num;
    }
    catch(...)
    {
        GERROR_STREAM("Error happened in ffdApprox(const CoordArrayType& pos, ValueArrayType& value, ValueArrayType& residual, real_value_type& totalResidual, size_t N, size_t& numOfRefinement, real_value_type thresResidual, size_t maxNumOfRefinement) ... ");
        return false;
    }

    return true;
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool FFDBase<T, CoordType, DIn, DOut>::evaluateFFDOnImage(ImageType& target) const
{
    GADGET_CHECK_RETURN_FALSE(DOut==1);
    return this->evaluateFFDOnImage(&target);
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool FFDBase<T, CoordType, DIn, DOut>::evaluateFFDOnImage(ImageType target[DOut]) const
{
    try
    {
        if ( DIn==2 )
        {
            size_t sx = target[0].get_size(0);
            size_t sy = target[0].get_size(1);

            long long y;

#pragma omp parallel private(y) shared(sx, sy, target)
            {
                coord_type px, py, pg[2];
                T v[DOut];
                unsigned int d;

#pragma omp for 
                for ( y=0; y<(long long)sy; y++ )
                {
                    for ( size_t x=0; x<sx; x++ )
                    {
                        size_t offset = x + y*sx;

                        // target to world
                        target[0].image_to_world(x, size_t(y), px, py);

                        // world to grid
                        this->world_to_grid(px, py, pg[0], pg[1]);

                        // evaluate the FFD
                        this->evaluateFFD(pg, v);

                        if ( DOut == 1 )
                        {
                            target[0](offset) = v[0];
                        }
                        else if ( DOut == 2 )
                        {
                            target[0](offset) = v[0];
                            target[1](offset) = v[1];
                        }
                        else if ( DOut == 3 )
                        {
                            target[0](offset) = v[0];
                            target[1](offset) = v[1];
                            target[2](offset) = v[2];
                        }
                        else
                        {
                            for ( d=0; d<DOut; d++ )
                            {
                                target[d](offset) = v[d];
                            }
                        }
                    }
                }
            }
        }
        else if ( DIn==3 )
        {
            size_t sx = target[0].get_size(0);
            size_t sy = target[0].get_size(1);
            size_t sz = target[0].get_size(2);

            long long z;

#pragma omp parallel private(z) shared(sx, sy, sz, target)
            {
                coord_type px, py, pz, pg[3];
                T v[DOut];
                unsigned int d;

#pragma omp for 
                for ( z=0; z<(long long)sz; z++ )
                {
                    for ( size_t y=0; y<sy; y++ )
                    {
                        size_t offset = y*sx + z*sx*sy;

                        for ( size_t x=0; x<sx; x++ )
                        {
                            // target to world
                            target[0].image_to_world(x, y, size_t(z), px, py, pz);

                            // world to grid
                            this->world_to_grid(px, py, pz, pg[0], pg[1], pg[2]);

                            // evaluate the FFD
                            this->evaluateFFD(pg, v);

                            if ( DOut == 1 )
                            {
                                target[0](offset+x) = v[0];
                            }
                            else if ( DOut == 2 )
                            {
                                target[0](offset+x) = v[0];
                                target[1](offset+x) = v[1];
                            }
                            else if ( DOut == 3 )
                            {
                                target[0](offset+x) = v[0];
                                target[1](offset+x) = v[1];
                                target[2](offset+x) = v[2];
                            }
                            else
                            {
                                for ( d=0; d<DOut; d++ )
                                {
                                    target[d](offset+x) = v[d];
                                }
                            }
                        }
                    }
                }
            }
        }
        else
        {
            size_t numOfPixels = target[0].get_number_of_elements();

            long long n;

#pragma omp parallel private(n) shared(numOfPixels, target)
            {
                size_t ind_target[DIn];
                coord_type pt_target[DIn];
                coord_type pt_grid[DIn];
                T v[DOut];
                unsigned int d;

#pragma omp for 
                for ( n=0; n<(long long)numOfPixels; n++ )
                {
                    // target to world
                    target[0].calculate_index( size_t(n), ind_target );

                    target[0].image_to_world(ind_target, pt_target);

                    // world to grid
                    this->world_to_grid(pt_target, pt_grid);

                    // evaluate the FFD
                    this->evaluateFFD(pt_grid, v);

                    if ( DOut == 1 )
                    {
                        target[0](n) = v[0];
                    }
                    else if ( DOut == 2 )
                    {
                        target[0](n) = v[0];
                        target[1](n) = v[1];
                    }
                    else if ( DOut == 3 )
                    {
                        target[0](n) = v[0];
                        target[1](n) = v[1];
                        target[2](n) = v[2];
                    }
                    else
                    {
                        for ( d=0; d<DOut; d++ )
                        {
                            target[d](n) = v[d];
                        }
                    }
                }
            }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Error happened in evaluateFFD(ImageType target[DOut]) const ... ");
        return false;
    }

    return true;
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
bool FFDBase<T, CoordType, DIn, DOut>::evaluateFFDOnArray(ArrayType& target) const
{
    GADGET_CHECK_RETURN_FALSE(DOut==1);
    return this->evaluateFFDOnArray(&target);
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
bool FFDBase<T, CoordType, DIn, DOut>::evaluateFFDOnArray(ArrayType target[DOut]) const
{
    try
    {
        if ( DIn==2 )
        {
            size_t sx = target[0].get_size(0);
            size_t sy = target[0].get_size(1);

            long long y;

#pragma omp parallel private(y) shared(sx, sy, target)
            {
                coord_type pg[2];
                T v[DOut];
                unsigned int d;

#pragma omp for 
                for ( y=0; y<(long long)sy; y++ )
                {
                    for ( size_t x=0; x<sx; x++ )
                    {
                        size_t offset = x + y*sx;

                        this->world_to_grid((CoordType)x, (CoordType)y, pg[0], pg[1]);

                        // evaluate the FFD
                        this->evaluateFFD(pg, v);

                        if ( DOut == 1 )
                        {
                            target[0](offset) = v[0];
                        }
                        else if ( DOut == 2 )
                        {
                            target[0](offset) = v[0];
                            target[1](offset) = v[1];
                        }
                        else if ( DOut == 3 )
                        {
                            target[0](offset) = v[0];
                            target[1](offset) = v[1];
                            target[2](offset) = v[2];
                        }
                        else
                        {
                            for ( d=0; d<DOut; d++ )
                            {
                                target[d](offset) = v[d];
                            }
                        }
                    }
                }
            }
        }
        else if ( DIn==3 )
        {
            size_t sx = target[0].get_size(0);
            size_t sy = target[0].get_size(1);
            size_t sz = target[0].get_size(2);

            long long z;

#pragma omp parallel private(z) shared(sx, sy, sz, target)
            {
                coord_type pg[3];
                T v[DOut];
                unsigned int d;

#pragma omp for 
                for ( z=0; z<(long long)sz; z++ )
                {
                    for ( size_t y=0; y<sy; y++ )
                    {
                        size_t offset = y*sx + z*sx*sy;

                        for ( size_t x=0; x<sx; x++ )
                        {
                            this->world_to_grid((CoordType)x, (CoordType)y, (CoordType)z, pg[0], pg[1], pg[2]);

                            // evaluate the FFD
                            this->evaluateFFD(pg, v);

                            if ( DOut == 1 )
                            {
                                target[0](offset+x) = v[0];
                            }
                            else if ( DOut == 2 )
                            {
                                target[0](offset+x) = v[0];
                                target[1](offset+x) = v[1];
                            }
                            else if ( DOut == 3 )
                            {
                                target[0](offset+x) = v[0];
                                target[1](offset+x) = v[1];
                                target[2](offset+x) = v[2];
                            }
                            else
                            {
                                for ( d=0; d<DOut; d++ )
                                {
                                    target[d](offset+x) = v[d];
                                }
                            }
                        }
                    }
                }
            }
        }
        else
        {
            size_t numOfPixels = target[0].get_number_of_elements();

            long long n;

#pragma omp parallel private(n) shared(numOfPixels, target)
            {
                std::vector<size_t> ind_target(DIn);
                coord_type pt_target[DIn];
                coord_type pt_grid[DIn];
                T v[DOut];
                unsigned int d;

#pragma omp for 
                for ( n=0; n<(long long)numOfPixels; n++ )
                {
                    ind_target = target[0].calculate_index( size_t(n) );

                    for ( d=0; d<DIn; d++ )
                    {
                        pt_target[d] = (CoordType)ind_target[d];
                    }

                    this->world_to_grid(pt_target, pt_grid);

                    // evaluate the FFD
                    this->evaluateFFD(pt_grid, v);

                    if ( DOut == 1 )
                    {
                        target[0](n) = v[0];
                    }
                    else if ( DOut == 2 )
                    {
                        target[0](n) = v[0];
                        target[1](n) = v[1];
                    }
                    else if ( DOut == 3 )
                    {
                        target[0](n) = v[0];
                        target[1](n) = v[1];
                        target[2](n) = v[2];
                    }
                    else
                    {
                        for ( d=0; d<DOut; d++ )
                        {
                            target[d](n) = v[d];
                        }
                    }
                }
            }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Error happened in evaluateFFD(ArrayType target[DOut]) const ... ");
        return false;
    }

    return true;
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
bool FFDBase<T, CoordType, DIn, DOut>::imageToFFDInputsW(ImageType target[DOut], CoordArrayType& pos, ValueArrayType& value)
{
    try
    {
        size_t N = target[0].get_number_of_elements();
        pos.create(DIn, N);
        value.create(DOut, N);

        if ( DIn==2 )
        {
            size_t sx = target[0].get_size(0);
            size_t sy = target[0].get_size(1);

            long long y;

#pragma omp parallel private(y) shared(sx, sy, target, pos, value)
            {
                coord_type px, py;
                unsigned int d;

#pragma omp for 
                for ( y=0; y<(long long)sy; y++ )
                {
                    for ( size_t x=0; x<sx; x++ )
                    {
                        size_t offset = x + y*sx;

                        // target to world
                        target[0].image_to_world(x, size_t(y), px, py);

                        pos(0, offset) = px;
                        pos(1, offset) = py;

                        if ( DOut == 1 )
                        {
                            value(0, offset) = target[0](offset);
                        }
                        else if ( DOut == 2 )
                        {
                            value(0, offset) = target[0](offset);
                            value(1, offset) = target[1](offset);
                        }
                        else if ( DOut == 3 )
                        {
                            value(0, offset) = target[0](offset);
                            value(1, offset) = target[1](offset);
                            value(2, offset) = target[2](offset);
                        }
                        else
                        {
                            for ( d=0; d<DOut; d++ )
                            {
                                value(d, offset) = target[d](offset);
                            }
                        }
                    }
                }
            }
        }
        else if ( DIn==3 )
        {
            size_t sx = target[0].get_size(0);
            size_t sy = target[0].get_size(1);
            size_t sz = target[0].get_size(2);

            long long z;

#pragma omp parallel private(z) shared(sx, sy, sz, target, pos, value)
            {
                coord_type px, py, pz;
                unsigned int d;

#pragma omp for 
                for ( z=0; z<(long long)sz; z++ )
                {
                    for ( size_t y=0; y<sy; y++ )
                    {
                        size_t offset = y*sx + z*sx*sy;

                        for ( size_t x=0; x<sx; x++ )
                        {
                            // target to world
                            target[0].image_to_world(x, y, size_t(z), px, py, pz);

                            pos(0, offset) = px;
                            pos(1, offset) = py;
                            pos(2, offset) = pz;

                            if ( DOut == 1 )
                            {
                                value(0, offset) = target[0](offset);
                            }
                            else if ( DOut == 2 )
                            {
                                value(0, offset) = target[0](offset);
                                value(1, offset) = target[1](offset);
                            }
                            else if ( DOut == 3 )
                            {
                                value(0, offset) = target[0](offset);
                                value(1, offset) = target[1](offset);
                                value(2, offset) = target[2](offset);
                            }
                            else
                            {
                                for ( d=0; d<DOut; d++ )
                                {
                                    value(d, offset) = target[d](offset);
                                }
                            }
                        }
                    }
                }
            }
        }
        else
        {
            long long n;

#pragma omp parallel private(n) shared(N, target)
            {
                size_t ind_target[DIn];
                coord_type pt_target[DIn];
                unsigned int d;

#pragma omp for 
                for ( n=0; n<(long long)N; n++ )
                {
                    // target to world
                    target[0].calculate_index( size_t(n), ind_target );

                    target[0].image_to_world(ind_target, pt_target);

                    for ( d=0; d<DIn; d++ )
                    {
                        pos(d, n) = pt_target[d];
                    }

                    if ( DOut == 1 )
                    {
                        value(0, n) = target[0](n);
                    }
                    else if ( DOut == 2 )
                    {
                        value(0, n) = target[0](n);
                        value(1, n) = target[1](n);
                    }
                    else if ( DOut == 3 )
                    {
                        value(0, n) = target[0](n);
                        value(1, n) = target[1](n);
                        value(2, n) = target[2](n);
                    }
                    else
                    {
                        for ( d=0; d<DOut; d++ )
                        {
                            value(d, n) = target[d](n);
                        }
                    }
                }
            }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Error happened in imageToFFDInputsW(ImageType target[DOut], CoordArrayType& pos, ValueArrayType& value) ... ");
        return false;
    }

    return true;
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
bool FFDBase<T, CoordType, DIn, DOut>::imageToFFDInputsW(ImageType target[DOut], const MaskArrayType& mask, CoordArrayType& pos, ValueArrayType& value)
{
    try
    {
        size_t N = target[0].get_number_of_elements();
        if ( mask.get_number_of_elements() != N ) return true;

        size_t n, d;
        size_t numOfPixels = 0;
        for ( n=0; n<N; n++ )
        {
            if ( mask(n)!=0 ) numOfPixels++;
        }

        CoordArrayType posTmp;
        ValueArrayType valueTmp;

        GADGET_CHECK_RETURN_FALSE(this->imageToFFDInputsW(target, posTmp, valueTmp));

        pos.create(DIn, numOfPixels);
        value.create(DOut, numOfPixels);

        numOfPixels = 0;
        for ( n=0; n<N; n++ )
        {
            if ( mask(n)!=0 )
            {
                memcpy(pos.begin()+numOfPixels*DIn, posTmp.begin()+n*DIn, sizeof(T)*DIn);

                if ( DOut == 1 )
                {
                    value(0, numOfPixels) = valueTmp(0, n);
                }
                else if ( DOut == 2 )
                {
                    value(0, numOfPixels) = valueTmp(0, n);
                    value(1, numOfPixels) = valueTmp(1, n);
                }
                else if ( DOut == 3 )
                {
                    value(0, numOfPixels) = valueTmp(0, n);
                    value(1, numOfPixels) = valueTmp(1, n);
                    value(2, numOfPixels) = valueTmp(2, n);
                }
                else
                {
                    for ( d=0; d<DOut; d++ )
                    {
                        value(d, numOfPixels) = valueTmp(d, n);
                    }
                }

                numOfPixels++;
            }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Error happened in imageToFFDInputsW(ImageType target[DOut], const MaskArrayType& mask, CoordArrayType& pos, ValueArrayType& value) ... ");
        return false;
    }

    return true;
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
bool FFDBase<T, CoordType, DIn, DOut>::arrayToFFDInputsW(ArrayType target[DOut], CoordArrayType& pos, ValueArrayType& value)
{
    try
    {
        size_t N = target[0].get_number_of_elements();
        pos.create(DIn, N);
        value.create(DOut, N);

        if ( DIn==2 )
        {
            size_t sx = target[0].get_size(0);
            size_t sy = target[0].get_size(1);

            long long y;

#pragma omp parallel private(y) shared(sx, sy, target, pos, value)
            {
                unsigned int d;

#pragma omp for 
                for ( y=0; y<(long long)sy; y++ )
                {
                    for ( size_t x=0; x<sx; x++ )
                    {
                        size_t offset = x + y*sx;

                        pos(0, offset) = (CoordType)x;
                        pos(1, offset) = (CoordType)y;

                        if ( DOut == 1 )
                        {
                            value(0, offset) = target[0](offset);
                        }
                        else if ( DOut == 2 )
                        {
                            value(0, offset) = target[0](offset);
                            value(1, offset) = target[1](offset);
                        }
                        else if ( DOut == 3 )
                        {
                            value(0, offset) = target[0](offset);
                            value(1, offset) = target[1](offset);
                            value(2, offset) = target[2](offset);
                        }
                        else
                        {
                            for ( d=0; d<DOut; d++ )
                            {
                                value(d, offset) = target[d](offset);
                            }
                        }
                    }
                }
            }
        }
        else if ( DIn==3 )
        {
            size_t sx = target[0].get_size(0);
            size_t sy = target[0].get_size(1);
            size_t sz = target[0].get_size(2);

            long long z;

#pragma omp parallel private(z) shared(sx, sy, sz, target, pos, value)
            {
                unsigned int d;

#pragma omp for 
                for ( z=0; z<(long long)sz; z++ )
                {
                    for ( size_t y=0; y<sy; y++ )
                    {
                        size_t offset = y*sx + z*sx*sy;

                        for ( size_t x=0; x<sx; x++ )
                        {
                            pos(0, offset) = (CoordType)x;
                            pos(1, offset) = (CoordType)y;
                            pos(2, offset) = (CoordType)z;

                            if ( DOut == 1 )
                            {
                                value(0, offset) = target[0](offset);
                            }
                            else if ( DOut == 2 )
                            {
                                value(0, offset) = target[0](offset);
                                value(1, offset) = target[1](offset);
                            }
                            else if ( DOut == 3 )
                            {
                                value(0, offset) = target[0](offset);
                                value(1, offset) = target[1](offset);
                                value(2, offset) = target[2](offset);
                            }
                            else
                            {
                                for ( d=0; d<DOut; d++ )
                                {
                                    value(d, offset) = target[d](offset);
                                }
                            }
                        }
                    }
                }
            }
        }
        else
        {
            long long n;

#pragma omp parallel private(n) shared(N, target)
            {
                std::vector<size_t> ind_target(DIn);
                unsigned int d;

#pragma omp for 
                for ( n=0; n<(long long)N; n++ )
                {
                    ind_target = target[0].calculate_index( size_t(n) );

                    for ( d=0; d<DIn; d++ )
                    {
                        pos(d, n) = (CoordType)ind_target[d];
                    }

                    if ( DOut == 1 )
                    {
                        value(0, n) = target[0](n);
                    }
                    else if ( DOut == 2 )
                    {
                        value(0, n) = target[0](n);
                        value(1, n) = target[1](n);
                    }
                    else if ( DOut == 3 )
                    {
                        value(0, n) = target[0](n);
                        value(1, n) = target[1](n);
                        value(2, n) = target[2](n);
                    }
                    else
                    {
                        for ( d=0; d<DOut; d++ )
                        {
                            value(d, n) = target[d](n);
                        }
                    }
                }
            }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Error happened in imageToFFDInputsW(ImageType target[DOut], CoordArrayType& pos, ValueArrayType& value) ... ");
        return false;
    }

    return true;
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
bool FFDBase<T, CoordType, DIn, DOut>::arrayToFFDInputsW(ArrayType target[DOut], const MaskArrayType& mask, CoordArrayType& pos, ValueArrayType& value)
{
    try
    {
        size_t N = target[0].get_number_of_elements();
        if ( mask.get_number_of_elements() != N ) return true;

        size_t n, d;
        size_t numOfPixels = 0;
        for ( n=0; n<N; n++ )
        {
            if ( mask(n)!=0 ) numOfPixels++;
        }

        CoordArrayType posTmp;
        ValueArrayType valueTmp;

        GADGET_CHECK_RETURN_FALSE(this->arrayToFFDInputsW(target, posTmp, valueTmp));

        pos.create(DIn, numOfPixels);
        value.create(DOut, numOfPixels);

        numOfPixels = 0;
        for ( n=0; n<N; n++ )
        {
            if ( mask(n)!=0 )
            {
                memcpy(pos.begin()+numOfPixels*DIn, posTmp.begin()+n*DIn, sizeof(T)*DIn);

                if ( DOut == 1 )
                {
                    value(0, numOfPixels) = valueTmp(0, n);
                }
                else if ( DOut == 2 )
                {
                    value(0, numOfPixels) = valueTmp(0, n);
                    value(1, numOfPixels) = valueTmp(1, n);
                }
                else if ( DOut == 3 )
                {
                    value(0, numOfPixels) = valueTmp(0, n);
                    value(1, numOfPixels) = valueTmp(1, n);
                    value(2, numOfPixels) = valueTmp(2, n);
                }
                else
                {
                    for ( d=0; d<DOut; d++ )
                    {
                        value(d, numOfPixels) = valueTmp(d, n);
                    }
                }

                numOfPixels++;
            }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Error happened in imageToFFDInputsW(ImageType target[DOut], const MaskArrayType& mask, CoordArrayType& pos, ValueArrayType& value) ... ");
        return false;
    }

    return true;
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool FFDBase<T, CoordType, DIn, DOut>::ffdApproxImage(ImageType target[DOut], real_value_type& totalResidual, size_t numOfRefinement)
{
    try
    {
        size_t N = target[0].get_number_of_elements();

        CoordArrayType pos(DIn, N);
        ValueArrayType value(DOut, N);
        ValueArrayType residual(DOut, N);

        GADGET_CHECK_RETURN_FALSE(this->imageToFFDInputsW(target, pos, value));
        GADGET_CHECK_RETURN_FALSE(this->ffdApproxW(pos, value, residual, totalResidual, N, numOfRefinement));
    }
    catch(...)
    {
        GERROR_STREAM("Error happened in ffdApproxImage(ImageType target[DOut], real_value_type& totalResidual, size_t numOfRefinement) ... ");
        return false;
    }

    return true;
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool FFDBase<T, CoordType, DIn, DOut>::ffdApproxImage(ImageType target[DOut], real_value_type& totalResidual, size_t& numOfRefinement, real_value_type thresResidual, size_t maxNumOfRefinement)
{
    try
    {
        size_t N = target[0].get_number_of_elements();

        CoordArrayType pos(DIn, N);
        ValueArrayType value(DOut, N);
        ValueArrayType residual(DOut, N);

        GADGET_CHECK_RETURN_FALSE(this->imageToFFDInputsW(target, pos, value));
        GADGET_CHECK_RETURN_FALSE(this->ffdApproxW(pos, value, residual, totalResidual, N, numOfRefinement, thresResidual, maxNumOfRefinement));
    }
    catch(...)
    {
        GERROR_STREAM("Error happened in ffdApproxImage(ImageType target[DOut], real_value_type& totalResidual, size_t& numOfRefinement, real_value_type thresResidual, size_t maxNumOfRefinement) ... ");
        return false;
    }

    return true;
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool FFDBase<T, CoordType, DIn, DOut>::ffdApproxImage(ImageType& target, real_value_type& totalResidual, size_t numOfRefinement)
{
    GADGET_CHECK_RETURN_FALSE(DOut==1);
    return this->ffdApproxImage(&target, totalResidual, numOfRefinement);
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool FFDBase<T, CoordType, DIn, DOut>::ffdApproxImage(ImageType& target, real_value_type& totalResidual, size_t& numOfRefinement, real_value_type thresResidual, size_t maxNumOfRefinement)
{
    GADGET_CHECK_RETURN_FALSE(DOut==1);
    return this->ffdApproxImage(&target, totalResidual, numOfRefinement, thresResidual, maxNumOfRefinement);
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool FFDBase<T, CoordType, DIn, DOut>::ffdApproxImage(ImageType target[DOut], const MaskArrayType& mask, real_value_type& totalResidual, size_t numOfRefinement)
{
    try
    {
        CoordArrayType pos;
        ValueArrayType value;
        ValueArrayType residual;

        GADGET_CHECK_RETURN_FALSE(this->imageToFFDInputsW(target, mask, pos, value));

        size_t N = pos.get_size(1);
        GADGET_CHECK_RETURN_FALSE(this->ffdApproxW(pos, value, residual, totalResidual, N, numOfRefinement));
    }
    catch(...)
    {
        GERROR_STREAM("Error happened in ffdApproxImage(ImageType target[DOut], const MaskArrayType& mask, real_value_type& totalResidual, size_t numOfRefinement) ... ");
        return false;
    }

    return true;
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool FFDBase<T, CoordType, DIn, DOut>::ffdApproxImage(ImageType target[DOut], const MaskArrayType& mask, real_value_type& totalResidual, size_t& numOfRefinement, real_value_type thresResidual, size_t maxNumOfRefinement)
{
    try
    {
        CoordArrayType pos;
        ValueArrayType value;
        ValueArrayType residual;

        GADGET_CHECK_RETURN_FALSE(this->imageToFFDInputsW(target, mask, pos, value));

        size_t N = pos.get_size(1);
        GADGET_CHECK_RETURN_FALSE(this->ffdApproxW(pos, value, residual, totalResidual, N, numOfRefinement, thresResidual, maxNumOfRefinement));
    }
    catch(...)
    {
        GERROR_STREAM("Error happened in ffdApproxImage(ImageType target[DOut], const MaskArrayType& mask, real_value_type& totalResidual, size_t& numOfRefinement, real_value_type thresResidual, size_t maxNumOfRefinement) ... ");
        return false;
    }

    return true;
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool FFDBase<T, CoordType, DIn, DOut>::ffdApproxArray(ArrayType target[DOut], real_value_type& totalResidual, size_t numOfRefinement)
{
    try
    {
        size_t N = target[0].get_number_of_elements();

        CoordArrayType pos(DIn, N);
        ValueArrayType value(DOut, N);
        ValueArrayType residual(DOut, N);

        GADGET_CHECK_RETURN_FALSE(this->arrayToFFDInputsW(target, pos, value));
        GADGET_CHECK_RETURN_FALSE(this->ffdApproxW(pos, value, residual, totalResidual, N, numOfRefinement));
    }
    catch(...)
    {
        GERROR_STREAM("Error happened in ffdApproxArray(ArrayType target[DOut], real_value_type& totalResidual, size_t numOfRefinement) ... ");
        return false;
    }

    return true;
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool FFDBase<T, CoordType, DIn, DOut>::ffdApproxArray(ArrayType target[DOut], real_value_type& totalResidual, size_t& numOfRefinement, real_value_type thresResidual, size_t maxNumOfRefinement)
{
    try
    {
        size_t N = target[0].get_number_of_elements();

        CoordArrayType pos(DIn, N);
        ValueArrayType value(DOut, N);
        ValueArrayType residual(DOut, N);

        GADGET_CHECK_RETURN_FALSE(this->arrayToFFDInputsW(target, pos, value));
        GADGET_CHECK_RETURN_FALSE(this->ffdApproxW(pos, value, residual, totalResidual, N, numOfRefinement, thresResidual, maxNumOfRefinement));
    }
    catch(...)
    {
        GERROR_STREAM("Error happened in ffdApproxArray(ArrayType target[DOut], real_value_type& totalResidual, size_t numOfRefinement) ... ");
        return false;
    }

    return true;
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool FFDBase<T, CoordType, DIn, DOut>::ffdApproxArray(ArrayType& target, real_value_type& totalResidual, size_t numOfRefinement)
{
    GADGET_CHECK_RETURN_FALSE(DOut==1);
    return this->ffdApproxArray(&target, totalResidual, numOfRefinement);
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool FFDBase<T, CoordType, DIn, DOut>::ffdApproxArray(ArrayType& target, real_value_type& totalResidual, size_t& numOfRefinement, real_value_type thresResidual, size_t maxNumOfRefinement)
{
    GADGET_CHECK_RETURN_FALSE(DOut==1);
    return this->ffdApproxArray(&target, totalResidual, numOfRefinement, thresResidual, maxNumOfRefinement);
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool FFDBase<T, CoordType, DIn, DOut>::ffdApproxArray(ArrayType target[DOut], const MaskArrayType& mask, real_value_type& totalResidual, size_t numOfRefinement)
{
    try
    {
        CoordArrayType pos;
        ValueArrayType value;
        ValueArrayType residual;

        GADGET_CHECK_RETURN_FALSE(this->arrayToFFDInputsW(target, mask, pos, value));

        size_t N = pos.get_size(1);
        GADGET_CHECK_RETURN_FALSE(this->ffdApproxW(pos, value, residual, totalResidual, N, numOfRefinement));
    }
    catch(...)
    {
        GERROR_STREAM("Error happened in ffdApproxArray(ArrayType target[DOut], const MaskArrayType& mask, real_value_type& totalResidual, size_t numOfRefinement) ... ");
        return false;
    }

    return true;
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool FFDBase<T, CoordType, DIn, DOut>::ffdApproxArray(ArrayType target[DOut], const MaskArrayType& mask, real_value_type& totalResidual, size_t& numOfRefinement, real_value_type thresResidual, size_t maxNumOfRefinement)
{
    try
    {
        CoordArrayType pos;
        ValueArrayType value, residual;

        GADGET_CHECK_RETURN_FALSE(this->arrayToFFDInputsW(target, mask, pos, value));

        size_t N = pos.get_size(1);
        GADGET_CHECK_RETURN_FALSE(this->ffdApproxW(pos, value, residual, totalResidual, N, numOfRefinement, thresResidual, maxNumOfRefinement));
    }
    catch(...)
    {
        GERROR_STREAM("Error happened in ffdApproxArray(ArrayType target[DOut], const MaskArrayType& mask, real_value_type& totalResidual, size_t numOfRefinement) ... ");
        return false;
    }

    return true;
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
bool FFDBase<T, CoordType, DIn, DOut>::clear(T v)
{
    try
    {
        unsigned int d;

        if ( std::abs(v) == 0 )
        {
            for ( d=0; d<DOut; d++ )
            {
                Gadgetron::clear(ctrl_pt_[d]);
            }
        }
        else
        {
            for ( d=0; d<DOut; d++ )
            {
                Gadgetron::fill(ctrl_pt_[d], v);
            }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Error happened in FFDBase<T, CoordType, DIn, DOut>::clear(T v) ... ");
        return false;
    }

    return true;
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
void FFDBase<T, CoordType, DIn, DOut>::print(std::ostream& os) const
{
    using namespace std;

    os << "---------------------- Free Form Deformation -------------------------" << endl;
    os << "Define the interface for Free Form Deformation (FFD) " << endl;
    os << "----------------------------------------------------------------------" << endl;
}

}
