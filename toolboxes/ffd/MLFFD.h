/** \file       MLFFD.h
    \brief      Implement multi-level FreeFormDeformation
                For every level, the fitting residual from previous level will be approximated
                The final fitted value is the sum of all levels

    \author     Hui Xue
*/

#pragma once

#include "FFDBase.h"

namespace Gadgetron { 

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut>
class MLFFD : public FFDBase<T, CoordType, DIn, DOut>
{
public:

    typedef FFDBase<T, CoordType, DIn, DOut> BaseClass;
    typedef FFDBase<T, CoordType, DIn, DOut> Self;

    typedef typename BaseClass::bspline_float_type real_value_type;
    typedef real_value_type bspline_float_type;

    typedef typename BaseClass::coord_type coord_type;

    enum { D = DIn };

    typedef typename BaseClass::CoordArrayType      CoordArrayType;
    typedef typename BaseClass::ValueArrayType      ValueArrayType;
    typedef typename BaseClass::ArrayType           ArrayType;
    typedef typename BaseClass::FFDCtrlPtGridType   FFDCtrlPtGridType;
    typedef typename BaseClass::PointType           PointType;
    typedef typename BaseClass::ImageType           ImageType;

    typedef std::vector<BaseClass*> FFDArrayType;

    MLFFD(bool delete_data_on_destruct=false);
    MLFFD(const FFDArrayType& a, bool delete_data_on_destruct=false);
    MLFFD(const Self& a);

    virtual ~MLFFD();

    size_t get_size() const { return ml_ffd_.size(); }

    /// get the FFD array
    FFDArrayType& getFFDArray();
    const FFDArrayType& getFFDArray() const;

    /// set the delete flag
    bool delete_data_on_destruct() const { return delete_data_on_destruct_; }
    void delete_data_on_destruct(bool flag) { delete_data_on_destruct_ = flag; }

    /// evaluate the FFD at a grid location
    /// the input points are in the FFD grid
    virtual bool evaluateFFD(const CoordType pt[D], T r[DOut]) const;

    /// evaluate the 1st order derivative of FFD at a grid location
    /// deriv: derivative for all D dimensions and all DOut values
    virtual bool evaluateFFDDerivative(const CoordType pt[D], T deriv[D][DOut]) const;

    virtual bool evaluateFFDDX(const CoordType pt[D], T dx[DOut]) const;
    virtual bool evaluateFFDDY(const CoordType pt[D], T dy[DOut]) const;
    virtual bool evaluateFFDDZ(const CoordType pt[D], T dz[DOut]) const;
    virtual bool evaluateFFDDS(const CoordType pt[D], T ds[DOut]) const;

    /// evaluate the 2nd order derivative of FFD at a grid location
    /// dderiv : D*D vector, stores dxx dxy dxz ...; dyx dyy dyz ...; dzx dzy dzz ...
    virtual bool evaluateFFDSecondOrderDerivative(const CoordType pt[D], T dderiv[D*D][DOut]) const;

    /// compute the FFD approximation once
    /// pos : the position of input points, DIn by N
    /// value : the value on input points, DOut by N
    /// residual : the approximation residual after computing FFD, DOut by N
    /// N : the number of points
    virtual bool ffdApprox(const CoordArrayType& pos, ValueArrayType& value, ValueArrayType& residual, T& totalResidual, size_t N);
    virtual bool ffdApprox(const CoordArrayType& pos, ValueArrayType& value, ValueArrayType& residual, real_value_type& totalResidual, size_t N, size_t& numOfRefinement, real_value_type thresResidual, size_t maxNumOfRefinement);

    /// refine the FFD
    virtual bool refine();

    /// general print function
    virtual void print(std::ostream& os) const;

protected: 

    FFDArrayType ml_ffd_;

    /// if true, all stored ffd will be deleted
    bool delete_data_on_destruct_;

};

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
MLFFD<T, CoordType, DIn, DOut>::MLFFD(bool delete_data_on_destruct) : delete_data_on_destruct_(delete_data_on_destruct)
{
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
MLFFD<T, CoordType, DIn, DOut>::MLFFD(const FFDArrayType& a, bool delete_data_on_destruct) : delete_data_on_destruct_(delete_data_on_destruct)
{
    ml_ffd_.resize(a.size());
    for ( size_t ii=0; ii<a.size(); ii++ )
    {
        if ( a[ii] != NULL )
        {
            ml_ffd_[ii] = a[ii];
        }
    }
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
MLFFD<T, CoordType, DIn, DOut>::
MLFFD(const Self& a)
{
    delete_data_on_destruct_ = false;
    ml_ffd_ = a.getFFDArray();
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
MLFFD<T, CoordType, DIn, DOut>::~MLFFD()
{
    if ( delete_data_on_destruct_ )
    {
        for ( size_t ii=0; ii<ml_ffd_.size(); ii++ )
        {
            if ( ml_ffd_[ii] != NULL )
            {
                delete ml_ffd_[ii];
                ml_ffd_[ii] = NULL;
            }
        }
    }
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool MLFFD<T, CoordType, DIn, DOut>::evaluateFFD(const CoordType pt[D], T r[DOut]) const
{
    unsigned int d;
    for ( d=0; d<DOut; d++ )
    {
        r[d] = 0;
    }

    T rLevel[DOut];

    size_t ii;
    for (ii=0; ii<ml_ffd_.size(); ii++)
    {
        if ( ml_ffd_[ii] == NULL ) continue;

        GADGET_CHECK_RETURN_FALSE(ml_ffd_[ii]->evaluateFFD(pt, rLevel));
        for ( d=0; d<DOut; d++ )
        {
            r[d] += rLevel[d];
        }
    }

    return true;
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool MLFFD<T, CoordType, DIn, DOut>::evaluateFFDDerivative(const CoordType pt[D], T deriv[D][DOut]) const
{
    unsigned int d, d2;
    for ( d=0; d<D; d++ )
    {
        for ( d2=0; d2<DOut; d2++ )
        {
            deriv[d][d2] = 0;
        }
    }

    T derivLevel[D][DOut];

    size_t ii;
    for (ii=0; ii<ml_ffd_.size(); ii++)
    {
        if ( ml_ffd_[ii] == NULL ) continue;

        GADGET_CHECK_RETURN_FALSE(ml_ffd_[ii]->evaluateFFDDerivative(pt, derivLevel));

        for ( d=0; d<D; d++ )
        {
            for ( d2=0; d2<DOut; d2++ )
            {
                deriv[d][d2] += derivLevel[d][d2];
            }
        }
    }

    return true;
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool MLFFD<T, CoordType, DIn, DOut>::evaluateFFDDX(const CoordType pt[D], T dx[DOut]) const
{
    unsigned int d;
    for ( d=0; d<DOut; d++ )
    {
        dx[d] = 0;
    }

    T dxLevel[D];

    size_t ii;
    for (ii=0; ii<ml_ffd_.size(); ii++)
    {
        if ( ml_ffd_[ii] == NULL ) continue;

        GADGET_CHECK_RETURN_FALSE(ml_ffd_[ii]->evaluateFFDDX(pt, dxLevel));

        for ( d=0; d<DOut; d++ )
        {
            dx[d] += dxLevel[d];
        }
    }

    return true;
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool MLFFD<T, CoordType, DIn, DOut>::evaluateFFDDY(const CoordType pt[D], T dy[DOut]) const
{
    unsigned int d;
    for ( d=0; d<DOut; d++ )
    {
        dy[d] = 0;
    }

    T dyLevel[D];

    size_t ii;
    for (ii=0; ii<ml_ffd_.size(); ii++)
    {
        if ( ml_ffd_[ii] == NULL ) continue;

        GADGET_CHECK_RETURN_FALSE(ml_ffd_[ii]->evaluateFFDDY(pt, dyLevel));

        for ( d=0; d<DOut; d++ )
        {
            dy[d] += dyLevel[d];
        }
    }

    return true;
}


template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool MLFFD<T, CoordType, DIn, DOut>::evaluateFFDDZ(const CoordType pt[D], T dz[DOut]) const
{
    unsigned int d;
    for ( d=0; d<DOut; d++ )
    {
        dz[d] = 0;
    }

    T dzLevel[D];

    size_t ii;
    for (ii=0; ii<ml_ffd_.size(); ii++)
    {
        if ( ml_ffd_[ii] == NULL ) continue;

        GADGET_CHECK_RETURN_FALSE(ml_ffd_[ii]->evaluateFFDDZ(pt, dzLevel));

        for ( d=0; d<DOut; d++ )
        {
            dz[d] += dzlevel[d];
        }
    }

    return true;
}


template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool MLFFD<T, CoordType, DIn, DOut>::evaluateFFDDS(const CoordType pt[D], T ds[DOut]) const
{
    unsigned int d;
    for ( d=0; d<DOut; d++ )
    {
        ds[d] = 0;
    }

    T dsLevel[D];

    size_t ii;
    for (ii=0; ii<ml_ffd_.size(); ii++)
    {
        if ( ml_ffd_[ii] == NULL ) continue;

        GADGET_CHECK_RETURN_FALSE(ml_ffd_[ii]->evaluateFFDDS(pt, dsLevel));

        for ( d=0; d<DOut; d++ )
        {
            ds[d] += dslevel[d];
        }
    }

    return true;
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool MLFFD<T, CoordType, DIn, DOut>::evaluateFFDSecondOrderDerivative(const CoordType pt[D], T dderiv[D*D][DOut]) const
{
    unsigned int d, d2;
    for ( d=0; d<D*D; d++ )
    {
        for ( d2=0; d2<DOut; d2++ )
        {
            dderiv[d][d2] = 0;
        }
    }

    T dderivLevel[D*D][DOut];

    size_t ii;
    for (ii=0; ii<ml_ffd_.size(); ii++)
    {
        if ( ml_ffd_[ii] == NULL ) continue;

        GADGET_CHECK_RETURN_FALSE(ml_ffd_[ii]->evaluateFFDSecondOrderDerivative(pt, dderivLevel));

        for ( d=0; d<D*D; d++ )
        {
            for ( d2=0; d2<DOut; d2++ )
            {
                dderiv[d][d2] += dderivLevel[d][d2];
            }
        }
    }

    return true;
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool MLFFD<T, CoordType, DIn, DOut>::ffdApprox(const CoordArrayType& pos, ValueArrayType& value, ValueArrayType& residual, T& totalResidual, size_t N)
{
    ValueArrayType valueLevel(value);

    size_t ii;
    for (ii=0; ii<ml_ffd_.size(); ii++)
    {
        if ( ml_ffd_[ii] == NULL ) continue;

        GADGET_CHECK_RETURN_FALSE(ml_ffd_[ii]->ffdApprox(pos, valueLevel, residual, totalResidual, N));
        valueLevel = residual;
    }

    return true;
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool MLFFD<T, CoordType, DIn, DOut>::ffdApprox(const CoordArrayType& pos, ValueArrayType& value, ValueArrayType& residual, real_value_type& totalResidual, size_t N, size_t& numOfRefinement, real_value_type thresResidual, size_t maxNumOfRefinement)
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(pos.get_size(0)==N);
        GADGET_CHECK_RETURN_FALSE(pos.get_size(1)==DIn);

        GADGET_CHECK_RETURN_FALSE(value.get_size(0)==N);
        GADGET_CHECK_RETURN_FALSE(value.get_size(1)==DOut);

        totalResidual = 0;

        if ( !residual.dimensions_equal(value) )
        {
            residual.create(value.get_dimensions());
            Gadgetron::clear(residual);
        }

        ValueArrayType valueLevel(value);
        size_t numOfRefinementLevel(0);

        size_t ii;
        for (ii=0; ii<ml_ffd_.size(); ii++)
        {
            if ( ml_ffd_[ii] == NULL ) continue;

            GADGET_CHECK_RETURN_FALSE(ml_ffd_[ii]->ffdApprox(pos, valueLevel, residual, totalResidual, N, numOfRefinementLevel, thresResidual, maxNumOfRefinement));
            numOfRefinement += numOfRefinementLevel;
            valueLevel = residual;
        }
    }
    catch(...)
    {
        GERROR_STREAM("Error happened in ffdApprox(const CoordArrayType& pos, ValueArrayType& value, ValueArrayType& residual, real_value_type& totalResidual, size_t N, size_t& numOfRefinement, real_value_type thresResidual, size_t maxNumOfRefinement) ... ");
        return false;
    }

    return true;
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
inline bool MLFFD<T, CoordType, DIn, DOut>::refine()
{
    size_t ii;
    for (ii=0; ii<ml_ffd_.size(); ii++)
    {
        if ( ml_ffd_[ii] == NULL ) continue;

        GADGET_CHECK_RETURN_FALSE(ml_ffd_[ii]->refine());
    }

    return true;
}

template <typename T, typename CoordType, unsigned int DIn, unsigned int DOut> 
void MLFFD<T, CoordType, DIn, DOut>::print(std::ostream& os) const
{
    using namespace std;

    os << "---------------------- Multi-level Free Form Deformation ------------------" << endl;
    os << "Number of level is : " << ml_ffd_.size() << endl;

    size_t ii;
    for (ii=0; ii<ml_ffd_.size(); ii++)
    {
        os << "Level " << ii << " : " << endl;
        if ( ml_ffd_[i]!=NULL )
        {
            ml_ffd_[i]->print(os);
        }
        else
        {
            os << "--> Pointer is NULL ... " << endl;
        }
    }
    os << "------------------------------------------------------------------------------" << endl;
}

}
