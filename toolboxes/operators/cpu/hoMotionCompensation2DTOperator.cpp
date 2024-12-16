
#include "hoMotionCompensation2DTOperator.h"

namespace Gadgetron {

template <typename T, typename CoordType>
hoMotionCompensation2DTOperator<T, CoordType>::hoMotionCompensation2DTOperator() : bg_value_(0), bc_(GT_BOUNDARY_CONDITION_BORDERVALUE)
{
}

template <typename T, typename CoordType>
hoMotionCompensation2DTOperator<T, CoordType>::~hoMotionCompensation2DTOperator()
{
}

template <typename T, typename CoordType>
void hoMotionCompensation2DTOperator<T, CoordType>::mult_M(ARRAY_TYPE* x, ARRAY_TYPE* y, bool accumulate)
{
    try
    {
        if(accumulate)
        {
            moco_im_ = *y;
        }

        this->warp_image(*x, dx_, dy_,*y, bg_value_, bc_);

        if (accumulate)
        {
            Gadgetron::add(moco_im_, *y, *y);
        }
    }
    catch(...)
    {
        GADGET_THROW("Errors in hoMotionCompensation2DTOperator<T>::mult_M(...) ... ");
    }
}

template <typename T, typename CoordType>
void hoMotionCompensation2DTOperator<T, CoordType>::mult_MH(ARRAY_TYPE* x, ARRAY_TYPE* y, bool accumulate)
{
    try
    {
        if (accumulate)
        {
            adj_moco_im_ = *y;
        }

        this->warp_image(*x, adj_dx_, adj_dy_, *y, bg_value_, bc_);

        if (accumulate)
        {
            Gadgetron::add(adj_moco_im_, *y, *y);
        }
    }
    catch (...)
    {
        GADGET_THROW("Errors in hoMotionCompensation2DTOperator<T>::mult_M(...) ... ");
    }
}

template <typename T, typename CoordType>
void hoMotionCompensation2DTOperator<T, CoordType>::gradient(ARRAY_TYPE* x, ARRAY_TYPE* g, bool accumulate)
{
    try
    {
        if (accumulate)
        {
            grad_im_ = *g;
        }

        // Mx
        this->warp_image(*x, dx_, dy_, moco_im_, bg_value_, bc_);

        // 1 ./ ( (Mx)'(Mx)+mu )^1/2
        Gadgetron::abs(moco_im_, moco_im_Norm_);

        if (!moco_im_Norm_approx_.dimensions_equal(moco_im_Norm_))
        {
            moco_im_Norm_approx_.create(moco_im_Norm_.dimensions());
        }

        const value_type* pCoeffNorm = moco_im_Norm_.begin();
        value_type* pBuf = moco_im_Norm_approx_.begin();

        long long N = (long long)moco_im_Norm_.get_number_of_elements();

        value_type mu = 10 * FLT_EPSILON;

        long long ii;
#pragma omp parallel for default(none) private(ii) shared(N, pBuf, pCoeffNorm, mu)
        for (ii = 0; ii<N; ii++)
        {
            pBuf[ii] = (value_type)(1.0 / std::sqrt(pCoeffNorm[ii] * pCoeffNorm[ii] + mu));
        }

        Gadgetron::multiply(moco_im_, moco_im_Norm_approx_, moco_im_);

        this->warp_image(moco_im_, adj_dx_, adj_dy_, *g, bg_value_, bc_);

        if (accumulate)
        {
            Gadgetron::add(grad_im_, *g, *g);
        }
    }
    catch (...)
    {
        GADGET_THROW("Errors in hoMotionCompensation2DTOperator<T>::gradient(...) ... ");
    }
}

template <typename T, typename CoordType>
typename hoMotionCompensation2DTOperator<T, CoordType>::value_type hoMotionCompensation2DTOperator<T, CoordType>::magnitude(ARRAY_TYPE* x)
{
    try
    {
        this->warp_image(*x, dx_, dy_, moco_im_, bg_value_, bc_);

        Gadgetron::abs(moco_im_, moco_im_Norm_);

        value_type L1CoeffNorm = Gadgetron::asum(&moco_im_Norm_);

        return L1CoeffNorm;
    }
    catch (...)
    {
        GADGET_THROW("Errors in hoMotionCompensation2DTOperator<T>::magnitude(...) ... ");
    }
}

template <typename T, typename CoordType>
void hoMotionCompensation2DTOperator<T, CoordType>::warp_image(const hoNDArray<T>& im, const hoNDArray<CoordType>& dx, const hoNDArray<CoordType>& dy, hoNDArray<T>& warpped, T bgValue, Gadgetron::GT_BOUNDARY_CONDITION bh)
{
    try
    {
        size_t RO = im.get_size(0);
        size_t E1 = im.get_size(1);
        size_t CHA = im.get_size(2);
        size_t N = im.get_size(3);

        GADGET_CHECK_THROW(dx.get_size(0)==RO);
        GADGET_CHECK_THROW(dx.get_size(1)==E1);
        GADGET_CHECK_THROW(dx.get_size(2)==N);

        GADGET_CHECK_THROW(dy.get_size(0)==RO);
        GADGET_CHECK_THROW(dy.get_size(1)==E1);
        GADGET_CHECK_THROW(dy.get_size(2)==N);

        warpped.create(RO, E1, CHA, N);

        typedef T ValueType;
        typedef hoNDImage<ValueType, 2> ImageType;
        typedef hoNDImage<value_type, 2> RealImageType;
        typedef hoNDImage<CoordType, 2> DeformType;

        ValueType* pIm = const_cast<ValueType*>(im.begin());
        CoordType* pDx = const_cast<CoordType*>(dx.begin());
        CoordType* pDy = const_cast<CoordType*>(dy.begin());
        ValueType* pWarpped = warpped.begin();

        std::vector<size_t> dim2D(2);
        dim2D[0] = RO;
        dim2D[1] = E1;

        size_t N2D = sizeof(ValueType)*RO*E1;

        long long ii;

        long long totalN = N*CHA;
        #pragma omp parallel private(ii) shared(RO, E1, CHA, totalN, pIm, pDx, pDy, pWarpped, bgValue, bh, dim2D, N2D)
        {
            hoImageRegDeformationField<CoordType, 2> deformTransform;
            hoNDBoundaryHandlerFixedValue< ImageType > bhFixedValue;
            hoNDBoundaryHandlerBorderValue< ImageType > bhBorderValue;
            hoNDBoundaryHandlerPeriodic< ImageType > bhPeriodic;
            hoNDBoundaryHandlerMirror< ImageType > bhMirror;

            hoNDInterpolatorBSpline<ImageType, 2> interpBSpline(5);

            hoImageRegWarper<ImageType, ImageType, CoordType> warper;
            warper.setBackgroundValue(bgValue);
            warper.setTransformation(deformTransform);
            warper.setInterpolator(interpBSpline);

            ImageType targetIm, sourceIm, warppedIm;
            DeformType dxIm, dyIm;

            #pragma omp for
            for ( ii=0; ii<totalN; ii++ )
            {
                long long n = ii/CHA;
                long long cha = ii - n*CHA;

                targetIm.create(dim2D, pIm+cha*RO*E1+n*RO*E1*CHA);
                sourceIm.create(dim2D, pIm+cha*RO*E1+n*RO*E1*CHA);

                warppedIm.create(dim2D, pWarpped+cha*RO*E1+n*RO*E1*CHA);

                dxIm.create(dim2D, pDx+n*RO*E1);
                dyIm.create(dim2D, pDy+n*RO*E1);

                bhFixedValue.setArray( sourceIm );
                interpBSpline.setArray( sourceIm );

                if ( bh == GT_BOUNDARY_CONDITION_FIXEDVALUE )
                    interpBSpline.setBoundaryHandler(bhFixedValue);
                else if ( bh == GT_BOUNDARY_CONDITION_BORDERVALUE )
                    interpBSpline.setBoundaryHandler(bhBorderValue);
                else if ( bh == GT_BOUNDARY_CONDITION_PERIODIC )
                    interpBSpline.setBoundaryHandler(bhPeriodic);
                else if ( bh == GT_BOUNDARY_CONDITION_MIRROR )
                    interpBSpline.setBoundaryHandler(bhMirror);
                else
                    interpBSpline.setBoundaryHandler(bhFixedValue);

                deformTransform.setDeformationField( dxIm, 0 );
                deformTransform.setDeformationField( dyIm, 1 );

                warper.warp(targetIm, sourceIm, false, warppedIm);
            }
        }
    }
    catch(...)
    {
        GADGET_THROW("Errors in hoMotionCompensation2DTOperator<T>::warp_image(...) ... ");
    }
}

// ------------------------------------------------------------
// Instantiation
// ------------------------------------------------------------

template class hoMotionCompensation2DTOperator< std::complex<float>, double >;
template class hoMotionCompensation2DTOperator< std::complex<double>, double >;

}
