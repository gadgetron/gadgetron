/** \file   hoImageRegDissimilarityMutualInformation.h
    \brief  Define the class to compute mutual information.

            The analytical derivatives are computed by using the formula proposed at:

            [1] Gerardo Hermosillo, Christophe Chefd'Hotel, Olivier Faugeras. Variational Methods for Multimodal Image Matching. 
            International Journal of Computer Vision. December 2002, Volume 50, Issue 3, pp 329-343.
            http://link.springer.com/article/10.1023%2FA%3A1020830525823

            [2] Gerardo Hermosillo. Variational Methods for Multimodal Image Matching. PhD Thesis, UNIVERSIT�E DE NICE - SOPHIA ANTIPOLIS. May 2002.
            http://webdocs.cs.ualberta.ca/~dana/readingMedIm/papers/hermosilloPhD.pdf

            This derivative computation code is based on the listed source code at page 172 - 174 in ref [2].

    \author Hui Xue
*/

#ifndef hoImageRegDissimilarityMutualInformation_H_
#define hoImageRegDissimilarityMutualInformation_H_

#pragma once

#include "hoImageRegDissimilarityHistogramBased.h"

namespace Gadgetron {

    template<typename ImageType> 
    class hoImageRegDissimilarityMutualInformation : public hoImageRegDissimilarityHistogramBased<ImageType>
    {
    public:

        typedef hoImageRegDissimilarityMutualInformation<ImageType> Self;
        typedef hoImageRegDissimilarityHistogramBased<ImageType> BaseClass;

        enum { D = ImageType::NDIM };

        typedef typename BaseClass::InterpolatorType InterpolatorType;
        typedef typename BaseClass::ValueType ValueType;
        typedef ValueType T;
        typedef ValueType element_type;
        typedef ValueType value_type;

        typedef typename BaseClass::coord_type coord_type;

        typedef typename BaseClass::hist_value_type hist_value_type;

        hoImageRegDissimilarityMutualInformation(ValueType betaArg=ValueType(2.0), unsigned int num_bin_target=64, unsigned int num_bin_warpped=64, ValueType bg_value=ValueType(0));
        virtual ~hoImageRegDissimilarityMutualInformation();

        virtual ValueType evaluate(ImageType& w);
        virtual bool evaluateDeriv(ImageType& w);

        virtual void print(std::ostream& os) const;

        /// kernel size for density estimation
        ValueType betaArg_[2];

        using BaseClass::num_bin_target_;
        using BaseClass::num_bin_warpped_;
        using BaseClass::pv_interpolation_;
        using BaseClass::step_size_ignore_pixel_;
        using BaseClass::gt_timer1_;
        using BaseClass::gt_timer2_;
        using BaseClass::gt_timer3_;
        using BaseClass::performTiming_;
        using BaseClass::gt_exporter_;
        using BaseClass::debugFolder_;

    protected:

        using BaseClass::target_;
        using BaseClass::warpped_;
        using BaseClass::deriv_;
        using BaseClass::bg_value_;
        using BaseClass::dissimilarity_;
        using BaseClass::target;
        using BaseClass::warped;
        using BaseClass::deriv;
        using BaseClass::image_dim_;
        using BaseClass::hist_;
        using BaseClass::min_target_;
        using BaseClass::max_target_;
        using BaseClass::min_warpped_;
        using BaseClass::max_warpped_;
        using BaseClass::num_samples_in_hist_;

        hoNDArray<hist_value_type> hist_target_;
        hoNDArray<hist_value_type> hist_warpped_;

        /// these variable names are kept same as the ref [2].
        ho2DArray<hist_value_type> Hy;
        hoNDArray<hist_value_type> hy;

        ho2DArray<hist_value_type> P;
        hoNDArray<hist_value_type> p;

        ho2DArray<hist_value_type> Dist;
    };

    template<typename ImageType> 
    hoImageRegDissimilarityMutualInformation<ImageType>::
    hoImageRegDissimilarityMutualInformation(ValueType betaArg, unsigned int num_bin_target, unsigned int num_bin_warpped, ValueType bg_value) 
        : BaseClass(num_bin_target, num_bin_warpped, bg_value) 
    {
        betaArg_[0] = betaArg;
        betaArg_[1] = betaArg;
    }

    template<typename ImageType> 
    hoImageRegDissimilarityMutualInformation<ImageType>::~hoImageRegDissimilarityMutualInformation()
    {
    }

    template<typename ImageType> 
    typename hoImageRegDissimilarityMutualInformation<ImageType>::ValueType hoImageRegDissimilarityMutualInformation<ImageType>::evaluate(ImageType& w)
    {
        try
        {
            BaseClass::evaluate(w);

            /// compute entorpy
            hist_target_.create(num_bin_target_);
            Gadgetron::clear(hist_target_);

            hist_warpped_.create(num_bin_warpped_);
            Gadgetron::clear(hist_warpped_);

            if ( betaArg_[0] > 0 )
            {
                Gadgetron::filterGaussian(hist_, betaArg_);
            }

            hist_value_type histSum=0;
            histSum = Gadgetron::asum(hist_);
            Gadgetron::scal( hist_value_type(1.0/histSum), hist_);

            hist_.sumOverRow(hist_target_);
            hist_.sumOverCol(hist_warpped_);

            dissimilarity_ = 0;

            size_t t, w;

            hist_value_type log2R = hist_value_type(1.0)/log( hist_value_type(2.0) );

            for ( w=0; w<num_bin_warpped_; w++ )
            {
                hist_value_type prob_w = hist_warpped_(w);

                for ( t=0; t<num_bin_target_; t++ )
                {
                    hist_value_type prob = hist_(t, w);
                    if ( prob > 0 )
                    {
                        hist_value_type prob_t = hist_target_(t);

                        dissimilarity_ -= prob * log( prob / (prob_t * prob_w) ) * log2R;
                    }
                }
            }
        }
        catch(...)
        {
            GERROR_STREAM("Errors happened in hoImageRegDissimilarityMutualInformation<ImageType>::evaluate(ImageType& t, ImageType& w) ... ");
        }

        return this->dissimilarity_;
    }

    template<typename ImageType> 
    bool hoImageRegDissimilarityMutualInformation<ImageType>::evaluateDeriv(ImageType& w)
    {
        try
        {
            this->evaluate(w);

            Hy.createArray(num_bin_target_, num_bin_warpped_);
            hy.create(num_bin_warpped_);

            P.createArray(num_bin_target_, num_bin_warpped_);
            p.create(num_bin_warpped_);

            Dist.createArray(num_bin_target_, num_bin_warpped_);

            // hoNDBoundaryHandlerFixedValue< hoMatrix<hist_value_type> > bh_hist(hist_, 0);

            size_t t, w;

            for ( t=0; t<num_bin_target_; t++ )
            {
                Hy(t, 0) = ( hist_(t, 1) - 0 );
                Hy(t, num_bin_warpped_-1) = ( 0 - hist_(t, num_bin_warpped_-2) );

                for ( w=1; w<num_bin_warpped_-1; w++ )
                {
                    Hy(t, w) = ( hist_(t, w+1) - hist_(t, w-1) );
                }
            }

            Gadgetron::scal( (hist_value_type)(0.5), Hy);

            hoNDBoundaryHandlerFixedValue< hoNDArray<hist_value_type> > bh_hist_warpped(hist_warpped_, 0);
            for ( w=0; w<num_bin_warpped_; w++ )
            {
                hy(w) = (hist_value_type)(0.5) * ( bh_hist_warpped(w+1) - bh_hist_warpped(w-1) );
            }

            P = Hy;
            p = hy;

            for ( t=0; t<num_bin_target_; t++ )
            {
                for ( w=0; w<num_bin_warpped_; w++ )
                {
                    hist_value_type v = hist_(t, w);

                    if ( v > 0 )
                    {
                        P(t, w) = Hy(t, w)/v;
                    }
                }
            }

            for ( w=0; w<num_bin_warpped_; w++ )
            {
                hist_value_type v = hist_warpped_(w);

                if ( v > 0 )
                {
                    p(w) = hy(w)/v;
                }
            }

            for ( t=0; t<num_bin_target_; t++ )
            {
                for ( w=0; w<num_bin_warpped_; w++ )
                {
                    Dist(t, w) = P(t, w) - p(w);
                }
            }

            if ( betaArg_[0] > 0 )
            {
                Gadgetron::filterGaussian(Dist, betaArg_);
            }

            hoNDBoundaryHandlerFixedValue< ho2DArray<hist_value_type> > bh_Dist(Dist, 0);
            hoNDInterpolatorLinear< ho2DArray<hist_value_type> > interp_Dist(Dist, bh_Dist);

            size_t N = target_->get_number_of_elements();

            ValueType range_t = ValueType(1.0)/(max_target_ - min_target_ + std::numeric_limits<ValueType>::epsilon());
            ValueType range_w = ValueType(1.0)/(max_warpped_ - min_warpped_ + std::numeric_limits<ValueType>::epsilon());

            long long n;

            ValueType v = (ValueType)(1.0/N);
            for ( n=0; n<(long long)N; n++ )
            {
                coord_type it = (coord_type)(range_t*(target(n)-min_target_)*(num_bin_target_-1));
                coord_type iw = (coord_type)(range_w*(warped(n)-min_warpped_)*(num_bin_warpped_-1));

                deriv_(n) = ValueType( interp_Dist(it, iw) ) * v;
            }

            // Gadgetron::math::scal(deriv_.get_number_of_elements(), ValueType(1.0/N), deriv_.begin());
        }
        catch(...)
        {
            GERROR_STREAM("Errors happened in hoImageRegDissimilarityMutualInformation<ImageType>::evaluate() ... ");
            return false;
        }

        return true;
    }

    template<typename ImageType> 
    void hoImageRegDissimilarityMutualInformation<ImageType>::print(std::ostream& os) const
    {
        using namespace std;
        os << "--------------Gagdgetron mutual information image dissimilarity meausre -------------" << endl;
        os << "Image dimension is : " << D << endl;

        std::string elemTypeName = std::string(typeid(ValueType).name());
        os << "Transformation data type is : " << elemTypeName << std::endl;

        os << "Number of intensity bins for target is : " << num_bin_target_ << endl;
        os << "Number of intensity bins for warped is : " << num_bin_warpped_ << endl;
        os << "PV interpolation for histogram is : " << pv_interpolation_ << endl;
        os << "Step size to ignore pixels when creating histogram is : " << step_size_ignore_pixel_ << endl;
        os << "Kernel size for probability density estimation is : " << betaArg_[0] << " x " << betaArg_[1] << endl;
    }
}
#endif // hoImageRegDissimilarityMutualInformation_H_
