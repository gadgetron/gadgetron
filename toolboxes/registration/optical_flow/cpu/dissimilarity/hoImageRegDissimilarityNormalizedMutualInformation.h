/** \file   hoImageRegDissimilarityNormalizedMutualInformation.h
    \brief  Define the class to compute normalized mutual information.

            C. Studholme, D.L.G. Hill, D.J. Hawkes. An overlap invariant entropy measure of 3D medical image alignment. Pattern Recognition, 32, 71-86, 1999.
            http://eecs.vanderbilt.edu/courses/cs359/other_links/papers/studholme_NMI_1999.pdf

    \author Hui Xue
*/

#ifndef hoImageRegDissimilarityNormalizedMutualInformation_H_
#define hoImageRegDissimilarityNormalizedMutualInformation_H_

#pragma once

#include "hoImageRegDissimilarityHistogramBased.h"

namespace Gadgetron {

    template<typename ImageType> 
    class hoImageRegDissimilarityNormalizedMutualInformation : public hoImageRegDissimilarityHistogramBased<ImageType>
    {
    public:

        typedef hoImageRegDissimilarityNormalizedMutualInformation<ImageType> Self;
        typedef hoImageRegDissimilarityHistogramBased<ImageType> BaseClass;

        enum { D = ImageType::NDIM };

        typedef typename BaseClass::InterpolatorType InterpolatorType;
        typedef typename BaseClass::ValueType ValueType;
        typedef ValueType T;
        typedef ValueType element_type;
        typedef ValueType value_type;

        typedef typename BaseClass::coord_type coord_type;

        typedef typename BaseClass::hist_value_type hist_value_type;

        hoImageRegDissimilarityNormalizedMutualInformation(unsigned int num_bin_target=64, unsigned int num_bin_warpped=64, ValueType bg_value=ValueType(0));
        virtual ~hoImageRegDissimilarityNormalizedMutualInformation();

        virtual ValueType evaluate(ImageType& w);

        virtual bool evaluateDeriv(ImageType& /*w*/) { return true; }

        virtual void print(std::ostream& os) const;

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
        using BaseClass::image_dim_;
        using BaseClass::bg_value_;
        using BaseClass::dissimilarity_;
        using BaseClass::target;
        using BaseClass::warped;
        using BaseClass::deriv;
        using BaseClass::hist_;
        using BaseClass::min_target_;
        using BaseClass::max_target_;
        using BaseClass::min_warpped_;
        using BaseClass::max_warpped_;
        using BaseClass::num_samples_in_hist_;

        hoNDArray<hist_value_type> hist_target_;
        hoNDArray<hist_value_type> hist_warpped_;
    };

    template<typename ImageType> 
    hoImageRegDissimilarityNormalizedMutualInformation<ImageType>::
    hoImageRegDissimilarityNormalizedMutualInformation(unsigned int num_bin_target, unsigned int num_bin_warpped, ValueType bg_value) 
        : BaseClass(num_bin_target, num_bin_warpped, bg_value)
    {
    }

    template<typename ImageType> 
    hoImageRegDissimilarityNormalizedMutualInformation<ImageType>::~hoImageRegDissimilarityNormalizedMutualInformation()
    {
    }

    template<typename ImageType> 
    typename hoImageRegDissimilarityNormalizedMutualInformation<ImageType>::ValueType hoImageRegDissimilarityNormalizedMutualInformation<ImageType>::evaluate(ImageType& w)
    {
        try
        {
            BaseClass::evaluate(w);

            // convert to probabilities
            if ( num_samples_in_hist_ > 0 )
            {
                Gadgetron::scal((hist_value_type)(1.0/num_samples_in_hist_), hist_);
            }

            /// compute entorpy
            hist_target_.create(num_bin_target_);
            Gadgetron::clear(hist_target_);

            hist_warpped_.create(num_bin_warpped_);
            Gadgetron::clear(hist_warpped_);

            hist_.sumOverRow(hist_target_);
            hist_.sumOverCol(hist_warpped_);

            hist_value_type entropy_t(0), entropy_w(0), joint_entropy(0);

            size_t t, w;

            hist_value_type log2 = hist_value_type(1.0)/log( hist_value_type(2.0) );

            for ( t=0; t<(size_t)num_bin_target_; t++ )
            {
                hist_value_type prob = hist_target_(t);
                if ( prob > 0 )
                {
                    entropy_t -= prob * log(prob) * log2;
                }
            }

            for ( w=0; w<num_bin_warpped_; w++ )
            {
                hist_value_type prob = hist_warpped_(w);
                if ( prob > 0 )
                {
                    entropy_w -= prob * log(prob) * log2;
                }
            }

            for ( w=0; w<num_bin_warpped_; w++ )
            {
                for ( t=0; t<num_bin_target_; t++ )
                {
                    hist_value_type prob = hist_(t, w);
                    if ( prob > 0 )
                    {
                        joint_entropy -= prob * log(prob) * log2;
                    }
                }
            }

            dissimilarity_ = - (entropy_t + entropy_w) / joint_entropy;
        }
        catch(...)
        {
            GERROR_STREAM("Errors happened in hoImageRegDissimilarityNormalizedMutualInformation<ImageType>::evaluate(ImageType& t, ImageType& w) ... ");
        }

        return this->dissimilarity_;
    }

    template<typename ImageType> 
    void hoImageRegDissimilarityNormalizedMutualInformation<ImageType>::print(std::ostream& os) const
    {
        using namespace std;
        os << "--------------Gagdgetron image dissimilarity with histogram -------------" << endl;
        os << "Image dimension is : " << D << endl;

        std::string elemTypeName = std::string(typeid(ValueType).name());
        os << "Transformation data type is : " << elemTypeName << std::endl;

        os << "Number of intensity bins for target is : " << num_bin_target_ << endl;
        os << "Number of intensity bins for warped is : " << num_bin_warpped_ << endl;
        os << "PV interpolation for histogram is : " << pv_interpolation_ << endl;
        os << "Step size to ignore pixels when creating histogram is : " << step_size_ignore_pixel_ << endl;
    }
}
#endif // hoImageRegDissimilarityNormalizedMutualInformation_H_
