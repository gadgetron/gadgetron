/** \file   hoImageRegDissimilarityHistogramBased.h
    \brief  Define the class to compute image dissimilarity based on histogram
    \author Hui Xue
*/

#ifndef hoImageRegDissimilarityHistogramBased_H_
#define hoImageRegDissimilarityHistogramBased_H_

#pragma once

#include <limits>
#include "hoMatrix.h"
#include "hoImageRegDissimilarity.h"

namespace Gadgetron {

    template<typename ImageType> 
    class hoImageRegDissimilarityHistogramBased : public hoImageRegDissimilarity<ImageType>
    {
    public:

        typedef hoImageRegDissimilarityHistogramBased<ImageType> Self;
        typedef hoImageRegDissimilarity<ImageType> BaseClass;

        enum { D = ImageType::NDIM };

        typedef typename BaseClass::InterpolatorType InterpolatorType;
        typedef typename BaseClass::ValueType ValueType;
        typedef ValueType T;
        typedef ValueType element_type;
        typedef ValueType value_type;

        typedef typename BaseClass::coord_type coord_type;

        typedef double hist_value_type;

        hoImageRegDissimilarityHistogramBased(unsigned int num_bin_target=64, unsigned int num_bin_warpped=64, ValueType bg_value=ValueType(0));
        virtual ~hoImageRegDissimilarityHistogramBased();

        virtual ValueType evaluate(ImageType& w);

        virtual bool evaluateDeriv(ImageType& w) = 0;

        virtual void print(std::ostream& os) const;

        /// number of intensity bins
        unsigned int num_bin_target_;
        unsigned int num_bin_warpped_;

        /// whether to perform partial interpolation for histogram
        bool pv_interpolation_;

        /// step size to ignore pixels when creating histogram
        size_t step_size_ignore_pixel_;

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

        /// store the 2D histogram
        hoMatrix<hist_value_type> hist_;

        /// min/max intensities of target and warped
        ValueType min_target_;
        ValueType max_target_;

        ValueType min_warpped_;
        ValueType max_warpped_;

        size_t num_samples_in_hist_;
    };

    template<typename ImageType> 
    hoImageRegDissimilarityHistogramBased<ImageType>::
    hoImageRegDissimilarityHistogramBased(unsigned int num_bin_target, unsigned int num_bin_warpped, ValueType bg_value) 
        : BaseClass(bg_value), num_bin_target_(num_bin_target), num_bin_warpped_(num_bin_warpped), pv_interpolation_(false), step_size_ignore_pixel_(1)
    {
    }

    template<typename ImageType> 
    hoImageRegDissimilarityHistogramBased<ImageType>::~hoImageRegDissimilarityHistogramBased()
    {
    }

    template<typename ImageType> 
    typename hoImageRegDissimilarityHistogramBased<ImageType>::ValueType hoImageRegDissimilarityHistogramBased<ImageType>::evaluate(ImageType& w)
    {
        try
        {
            BaseClass::evaluate(w);

            // allocate histogram
            hist_.createMatrix(num_bin_target_, num_bin_warpped_);
            Gadgetron::clear(hist_);

            // intensity range
            min_target_ = std::numeric_limits<ValueType>::max();
            max_target_ = std::numeric_limits<ValueType>::min();

            min_warpped_ = min_target_;
            max_warpped_ = max_target_;

            size_t N = target_->get_number_of_elements();

            long long n;
            for ( n=0; n<(long long)N; n++ )
            {
                ValueType vt = target(n);
                if ( vt < min_target_ ) min_target_ = vt;
                if ( vt > max_target_ ) max_target_ = vt;

                ValueType vw = warped(n);
                if ( vw < min_warpped_ ) min_warpped_ = vw;
                if ( vw > max_warpped_ ) max_warpped_ = vw;
            }

            ValueType range_t = ValueType(1.0)/(max_target_ - min_target_ + std::numeric_limits<ValueType>::epsilon());
            ValueType range_w = ValueType(1.0)/(max_warpped_ - min_warpped_ + std::numeric_limits<ValueType>::epsilon());

            num_samples_in_hist_ = 0;

            if ( pv_interpolation_ )
            {
                #pragma omp parallel for default(none) private(n) shared(N, range_t, range_w)
                for ( n=0; n<(long long)N; n+=(long long)step_size_ignore_pixel_ )
                {
                    ValueType vt = target(n);
                    ValueType vw = warped(n);

                    if ( std::abs(vt-bg_value_)<FLT_EPSILON 
                        && std::abs(vw-bg_value_)<FLT_EPSILON )
                    {
                        continue;
                    }

                    ValueType xT = range_t*(vt-min_target_)*(num_bin_target_-1);
                    ValueType xW = range_w*(vw-min_warpped_)*(num_bin_warpped_-1);

                    size_t indT = static_cast<size_t>(xT);
                    size_t indW = static_cast<size_t>(xW);

                    ValueType sT, s1T, sW, s1W;

                    sT = xT - indT; s1T = 1 - sT;
                    sW = xW - indW; s1W = 1 - sW;

                    #pragma omp critical
                    {
                        hist_(indT, indW) += s1T*s1W;

                        if ( indT<num_bin_target_-1 && indW<num_bin_warpped_-1 )
                        {
                            hist_(indT, indW+1) += s1T*sW;
                            hist_(indT+1, indW) += sT*s1W;
                            hist_(indT+1, indW+1) += sT*sW;
                        }
                    }

                    #pragma omp atomic
                    num_samples_in_hist_++;
                }
            }
            else
            {
                #pragma omp parallel for default(none) private(n) shared(N, range_t, range_w)
                for ( n=0; n<(long long)N; n+=(long long)step_size_ignore_pixel_ )
                {
                    ValueType vt = target(n);
                    ValueType vw = warped(n);

                    if ( std::abs(vt-bg_value_)<FLT_EPSILON 
                        && std::abs(vw-bg_value_)<FLT_EPSILON )
                    {
                        continue;
                    }

                    size_t indT = static_cast<size_t>( range_t*(vt-min_target_)*(num_bin_target_-1) + 0.5 );
                    size_t indW = static_cast<size_t>( range_w*(vw-min_warpped_)*(num_bin_warpped_-1) + 0.5 );

                    #pragma omp critical
                    {
                        hist_(indT, indW)++;
                    }

                    #pragma omp atomic
                    num_samples_in_hist_++;
                }
            }

            if ( !debugFolder_.empty() ) {  gt_exporter_.export_array(hist_, debugFolder_+"hist2D"); }
        }
        catch(...)
        {
            GERROR_STREAM("Errors happened in hoImageRegDissimilarityHistogramBased<ImageType>::evaluate(ImageType& t, ImageType& w) ... ");
        }

        return this->dissimilarity_;
    }

    template<typename ImageType> 
    void hoImageRegDissimilarityHistogramBased<ImageType>::print(std::ostream& os) const
    {
        using namespace std;
        os << "--------------Gagdgetron image dissimilarity with histogram -------------" << endl;
        os << "Image dimension is : " << D << endl;

        std::string elemTypeName = std::string(typeid(ValueType).name());
        os << "Transformation data type is : " << elemTypeName << std::endl;

        os << "Number of intensity bins for target is : " << num_bin_target_ << endl;
        os << "Number of intensity bins for warped is : " << num_bin_warpped_ << endl;
        os << "PV interpolation for histogram is : " << pv_interpolation_ << endl;
        os << "Step size to ignore pixels when creating histogram is : " << step_size_ignore_pixel_ << endl << ends;
    }
}
#endif // hoImageRegDissimilarityHistogramBased_H_
