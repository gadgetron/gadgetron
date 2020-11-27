/** \file   hoImageRegDissimilarity.h
    \brief  Define the class to compute image dissimilarity in gadgetron registration

            Four different types of image dissimilarity measures are implemented here:

            SSD: sum-of-square difference
            LocalCCR: localized cross-correlation
            MI: mutual information
            NMI: normalized mutual information

            For  SSD, LocalCCR and MI, the analytical derivatives are computed.

            The analytical derivatives are computed by using the formula proposed at:

            [1] Gerardo Hermosillo, Christophe Chefd'Hotel, Olivier Faugeras. Variational Methods for Multimodal Image Matching. 
            International Journal of Computer Vision. December 2002, Volume 50, Issue 3, pp 329-343.
            http://link.springer.com/article/10.1023%2FA%3A1020830525823

            [2] Gerardo Hermosillo. Variational Methods for Multimodal Image Matching. PhD Thesis, UNIVERSIT´E DE NICE - SOPHIA ANTIPOLIS. May 2002.
            http://webdocs.cs.ualberta.ca/~dana/readingMedIm/papers/hermosilloPhD.pdf

            The derivative computation code is based on the listed source code at page 179 - 185 in ref [2] and extended.

            [3] Christophe Chefd'Hotel, Gerardo Hermosillo, Olivier D. Faugeras: Flows of diffeomorphisms for multimodal image registration. ISBI 2002: 753-756.
            http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=1029367&tag=1

            [4] C. Studholme, D.L.G. Hill, D.J. Hawkes. An overlap invariant entropy measure of 3D medical image alignment. Pattern Recognition, 32, 71-86, 1999.
            http://eecs.vanderbilt.edu/courses/cs359/other_links/papers/studholme_NMI_1999.pdf

    \author Hui Xue
*/

#ifndef hoImageRegDissimilarity_H_
#define hoImageRegDissimilarity_H_

#pragma once

#include "hoNDArray.h"
#include "hoNDImage.h"
#include "hoNDInterpolator.h"
#include "hoNDBoundaryHandler.h"
#include "hoMatrix.h"
#include "hoNDArray_utils.h"
#include "hoNDArray_elemwise.h"
#include "hoNDImage_util.h"

#ifdef USE_OMP
    #include <omp.h>
#endif // USE_OMP

namespace Gadgetron {

    // define the image dissimilarity type
    enum GT_IMAGE_DISSIMILARITY
    {
        GT_IMAGE_DISSIMILARITY_SSD,
        GT_IMAGE_DISSIMILARITY_LocalCCR,
        GT_IMAGE_DISSIMILARITY_MI,
        GT_IMAGE_DISSIMILARITY_NMI
    };

    inline std::string getDissimilarityName(GT_IMAGE_DISSIMILARITY v)
    {
        std::string name;

        switch (v)
        {
            case GT_IMAGE_DISSIMILARITY_SSD:
                name = "SSD";
                break;

            case GT_IMAGE_DISSIMILARITY_LocalCCR:
                name = "LocalCCR";
                break;

            case GT_IMAGE_DISSIMILARITY_MI:
                name = "MutualInformation";
                break;

            case GT_IMAGE_DISSIMILARITY_NMI:
                name = "NormalizedMutualInformation";
                break;

            default:
                GERROR_STREAM("Unrecognized image dissimilarity type : " << v);
        }

        return name;
    }

    inline GT_IMAGE_DISSIMILARITY getDissimilarityType(const std::string& name)
    {
        GT_IMAGE_DISSIMILARITY v;

        if ( name == "SSD" )
        {
            v = GT_IMAGE_DISSIMILARITY_SSD;
        }
        else if ( name == "LocalCCR" )
        {
            v = GT_IMAGE_DISSIMILARITY_LocalCCR;
        }
        else if ( name == "MutualInformation" )
        {
            v = GT_IMAGE_DISSIMILARITY_MI;
        }
        else if ( name == "NormalizedMutualInformation" )
        {
            v = GT_IMAGE_DISSIMILARITY_NMI;
        }
        else
        {
            GERROR_STREAM("Unrecognized image dissimilarity name : " << name);
        }

        return v;
    }

    /// compute the image dissimilarity measures
    /// if possible, compute the analytical derivatives
    template<typename ImageType> 
    class hoImageRegDissimilarity
    {
    public:

        typedef hoImageRegDissimilarity<ImageType> Self;
        typedef hoNDInterpolator<ImageType> InterpolatorType;

        typedef typename ImageType::value_type ValueType;
        enum { D = ImageType::NDIM };

        typedef ValueType T;
        typedef ValueType element_type;
        typedef ValueType value_type;

        typedef float coord_type;

        hoImageRegDissimilarity(ValueType bg_value=ValueType(0));
        virtual ~hoImageRegDissimilarity();

        /// initialize the dissimilarity
        virtual void initialize(ImageType& t);

        const ImageType& getDeriv() const { return deriv_; }

        ValueType getDissimilarity() const { return dissimilarity_; }

        void setBackgroundValue(ValueType bg_value) { bg_value_ = bg_value; }

        /// compute the dissimilarity value
        virtual ValueType evaluate(ImageType& w);

        /// compute the derivative and dissimilarity value
        virtual bool evaluateDeriv(ImageType& w) = 0;

        virtual void print(std::ostream& os) const;

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

    protected:

        ImageType* target_;
        ImageType* warpped_;

        std::vector<size_t> image_dim_;

        /// background pixels
        ValueType bg_value_;

        /// derivative to spatial locations
        ImageType deriv_;

        /// dissimilarity value
        ValueType dissimilarity_;

        hoNDArray<ValueType> target;
        hoNDArray<ValueType> warped;
        hoNDArray<ValueType> deriv;
    };

    template<typename ImageType> 
    hoImageRegDissimilarity<ImageType>::hoImageRegDissimilarity(ValueType bg_value) 
        : target_(NULL), warpped_(NULL), bg_value_(bg_value), dissimilarity_(0), performTiming_(false)
    {
        gt_timer1_.set_timing_in_destruction(false);
        gt_timer2_.set_timing_in_destruction(false);
        gt_timer3_.set_timing_in_destruction(false);
    }

    template<typename ImageType> 
    hoImageRegDissimilarity<ImageType>::~hoImageRegDissimilarity()
    {
    }

    template<typename ImageType> 
    void hoImageRegDissimilarity<ImageType>::initialize(ImageType& t)
    {
        target_ = &t;

        if ( !deriv_.dimensions_equal(*target_) )
        {
            deriv_.create(target_->get_dimensions());
        }
        memset( deriv_.get_data_ptr(), 0, deriv_.get_number_of_elements()*sizeof(ValueType));

        target_->get_dimensions(image_dim_);

        /// these conversion can be removed if more utility functions are added for hoNDImage
        target.create(image_dim_, target_->begin(), false);
        deriv.create(image_dim_, deriv_.begin(), false);
    }

    template<typename ImageType> 
    typename hoImageRegDissimilarity<ImageType>::ValueType hoImageRegDissimilarity<ImageType>::evaluate(ImageType& w)
    {
        if ( warpped_ != &w )
        {
            warpped_ = &w;
            GADGET_CHECK_THROW(warpped_->dimensions_equal(*target_));
            warped.create(image_dim_, warpped_->begin(), false);
        }

        this->dissimilarity_ = 0;

        return this->dissimilarity_;
    }

    template<typename ImageType> 
    void hoImageRegDissimilarity<ImageType>::print(std::ostream& os) const
    {
        using namespace std;
        os << "--------------Gagdgetron image dissimilarity measure -------------" << endl;
        os << "Image dimension is : " << D << endl;

        std::string elemTypeName = std::string(typeid(ValueType).name());
        os << "Transformation data type is : " << elemTypeName << endl << ends;
    }
}
#endif // hoImageRegDissimilarity_H_
