/** \file   hoImageRegDissimilaritySSD.h
    \brief  Define the class to compute image sum-of-square difference (SSD ) in gadgetron registration

            The analytical derivatives are computed by using the formula proposed at:

            [1] Gerardo Hermosillo, Christophe Chefd'Hotel, Olivier Faugeras. Variational Methods for Multimodal Image Matching. 
            International Journal of Computer Vision. December 2002, Volume 50, Issue 3, pp 329-343.
            http://link.springer.com/article/10.1023%2FA%3A1020830525823

            [2] Gerardo Hermosillo. Variational Methods for Multimodal Image Matching. PhD Thesis, UNIVERSIT�E DE NICE - SOPHIA ANTIPOLIS. May 2002.
            http://webdocs.cs.ualberta.ca/~dana/readingMedIm/papers/hermosilloPhD.pdf

            The derivative computation code is modified from the listed source code at page 179 - 185 in ref [2].

    \author Hui Xue
*/

#ifndef hoImageRegDissimilaritySSD_H_
#define hoImageRegDissimilaritySSD_H_

#pragma once

#include "hoImageRegDissimilarity.h"

namespace Gadgetron {

    template<typename ImageType> 
    class hoImageRegDissimilaritySSD : public hoImageRegDissimilarity<ImageType>
    {
    public:

        typedef hoImageRegDissimilaritySSD<ImageType> Self;
        typedef hoImageRegDissimilarity<ImageType> BaseClass;

        enum { D = ImageType::NDIM };

        typedef typename BaseClass::InterpolatorType InterpolatorType;
        typedef typename BaseClass::ValueType ValueType;
        typedef ValueType T;
        typedef ValueType element_type;
        typedef ValueType value_type;

        typedef typename BaseClass::coord_type coord_type;

        hoImageRegDissimilaritySSD();
        virtual ~hoImageRegDissimilaritySSD();

        virtual ValueType evaluate(ImageType& w);
        virtual bool evaluateDeriv(ImageType& w) { this->evaluate(w); return true; }

        virtual void print(std::ostream& os) const;

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
    };

    template<typename ImageType> 
    hoImageRegDissimilaritySSD<ImageType>::hoImageRegDissimilaritySSD() : BaseClass()
    {
    }

    template<typename ImageType> 
    hoImageRegDissimilaritySSD<ImageType>::~hoImageRegDissimilaritySSD()
    {
    }

    template<typename ImageType> 
    typename hoImageRegDissimilaritySSD<ImageType>::ValueType hoImageRegDissimilaritySSD<ImageType>::evaluate(ImageType& w)
    {
        try
        {
            BaseClass::evaluate(w);

            Gadgetron::subtract(target, warped, deriv);
            dissimilarity_ = Gadgetron::nrm2(deriv);

            dissimilarity_ = (dissimilarity_*dissimilarity_) / (ValueType)(target.get_number_of_elements());
        }
        catch(...)
        {
            GERROR_STREAM("Errors happened in hoImageRegDissimilaritySSD<ImageType>::evaluate(w) ... ");
        }

        return this->dissimilarity_;
    }

    template<typename ImageType> 
    void hoImageRegDissimilaritySSD<ImageType>::print(std::ostream& os) const
    {
        using namespace std;
        os << "--------------Gagdgetron image dissimilarity SSD measure -------------" << endl;
        os << "Image dimension is : " << D << endl;

        std::string elemTypeName = std::string(typeid(ValueType).name());
        os << "Transformation data type is : " << elemTypeName << std::endl;
    }
}
#endif // hoImageRegDissimilaritySSD_H_
