/** \file   hoImageRegDissimilaritySSD.h
    \brief  Define the class to compute image sum-of-square difference (SSD ) in gadgetron registration

            The analytical derivatives are computed by using the formula proposed at:

            [1] Gerardo Hermosillo, Christophe Chefd'Hotel, Olivier Faugeras. Variational Methods for Multimodal Image Matching. 
            International Journal of Computer Vision. December 2002, Volume 50, Issue 3, pp 329-343.
            http://link.springer.com/article/10.1023%2FA%3A1020830525823

            [2] Gerardo Hermosillo. Variational Methods for Multimodal Image Matching. PhD Thesis, UNIVERSIT´E DE NICE - SOPHIA ANTIPOLIS. May 2002.
            http://webdocs.cs.ualberta.ca/~dana/readingMedIm/papers/hermosilloPhD.pdf

            The derivative computation code is modified from the listed source code at page 179 - 185 in ref [2].

    \author Hui Xue
*/

#pragma once

#include "hoImageRegDissimilarity.h"

namespace Gadgetron
{
    template<typename ValueType, unsigned int D> 
    class hoImageRegDissimilaritySSD : public hoImageRegDissimilarity<ValueType, D>
    {
    public:

        typedef hoImageRegDissimilaritySSD<ValueType, D> Self;
        typedef hoImageRegDissimilarity<ValueType, D> BaseClass;

        typedef typename BaseClass::ImageType ImageType;
        typedef typename BaseClass::InterpolatorType InterpolatorType;

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

    template<typename ValueType, unsigned int D> 
    hoImageRegDissimilaritySSD<ValueType, D>::hoImageRegDissimilaritySSD() : BaseClass()
    {
    }

    template<typename ValueType, unsigned int D> 
    hoImageRegDissimilaritySSD<ValueType, D>::~hoImageRegDissimilaritySSD()
    {
    }

    template<typename ValueType, unsigned int D> 
    ValueType hoImageRegDissimilaritySSD<ValueType, D>::evaluate(ImageType& w)
    {
        try
        {
            BaseClass::evaluate(w);

            GADGET_CHECK_RETURN_FALSE(Gadgetron::subtract(target, warped, deriv));
            GADGET_CHECK_RETURN_FALSE(Gadgetron::norm2(deriv, dissimilarity_));

            dissimilarity_ = (dissimilarity_*dissimilarity_) / (ValueType)(target.get_number_of_elements());
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors happened in hoImageRegDissimilaritySSD<ValueType, D>::evaluate(w) ... ");
        }

        return this->dissimilarity_;
    }

    template<typename ValueType, unsigned int D> 
    void hoImageRegDissimilaritySSD<ValueType, D>::print(std::ostream& os) const
    {
        using namespace std;
        os << "--------------Gagdgetron image dissimilarity SSD measure -------------" << endl;
        os << "Image dimension is : " << D << endl;

        std::string elemTypeName = std::string(typeid(ValueType).name());
        os << "Transformation data type is : " << elemTypeName << std::endl;
    }
}
