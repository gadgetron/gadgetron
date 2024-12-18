/** \file   hoImageRegDissimilarityLocalCCR.h
    \brief  Define the class to compute image Local Cross CorRelation (LocalCCR) in gadgetron registration

            The analytical derivatives are computed by using the formula proposed at:

            [1] Gerardo Hermosillo, Christophe Chefd'Hotel, Olivier Faugeras. Variational Methods for Multimodal Image Matching. 
            International Journal of Computer Vision. December 2002, Volume 50, Issue 3, pp 329-343.
            http://link.springer.com/article/10.1023%2FA%3A1020830525823

            [2] Gerardo Hermosillo. Variational Methods for Multimodal Image Matching. PhD Thesis, UNIVERSIT´E DE NICE - SOPHIA ANTIPOLIS. May 2002.
            http://webdocs.cs.ualberta.ca/~dana/readingMedIm/papers/hermosilloPhD.pdf

            This derivative computation code is based on the listed source code at page 183 - 185 in ref [2].

    \author Hui Xue
*/

#ifndef hoImageRegDissimilarityLocalCCR_H_
#define hoImageRegDissimilarityLocalCCR_H_

#pragma once

#include <limits>
#include "hoImageRegDissimilarity.h"

namespace Gadgetron {

    template<typename ImageType> 
    class hoImageRegDissimilarityLocalCCR : public hoImageRegDissimilarity<ImageType>
    {
    public:

        typedef hoImageRegDissimilarityLocalCCR<ImageType> Self;
        typedef hoImageRegDissimilarity<ImageType> BaseClass;

        enum { D = ImageType::NDIM };

        typedef typename BaseClass::InterpolatorType InterpolatorType;
        typedef typename BaseClass::ValueType ValueType;
        typedef ValueType T;
        typedef ValueType element_type;
        typedef ValueType value_type;

        typedef double computing_value_type;

        typedef typename BaseClass::coord_type coord_type;

        hoImageRegDissimilarityLocalCCR(computing_value_type betaArg=std::numeric_limits<ValueType>::epsilon() );
        hoImageRegDissimilarityLocalCCR(ValueType sigmaArg[D], computing_value_type betaArg=std::numeric_limits<ValueType>::epsilon() );
        virtual ~hoImageRegDissimilarityLocalCCR();

        void initialize(ImageType& t);

        virtual ValueType evaluate(ImageType& w);
        virtual bool evaluateDeriv(ImageType& w);

        virtual void print(std::ostream& os) const;

        /// these parameter names are kept same as the source code on page 183 - 185 in ref [2]
        computing_value_type sigmaArg_[D]; // kernel size of local weighting function

        computing_value_type betaArg_;

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

        /// these parameter names are kept same as the source code on page 183 - 185 in ref [2]
        hoNDArray<computing_value_type> cc; computing_value_type* p_cc;
        hoNDArray<computing_value_type> mu1; computing_value_type* p_mu1;
        hoNDArray<computing_value_type> mu2; computing_value_type* p_mu2;
        hoNDArray<computing_value_type> v1; computing_value_type* p_v1;
        hoNDArray<computing_value_type> v2; computing_value_type* p_v2;
        hoNDArray<computing_value_type> v12; computing_value_type* p_v12;

        //hoNDArray<computing_value_type> vv1; computing_value_type* p_vv1;
        //hoNDArray<computing_value_type> vv2; computing_value_type* p_vv2;
        //hoNDArray<computing_value_type> vv12; computing_value_type* p_vv12;

        hoNDArray<computing_value_type> mem_;

        computing_value_type eps_;
    };

    template<typename ImageType> 
    hoImageRegDissimilarityLocalCCR<ImageType>::hoImageRegDissimilarityLocalCCR(computing_value_type betaArg) 
        : BaseClass(), betaArg_(betaArg)
    {
        unsigned int ii;
        for ( ii=0; ii<D; ii++ )
        {
            sigmaArg_[ii] = (computing_value_type)(2.0);
        }
    }

    template<typename ImageType> 
    hoImageRegDissimilarityLocalCCR<ImageType>::hoImageRegDissimilarityLocalCCR(ValueType sigmaArg[D], computing_value_type betaArg) 
        : BaseClass(), betaArg_(betaArg)
    {
        unsigned int ii;
        for ( ii=0; ii<D; ii++ )
        {
            sigmaArg_[ii] = (computing_value_type)(sigmaArg[ii]);
        }
    }

    template<typename ImageType> 
    hoImageRegDissimilarityLocalCCR<ImageType>::~hoImageRegDissimilarityLocalCCR()
    {
    }

    template<typename ImageType> 
    void hoImageRegDissimilarityLocalCCR<ImageType>::initialize(ImageType& t)
    {
        BaseClass::initialize(t);

        // allocate arrays for the computation
        cc.create(image_dim_); p_cc = cc.begin();
        mu1.create(image_dim_); p_mu1 = mu1.begin();
        mu2.create(image_dim_); p_mu2 = mu2.begin();
        v1.create(image_dim_); p_v1 = v1.begin();
        v2.create(image_dim_); p_v2 = v2.begin();
        v12.create(image_dim_); p_v12 = v12.begin();

        //vv1.create(image_dim_); p_vv1 = vv1.begin();
        //vv2.create(image_dim_); p_vv2 = vv2.begin();
        //vv12.create(image_dim_); p_vv12 = vv12.begin();

        eps_ = std::numeric_limits<computing_value_type>::epsilon();
    }

    template<typename ImageType> 
    typename hoImageRegDissimilarityLocalCCR<ImageType>::ValueType hoImageRegDissimilarityLocalCCR<ImageType>::evaluate(ImageType& w)
    {
        try
        {
            /// in the ref [2], the code are:
            /*
            Image<float>
                mu1(I1.domain()), mu2(I1.domain()),
                v1(I1.domain()), v2(I1.domain()),
                v12(I1.domain()), f1(I1.domain()),
                f2(I1.domain()), f3(I1.domain());
                Map(I1,x) {
                const real i1 = I1[x];
                const real i2 = I2[x];
                mu1[x] = i1; v1[x] = i1 * i1;
                mu2[x] = i2; v12[x] = i1 * i2;
                v2[x] = i2 * i2;
                }
                mu1.SelfRecSmoothZeroBC(sigma); v1.SelfRecSmoothZeroBC(sigma);
                mu2.SelfRecSmoothZeroBC(sigma); v2.SelfRecSmoothZeroBC(sigma);
                v12.SelfRecSmoothZeroBC(sigma);

                criter = 0;
                Map(v1,x) {
                const real u1 = mu1[x];
                const real u2 = mu2[x];
                const real vv1 = v1[x] + beta - u1 * u1;
                const real vv2 = v2[x] + beta - u2 * u2;
                const real vv12 = v12[x] - u1 * u2;
                const real ff1 = vv12 / (vv1 * vv2);
                const real CC = vv12 * ff1;
                const real ff2 = - CC / vv2;
                const real ff3 =  - (ff2 * u2 + ff1 * u1);
                f1[x] = ff1; f2[x] = ff2; f3[x] = ff3;
                cc[x] = -CC;
                criter += -CC;
                }
                f1.SelfRecSmoothZeroBC(sigma);
                f2.SelfRecSmoothZeroBC(sigma);
                f3.SelfRecSmoothZeroBC(sigma);

                norm = 0;
                Map(f1,x) {
                const float val = 2.0 * ( f1[x] * I1[x] + f2[x] * I2[x] + f3[x] ) ;
                dist[x] = val;
                norm += val * val;
                }
            */

            /// we rewrite these code for gadgetron

            //if ( performTiming_ ) { gt_timer1_.start("1"); }
            BaseClass::evaluate(w);
            //if ( performTiming_ ) { gt_timer1_.stop(); }

            long long N = (long long)target.get_number_of_elements();

            //if ( performTiming_ ) { gt_timer1_.start("2"); }
            //mu1.copyFrom(target);
            //mu2.copyFrom(warped);
            //Gadgetron::multiply(mu1, mu1, v1);
            //Gadgetron::multiply(mu2, mu2, v2);
            //Gadgetron::multiply(mu1, mu2, v12);

            long long n;

            ValueType* pT = target.begin();
            ValueType* pW = warped.begin();

            for ( n=0; n<N; ++n )
            {
                const computing_value_type v1 = (computing_value_type)pT[n];
                const computing_value_type v2 = (computing_value_type)pW[n];

                p_mu1[n] = v1;
                p_mu2[n] = v2;
                p_v1[n] = v1*v1;
                p_v2[n] = v2*v2;
                p_v12[n] = v1*v2;
            }

            Gadgetron::filterGaussian(mu1, sigmaArg_, mem_.begin());
            Gadgetron::filterGaussian(mu2, sigmaArg_, mem_.begin());
            Gadgetron::filterGaussian(v1, sigmaArg_, mem_.begin());
            Gadgetron::filterGaussian(v2, sigmaArg_, mem_.begin());
            Gadgetron::filterGaussian(v12, sigmaArg_, mem_.begin());

            //if ( 0 )
            //{
            //    //#pragma omp parallel sections if ( D==2 )
            //    {
            //        //#pragma omp section
            //        {
            //            Gadgetron::multiply(mu1, mu1, vv1);
            //            Gadgetron::subtract(v1, vv1, vv1);
            //            Gadgetron::addEpsilon(vv1);
            //        }

            //        //#pragma omp section
            //        {
            //            Gadgetron::multiply(mu2, mu2, vv2);
            //            Gadgetron::subtract(v2, vv2, vv2);
            //            Gadgetron::addEpsilon(vv2);
            //        }

            //        //#pragma omp section
            //        {
            //            Gadgetron::multiply(mu1, mu2, vv12);
            //            Gadgetron::subtract(v12, vv12, vv12);
            //        }
            //    }

            //    Gadgetron::multiply(vv1, vv2, vv1);
            //    Gadgetron::divide(vv12, vv1, v1); // ff1

            //    Gadgetron::multiply(vv12, v1, cc); // cc

            //    Gadgetron::divide(cc, vv2, v2); // ff2
            //    Gadgetron::scal( (computing_value_type)(-1), v2);

            //    Gadgetron::multiply(v2, mu2, v12);
            //    Gadgetron::multiply(v1, mu1, vv12);
            //    Gadgetron::add(v12, vv12, v12);

            //    computing_value_type v=0;
            //    Gadgetron::norm1(cc, v);

            //    dissimilarity_ = static_cast<T>(-v/N);
            //}

            dissimilarity_ = 0;
            computing_value_type v=0;

            //#pragma omp parallel for private(n)
            for ( n=0; n<N; ++n )
            {
                const computing_value_type u1 = p_mu1[n];
                const computing_value_type u2 = p_mu2[n];

                const computing_value_type vv1 = p_v1[n] - u1 * u1;
                const computing_value_type vv2 = p_v2[n] - u2 * u2;
                const computing_value_type vv12 = p_v12[n] - u1 * u2;

                const computing_value_type ff1 = vv12 / (vv1 * vv2);
                const computing_value_type lcc = vv12 * ff1;

                const computing_value_type ff2 = - lcc / vv2;
                const computing_value_type ff3 = ff2 * u2 + ff1 * u1;

                p_v1[n] = ff1; p_v2[n] = ff2; p_v12[n] = ff3;

                p_cc[n] = lcc;
            }

            computing_value_type lcc = 0;

            // #pragma omp parallel for reduction(+:lcc)
            for (n=0; n<N; n++)
            {
                lcc += cc[n];
            }

            dissimilarity_ = -lcc/N;
        }
        catch(...)
        {
            GERROR_STREAM("Errors happened in hoImageRegDissimilarityLocalCCR<ImageType>::evaluate(w) ... ");
        }

        return this->dissimilarity_;
    }

    template<typename ImageType> 
    bool hoImageRegDissimilarityLocalCCR<ImageType>::evaluateDeriv(ImageType& w)
    {
        try
        {
            this->evaluate(w);

            size_t N = target.get_number_of_elements();

            long long n;

            //#pragma omp parallel sections if ( D==2 )
            {
                //#pragma omp section
                {
                    Gadgetron::filterGaussian(v1, sigmaArg_, mem_.begin());
                }

                //#pragma omp section
                {
                    Gadgetron::filterGaussian(v2, sigmaArg_, mem_.begin());
                }

                //#pragma omp section
                {
                    Gadgetron::filterGaussian(v12, sigmaArg_, mem_.begin());
                }
            }

            // deriv = f1*i1 + f2*i2 + f3, we don't need to multiply this by 2.0

            //if ( typeid(ValueType) == typeid(computing_value_type) )
            //{
                //Gadgetron::multiply(v1, target, mu1);
                //Gadgetron::multiply(v2, warped, mu2);
                //Gadgetron::add(mu1, mu2, deriv);
                //Gadgetron::subtract(deriv, v12, deriv);
            //}
            //else
            //{
                T* pT = target.begin();
                T* pW = warped.begin();

                // #pragma omp parallel for default(none) shared(N, pT, pW)
                for ( n=0; n<(long long)N; n++ )
                {
                    deriv(n) = static_cast<T>( p_v1[n]* (computing_value_type)pT[n] + ( p_v2[n]*(computing_value_type)pW[n] - p_v12[n] ) );
                }
            //}
        }
        catch(...)
        {
            GERROR_STREAM("Errors happened in hoImageRegDissimilarityLocalCCR<ImageType>::evaluateDeriv(w) ... ");
            return false;
        }

        return true;
    }

    template<typename ImageType> 
    void hoImageRegDissimilarityLocalCCR<ImageType>::print(std::ostream& os) const
    {
        using namespace std;
        os << "--------------Gagdgetron image dissimilarity LocalCCR measure -------------" << endl;
        os << "Image dimension is : " << D << endl;

        std::string elemTypeName = std::string(typeid(ValueType).name());
        os << "Transformation data type is : " << elemTypeName << endl << ends;
    }
}
#endif // hoImageRegDissimilarityLocalCCR_H_
