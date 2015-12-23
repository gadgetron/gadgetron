/** \file       gtPlusWavelet2DOperator.h
    \brief      Implement 2D wavelet operator for L1 regularization
    \author     Hui Xue

    Redundant haar wavelet transformation is implemented here.
*/

#pragma once

#include "gtPlusWaveletOperator.h"

namespace Gadgetron { namespace gtPlus {

template <typename T> 
class gtPlusWavelet2DOperator : public gtPlusWaveletOperator<T>
{
public:

    typedef gtPlusWaveletOperator<T> BaseClass;
    typedef typename BaseClass::value_type value_type;

    gtPlusWavelet2DOperator();
    virtual ~gtPlusWavelet2DOperator();

    virtual void printInfo(std::ostream& os);

    // forward operator, perform wavelet transform
    // x: [RO E1 ...]
    // y: [RO E1 W ...]
    virtual bool forwardOperator(const hoNDArray<T>& x, hoNDArray<T>& y);

    // adjoint operator, perform inverse transform
    // x: [RO E1 W ...]
    // y: [RO E1 ...]
    virtual bool adjointOperator(const hoNDArray<T>& x, hoNDArray<T>& y);

    virtual bool unitary() const { return true; }

    using BaseClass::scale_factor_first_dimension_;
    using BaseClass::scale_factor_second_dimension_;
    using BaseClass::numOfWavLevels_;
    using BaseClass::with_approx_coeff_;
    using BaseClass::gt_timer1_;
    using BaseClass::gt_timer2_;
    using BaseClass::gt_timer3_;
    using BaseClass::performTiming_;
    using BaseClass::gt_exporter_;
    using BaseClass::debugFolder_;
    using BaseClass::gtPlus_util_;
    using BaseClass::gtPlus_util_complex_;

protected:

    using BaseClass::acquired_points_;
    using BaseClass::acquired_points_indicator_;
    using BaseClass::unacquired_points_indicator_;
    using BaseClass::coil_senMap_;

    // helper memory
    using BaseClass::kspace_;
    using BaseClass::complexIm_;
    using BaseClass::res_after_apply_kernel_;
    using BaseClass::res_after_apply_kernel_sum_over_;

    using BaseClass::wav_coeff_norm_;
    using BaseClass::kspace_wav_;
    using BaseClass::complexIm_wav_;

    using BaseClass::kspace_Managed_;
    using BaseClass::complexIm_Managed_;
    using BaseClass::res_after_apply_kernel_Managed_;
    using BaseClass::res_after_apply_kernel_sum_over_Managed_;
};

template <typename T> 
gtPlusWavelet2DOperator<T>::gtPlusWavelet2DOperator() : BaseClass()
{

}

template <typename T> 
gtPlusWavelet2DOperator<T>::~gtPlusWavelet2DOperator()
{
}

template <typename T> 
bool gtPlusWavelet2DOperator<T>::
forwardOperator(const hoNDArray<T>& x, hoNDArray<T>& y)
{
    try
    {
        boost::shared_ptr< std::vector<size_t> > dims = x.get_dimensions();
        size_t NDim = dims->size();

        size_t RO = (*dims)[0];
        size_t E1 = (*dims)[1];
        size_t W = 1+3*numOfWavLevels_;

        std::vector<size_t> dimR(NDim+1);
        dimR[0] = RO;
        dimR[1] = E1;
        dimR[2] = W;

        size_t n;
        for ( n=2; n<NDim; n++ )
        {
            dimR[n+1] = (*dims)[n];
        }

        if ( !y.dimensions_equal(&dimR) )
        {
            y.create(&dimR);
        }

        size_t num = x.get_number_of_elements()/(RO*E1);

        T* pX = const_cast<T*>(x.begin());
        T* pY = y.begin();

        long long t;

        #pragma omp parallel for default(none) private(t) shared(num, RO, E1, W, pX, pY)
        for ( t=0; t<(long long)num; t++ )
        {
            hoNDArray<T> in(RO, E1, pX+t*RO*E1);
            hoNDArray<T> out(RO, E1, W, pY+t*RO*E1*W);
            // this->dwtRedundantHaar(in, out, numOfWavLevels_);

            Gadgetron::hoNDHarrWavelet<T> wav;
            wav.transform(in, out, 2, numOfWavLevels_, true);
        }
    }
    catch (...)
    {
        GERROR_STREAM("Errors in gtPlusWavelet2DOperator<T>::forwardOperator(const hoNDArray<T>& x, hoNDArray<T>& y) ... ");
        return false;
    }
    return true;
}

template <typename T> 
bool gtPlusWavelet2DOperator<T>::
adjointOperator(const hoNDArray<T>& x, hoNDArray<T>& y)
{
    try
    {
        boost::shared_ptr< std::vector<size_t> > dims = x.get_dimensions();
        size_t NDim = dims->size();

        size_t RO = (*dims)[0];
        size_t E1 = (*dims)[1];
        size_t W = (*dims)[2];

        std::vector<size_t> dimR(NDim-1);
        dimR[0] = RO;
        dimR[1] = E1;

        size_t n;
        for ( n=2; n<NDim-1; n++ )
        {
            dimR[n] = (*dims)[n+1];
        }

        if ( !y.dimensions_equal(&dimR) )
        {
            y.create(&dimR);
        }

        size_t num = x.get_number_of_elements()/(RO*E1*W);

        T* pX = const_cast<T*>(x.begin());
        T* pY = y.begin();

        long long t;

        #pragma omp parallel for default(none) private(t) shared(num, RO, E1, W, pX, pY)
        for ( t=0; t<(long long)num; t++ )
        {
            hoNDArray<T> in(RO, E1, W, pX+t*RO*E1*W);
            hoNDArray<T> out(RO, E1, pY+t*RO*E1);
            // this->idwtRedundantHaar(in, out, numOfWavLevels_);

            Gadgetron::hoNDHarrWavelet<T> wav;
            wav.transform(in, out, 2, numOfWavLevels_, false);
        }
    }
    catch (...)
    {
        GERROR_STREAM("Errors in gtPlusWavelet2DOperator<T>::adjointOperator(const hoNDArray<T>& x, hoNDArray<T>& y) ... ");
        return false;
    }
    return true;
}

template <typename T> 
void gtPlusWavelet2DOperator<T>::printInfo(std::ostream& os)
{
    using namespace std;

    os << "-------------- GTPlus ISMRMRD wavelet 2D operator --------------------" << endl;
    os << "Wavelet 2D operator for gtPlus ISMRMRD package" << endl;
    os << "----------------------------------------------------------------------" << endl;
}

}}
