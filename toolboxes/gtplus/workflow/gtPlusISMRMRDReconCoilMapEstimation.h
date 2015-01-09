/** \file   gtPlusISMRMRDReconCoilMapEstimation.h
    \brief  Implement coil map estimation methods.
    \author Hui Xue
*/

#pragma once

#include "GtPlusExport.h"
#include "gtPlusISMRMRDReconUtil.h"
#include "gtPlusSPIRIT.h"

namespace Gadgetron { namespace gtPlus {

// ================================================================================================== //

template <typename T> 
class gtPlusISMRMRDReconCoilMapEstimation
{
public:

    typedef typename realType<T>::Type value_type;

    gtPlusISMRMRDReconCoilMapEstimation();
    virtual ~gtPlusISMRMRDReconCoilMapEstimation();

    void printInfo(std::ostream& os);

    // compute dual coil map
    // data : ref kspace [RO E1 CHA]
    // coilMap : [RO E1 CHA 2]
    bool coilMap2DSPIRIT(const hoNDArray<T>& data, hoNDArray<T>& coilMap, hoNDArray<value_type>& eigD, size_t kRO, size_t kE1, value_type thres=0.01);
};

template <typename T> 
gtPlusISMRMRDReconCoilMapEstimation<T>::gtPlusISMRMRDReconCoilMapEstimation()
{
}

template <typename T> 
gtPlusISMRMRDReconCoilMapEstimation<T>::~gtPlusISMRMRDReconCoilMapEstimation()
{
}

template <typename T> 
bool gtPlusISMRMRDReconCoilMapEstimation<T>::coilMap2DSPIRIT(const hoNDArray<T>& data, hoNDArray<T>& coilMap, hoNDArray<value_type>& eigD, size_t kRO, size_t kE1, value_type thres)
{
    try
    {
        gtPlusSPIRIT<T> spirit;

        size_t oRO = 1;
        size_t oE1 = 1;

        size_t RO = coilMap.get_size(0);
        size_t E1 = coilMap.get_size(1);
        size_t CHA = data.get_size(2);

        ho3DArray<T> acsSrc(data.get_size(0), data.get_size(1), CHA, const_cast<T*>(data.begin()));
        ho3DArray<T> acsDst(data.get_size(0), data.get_size(1), CHA, const_cast<T*>(data.begin()));

        ho6DArray<T> ker(kRO, kE1, CHA, CHA, oRO, oE1);

        GADGET_CHECK_RETURN_FALSE(spirit.calib(acsSrc, acsDst, thres, kRO, kE1, oRO, oE1, ker));

        // std::string debugFolder_ = "D:/gtuser/mrprogs/gadgetron/toolboxes/gtplus/ut/result/";
        // Gadgetron::gtPlus::gtPlusIOAnalyze gt_exporter_;
        // if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(ker, debugFolder_+"ker"); }

        bool minusI = false;
        hoNDArray<T> kIm(RO, E1, CHA, CHA);
        GADGET_CHECK_RETURN_FALSE(spirit.imageDomainKernel(ker, kRO, kE1, oRO, oE1, RO, E1, kIm, minusI));
        T* pkIm = kIm.begin();

        // if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(kIm, debugFolder_+"kIm"); }

        coilMap.create(RO, E1, CHA, 2);
        eigD.create(RO, E1, 2);

        long long ro, e1, scha, dcha;

        #pragma omp parallel default(none) private(ro, e1, scha, dcha) shared(RO, E1, CHA, pkIm, coilMap, eigD)
        {
            hoMatrix<T> R(CHA, CHA), RC(CHA, CHA), RRT(CHA, CHA);
            Gadgetron::clear(RRT);

            hoMatrix<value_type> eigenValue;

            #pragma omp for 
            for ( e1=0; e1<E1; e1++ )
            {
                for ( ro=0; ro<RO; ro++ )
                {
                    const size_t offset = e1*RO + ro;

                    for ( dcha=0; dcha<CHA; dcha++ )
                    {
                        for ( scha=0; scha<CHA; scha++ )
                        {
                            // T v = kIm(ro, e1, scha, dcha);
                            T v = pkIm[dcha*RO*E1*CHA + scha*RO*E1 + offset];
                            if ( scha == dcha )
                            {
                                v -= 1;
                            }

                            R(scha, dcha) = v;
                        }
                    }

                    memcpy(RC.begin(), R.begin(), sizeof(T)*CHA*CHA);
                    Gadgetron::gemm(RRT, RC, false, R, true);

                    Gadgetron::heev(RRT, eigenValue);

                    for ( scha=0; scha<CHA; scha++ )
                    {
                        coilMap(ro, e1, scha, 0) = RRT(scha, 0);
                        coilMap(ro, e1, scha, 1) = RRT(scha, 1);
                        eigD(ro, e1, 0) = 1.0 - eigenValue(0, 0);
                        eigD(ro, e1, 1) = 1.0 - eigenValue(1, 0);
                    }
                }
            }
        }

        Gadgetron::conjugate(coilMap, coilMap);
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconCoilMapEstimation<T>::coilMap2DSPIRIT(...) ... ");
        return false;
    }

    return true;
}

}}
