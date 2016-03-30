
#include "gtPlusISMRMRDReconUtil.h"
#include "hoNDArray_elemwise.h"
#include <algorithm>

#define GT_IMAGING_GEOMETRY_DELTA 0.001

namespace Gadgetron { namespace gtPlus {

template <typename T> 
gtPlusISMRMRDReconUtil<T>::gtPlusISMRMRDReconUtil() {}

template <typename T> 
gtPlusISMRMRDReconUtil<T>::~gtPlusISMRMRDReconUtil() {}

template <typename T> 
void gtPlusISMRMRDReconUtil<T>::printInfo(std::ostream& os)
{
    using namespace std;
    os << "-------------- GTPlus ISMRMRD Recon Util -------------" << endl;
    os << "Implementation of recon utilities for ISMRMRD format" << endl;
    os << "------------------------------------------------------" << endl;
}

// ------------------------------------------------------------------------
// detect sampled region
// ------------------------------------------------------------------------

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::
detectSampledRegion2D(const hoNDArray<T>& data, size_t& startRO, size_t& endRO, size_t& startE1, size_t& endE1)
{
    try
    {
        size_t NDim = data.get_number_of_dimensions();

        hoNDArray<typename realType<T>::Type> mag(data.get_dimensions()), magSum, magSumE1, magSumRO;
        Gadgetron::abs(data, mag);

        if ( NDim > 2 )
        {
            size_t ii;
            std::vector<size_t> dim;
            for ( ii=0; ii<NDim-2; ii++ )
            {
                mag.get_dimensions(dim);

                // GADGET_CHECK_RETURN_FALSE(Gadgetron::sumOverLastDimension(mag, magSum));
                GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::sum_over_dimension(mag, magSum, mag.get_number_of_dimensions() - 1));

                std::vector<size_t> dimSum(dim.size()-1);
                memcpy(&dimSum[0], &dim[0], sizeof(size_t)*dimSum.size());
                magSum.reshape(dimSum);

                mag = magSum;
            }
        }

        size_t RO = mag.get_size(0);
        size_t E1 = mag.get_size(1);

        startRO = RO-1;
        endRO = 0;

        startE1 = E1-1;
        endE1 = 0;

        size_t ro, e1;

        GADGET_CATCH_THROW(Gadgetron::sum_over_dimension(mag, magSumE1, 1));
        GADGET_CATCH_THROW(Gadgetron::sum_over_dimension(mag, magSumRO, 0));

        for ( ro=0; ro<RO; ro++ )
        {
            if ( magSumE1(ro) > 0 )
            {
                if ( ro < startRO ) startRO = ro;
                if ( ro > endRO ) endRO = ro;
            }
        }

        for ( e1=0; e1<E1; e1++ )
        {
            if ( magSumRO(e1) > 0 )
            {
                if ( e1 < startE1 ) startE1 = e1;
                if ( e1 > endE1 ) endE1 = e1;
            }
        }

        if ( startRO > endRO )
        {
            startRO = 0;
            endRO = RO-1;
        }

        if ( startE1 > endE1 )
        {
            startE1 = 0;
            endE1 = E1-1;
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::detectSampledRegion2D(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::
detectSampledRegion3D(const hoNDArray<T>& data, size_t& startRO, size_t& endRO, size_t& startE1, size_t& endE1, size_t& startE2, size_t& endE2)
{
    try
    {
        size_t NDim = data.get_number_of_dimensions();

        hoNDArray<typename realType<T>::Type> mag(data.get_dimensions()), magSum, magSum2, magSumRO, magSumE1, magSumE2;
        Gadgetron::abs(data, mag);

        if ( NDim > 5 )
        {
            std::vector<size_t> dim;

            size_t ii;
            for ( ii=0; ii<NDim-5; ii++ )
            {
                mag.get_dimensions(dim);

                // GADGET_CHECK_RETURN_FALSE(Gadgetron::sumOverLastDimension(mag, magSum));
                GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::sum_over_dimension(mag, magSum, mag.get_number_of_dimensions() - 1));

                std::vector<size_t> dimSum(dim.size() - 1);
                memcpy(&dimSum[0], &dim[0], sizeof(size_t)*dimSum.size());
                magSum.reshape(dimSum);

                mag = magSum;
            }
        }

        size_t RO = mag.get_size(0);
        size_t E1 = mag.get_size(1);
        size_t E2 = mag.get_size(4);

        startRO = RO-1;
        endRO = 0;

        startE1 = E1-1;
        endE1 = 0;

        startE2 = E2-1;
        endE2 = 0;

        size_t ro, e1, e2;

        GADGET_CATCH_THROW(Gadgetron::sum_over_dimension(mag, magSum2, 4));
        GADGET_CATCH_THROW(Gadgetron::sum_over_dimension(magSum2, magSum, 3));
        GADGET_CATCH_THROW(Gadgetron::sum_over_dimension(magSum, magSum2, 2));

        GADGET_CATCH_THROW(Gadgetron::sum_over_dimension(magSum2, magSumE1, 1));
        GADGET_CATCH_THROW(Gadgetron::sum_over_dimension(magSum2, magSumRO, 0));

        GADGET_CATCH_THROW(Gadgetron::sum_over_dimension(mag, magSum2, 3));
        GADGET_CATCH_THROW(Gadgetron::sum_over_dimension(magSum2, magSum, 2));

        GADGET_CATCH_THROW(Gadgetron::sum_over_dimension(magSum, magSum2, 1));
        GADGET_CATCH_THROW(Gadgetron::sum_over_dimension(magSum2, magSumE2, 0));

        for ( ro=0; ro<RO; ro++ )
        {
            if ( magSumE1(ro) > 0 )
            {
                if ( ro < startRO ) startRO = ro;
                if ( ro > endRO ) endRO = ro;
            }
        }

        for ( e1=0; e1<E1; e1++ )
        {
            if ( magSumRO(e1) > 0 )
            {
                if ( e1 < startE1 ) startE1 = e1;
                if ( e1 > endE1 ) endE1 = e1;
            }
        }

        for ( e2=0; e2<E2; e2++ )
        {
            if ( magSumE2(e2) > 0 )
            {
                if ( e2 < startE2 ) startE2 = e2;
                if ( e2 > endE2 ) endE2 = e2;
            }
        }

        if ( startRO > endRO )
        {
            startRO = 0;
            endRO = RO-1;
        }

        if ( startE1 > endE1 )
        {
            startE1 = 0;
            endE1 = E1-1;
        }

        if ( startE2 > endE2 )
        {
            startE2 = 0;
            endE2 = E2-1;
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::detectSampledRegion3D(...) ... ");
        return false;
    }

    return true;
}

// ------------------------------------------------------------------------
// coil sensitivity
// ------------------------------------------------------------------------
template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::
averageKSpace4D(const hoNDArray<T>& data, hoNDArray<T>& ave)
{
    try
    {
        GADGET_CATCH_THROW(Gadgetron::sum_over_dimension(data, ave, 3));
        Gadgetron::scal( (typename realType<T>::Type)(1.0/data.get_size(3)), ave);
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::averageKSpace4D(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::
averageKSpace4D(const hoNDArray<T>& data, hoNDArray<T>& ave, std::vector<size_t>& sampledTimes)
{
    try
    {
        size_t NDim = data.get_number_of_dimensions();
        GADGET_CHECK_RETURN_FALSE(NDim>=2);

        if ( NDim < 4 )
        {
            ave = data;
            GADGET_CHECK_RETURN_FALSE(detectSampledTimesE1(data, sampledTimes));
            return true;
        }

        size_t RO = data.get_size(0);
        size_t E1 = data.get_size(1);
        size_t CHA = data.get_size(2);
        size_t N = data.get_size(3);

        hoNDArray<T> data4D(RO, E1, CHA, N, const_cast<T*>(data.begin()));
        GADGET_CHECK_RETURN_FALSE(detectSampledTimesE1(data4D, sampledTimes));
        GADGET_CATCH_THROW(Gadgetron::sum_over_dimension(data, ave, 3));

        boost::shared_ptr< std::vector<size_t> > dim = ave.get_dimensions();

        if ( dim->size() != NDim )
        {
            (*dim).insert((*dim).begin()+3, 1);
            ave.reshape(dim.get());
        }

        hoNDArray<T> sampledTimes2D(RO, E1);
        T* pTimes = sampledTimes2D.begin();
        size_t ro, e1;
        for ( e1=0; e1<E1; e1++ )
        {
            double t = (double)sampledTimes[e1];
            if ( t == 0 ) t = 1;

            for ( ro=0; ro<RO; ro++ )
            {
                pTimes[e1*RO+ro] = (value_type)(1.0/t);
            }
        }

        // GADGET_CHECK_RETURN_FALSE(multipleMultiply(sampledTimes2D, ave, ave));
        GADGET_CHECK_EXCEPTION_RETURN_FALSE(multiply(ave, sampledTimes2D, ave));
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::averageKSpace4D(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::
averageKSpace5D(const hoNDArray<T>& data, hoNDArray<T>& ave)
{
    try
    {
        GADGET_CATCH_THROW(Gadgetron::sum_over_dimension(data, ave, 4));
        Gadgetron::scal( (typename realType<T>::Type)(1.0/data.get_size(4)), ave);
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::averageKSpace5D(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::
averageKSpace5D(const hoNDArray<T>& data, hoNDArray<T>& ave, hoNDArray<size_t>& sampledTimes)
{
    try
    {
        size_t NDim = data.get_number_of_dimensions();
        GADGET_CHECK_RETURN_FALSE(NDim>=3);

        if ( NDim < 5 )
        {
            ave = data;
            GADGET_CHECK_RETURN_FALSE(detectSampledTimesE1E2(data, sampledTimes));
            return true;
        }

        size_t RO = data.get_size(0);
        size_t E1 = data.get_size(1);
        size_t E2 = data.get_size(2);
        size_t CHA = data.get_size(3);
        size_t N = data.get_size(4);

        hoNDArray<T> data5D(RO, E1, E2, CHA, N, const_cast<T*>(data.begin()));
        GADGET_CHECK_RETURN_FALSE(detectSampledTimesE1E2(data5D, sampledTimes));
        GADGET_CATCH_THROW(Gadgetron::sum_over_dimension(data, ave, 4));

        hoNDArray<T> sampledTimes3D(RO, E1, E2);
        T* pTimes = sampledTimes3D.begin();
        size_t ro, e1, e2;
        for ( e2=0; e2<E2; e2++ )
        {
            for ( e1=0; e1<E1; e1++ )
            {
                double t = (double)sampledTimes(e1+e2*E1);
                if ( t == 0 ) t = 1;

                for ( ro=0; ro<RO; ro++ )
                {
                    pTimes[e2*RO*E1+e1*RO+ro] = (value_type)(1.0/t);
                }
            }
        }

        // GADGET_CHECK_RETURN_FALSE(multipleMultiply(sampledTimes3D, ave, ave));
        GADGET_CHECK_EXCEPTION_RETURN_FALSE(multiply(ave, sampledTimes3D, ave));
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::averageKSpace5D(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::
detectSampledTimesE1(const hoNDArray<T>& data4D, std::vector<size_t>& sampledTimes)
{
    try
    {
        size_t NDim = data4D.get_number_of_dimensions();
        GADGET_CHECK_RETURN_FALSE(NDim>=2);

        size_t RO = data4D.get_size(0);
        size_t E1 = data4D.get_size(1);
        size_t CHA = data4D.get_size(2);
        size_t N = data4D.get_size(3);

        hoNDArray<typename realType<T>::Type> mag(data4D.get_dimensions());
        Gadgetron::abs(data4D, mag);

        hoNDArray<typename realType<T>::Type> mag3D(RO, E1, 1, N);
        GADGET_CATCH_THROW(Gadgetron::sum_over_dimension(mag, mag3D, 2));

        hoNDArray<typename realType<T>::Type> mag2D(1, E1, 1, N);
        GADGET_CATCH_THROW(Gadgetron::sum_over_dimension(mag3D, mag2D, 0));
        typename realType<T>::Type* pMag2D = mag2D.begin();

        sampledTimes.resize(E1, 0);

        size_t e1, n;
        for ( e1=0; e1<E1; e1++ )
        {
            for ( n=0; n<N; n++ )
            {
                if ( pMag2D[e1+n*E1] > 0 )
                {
                    sampledTimes[e1]++;
                }
            }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::detectSampledTimesE1(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::
detectSampledRegionE1(const hoNDArray<T>& data, size_t& startE1, size_t& endE1)
{
    try
    {
        std::vector<size_t> sampledTimes;
        GADGET_CHECK_RETURN_FALSE(detectSampledTimesE1(data, sampledTimes));

        size_t E1 = sampledTimes.size();

        startE1 = E1-1;
        endE1 = 0;

        for ( size_t e1=0; e1<E1; e1++ )
        {
            if ( sampledTimes[e1] > 0 )
            {
                if ( e1 > endE1 ) endE1 = e1;
                if ( e1 < startE1 ) startE1 = e1;
            }
        }

        if ( endE1 < startE1 )
        {
            startE1 = 0;
            endE1 = E1-1;
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::detectSampledRegionE1(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::
detectSampledTimesE1E2(const hoNDArray<T>& data5D, hoNDArray<size_t>& sampledTimes)
{
    try
    {
        size_t NDim = data5D.get_number_of_dimensions();
        GADGET_CHECK_RETURN_FALSE(NDim>=3);

        size_t RO = data5D.get_size(0);
        size_t E1 = data5D.get_size(1);
        size_t E2 = data5D.get_size(2);
        size_t CHA = data5D.get_size(3);
        size_t N = data5D.get_size(4);

        hoNDArray<typename realType<T>::Type> mag(RO, E1, E2);

        hoNDArray<T> dataFirstChannel(RO, E1, E2, const_cast<T*>(data5D.begin()));
        Gadgetron::abs(dataFirstChannel, mag);

        hoNDArray<typename realType<T>::Type> mag3D(1, E1, E2);
        GADGET_CATCH_THROW(Gadgetron::sum_over_dimension(mag, mag3D, 0));

        typename realType<T>::Type* pMag3D = mag3D.begin();

        sampledTimes.create(E1, E2);
        Gadgetron::clear(sampledTimes);
        size_t* pTimes = sampledTimes.get_data_ptr();

        size_t e1, e2, n;
        for ( e2=0; e2<E2; e2++ )
        {
            for ( e1=0; e1<E1; e1++ )
            {
                for ( n=0; n<N; n++ )
                {
                    if ( pMag3D[e1+e2*E1+n*E1*E2] > 0 )
                    {
                        pTimes[e1+e2*E1]++;
                    }
                }
            }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::detectSampledTimesE1E2(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::
detectSampledRegionE1E2(const hoNDArray<T>& data, size_t& startE1, size_t& endE1, size_t& startE2, size_t& endE2)
{
    try
    {
        hoNDArray<size_t> sampledTimes;
        GADGET_CHECK_RETURN_FALSE(detectSampledTimesE1E2(data, sampledTimes));

        size_t E1 = sampledTimes.get_size(0);
        size_t E2 = sampledTimes.get_size(1);

        startE1 = E1-1;
        endE1 = 0;

        startE2 = E2-1;
        endE2 = 0;

        size_t e1, e2;
        for ( e2=0; e2<E2; e2++ )
        {
            for ( e1=0; e1<E1; e1++ )
            {
                if ( sampledTimes(e1+e2*E1) > 0 )
                {
                    if ( e1 > endE1 ) endE1 = e1;
                    if ( e1 < startE1 ) startE1 = e1;

                    if ( e2 > endE2 ) endE2 = e2;
                    if ( e2 < startE2 ) startE2 = e2;
                }
            }
        }

        if ( endE1 < startE1 )
        {
            startE1 = 0;
            endE1 = E1-1;
        }

        if ( endE2 < startE2 )
        {
            startE2 = 0;
            endE2 = E2-1;
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::detectSampledRegionE1E2(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::
copyAlongE1(const hoNDArray<T>& src, hoNDArray<T>& dst, size_t startE1, size_t endE1)
{
    try
    {
        size_t NDim = src.get_number_of_dimensions();
        GADGET_CHECK_RETURN_FALSE(NDim>=2);

        size_t RO = dst.get_size(0);
        size_t E1 = dst.get_size(1);

        size_t RO_src = src.get_size(0);
        size_t E1_src = src.get_size(1);

        GADGET_CHECK_RETURN_FALSE(RO==RO_src);
        GADGET_CHECK_RETURN_FALSE(E1==E1_src);
        GADGET_CHECK_RETURN_FALSE(src.get_number_of_elements()==dst.get_number_of_elements());

        if ( (startE1>=E1) || (endE1>=E1) || (startE1>endE1) )
        {
            dst = src;
            GWARN_STREAM("copyAlongE1(...) : (startE1>=E1) || (endE1>=E1) || (startE1>endE1) ... ");
            return true;
        }

        size_t N = dst.get_number_of_elements()/(RO*E1);

        size_t n, e1;
        for ( n=0; n<N; n++ )
        {
            for ( e1=startE1; e1<=endE1; e1++ )
            {
                memcpy(dst.begin()+n*RO*E1+e1*RO, src.begin()+n*RO*E1+e1*RO, sizeof(T)*RO);
            }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::copyAlongE1(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::
copyAlongROE1(const hoNDArray<T>& src, hoNDArray<T>& dst, size_t startRO, size_t endRO, size_t startE1, size_t endE1)
{
    try
    {
        size_t NDim = src.get_number_of_dimensions();
        GADGET_CHECK_RETURN_FALSE(NDim>=2);

        size_t RO = dst.get_size(0);
        size_t E1 = dst.get_size(1);

        size_t RO_src = src.get_size(0);
        size_t E1_src = src.get_size(1);

        GADGET_CHECK_RETURN_FALSE(RO==RO_src);
        GADGET_CHECK_RETURN_FALSE(E1==E1_src);
        GADGET_CHECK_RETURN_FALSE(src.get_number_of_elements()==dst.get_number_of_elements());

        if ( (startRO>=RO) || (endRO>=RO) || (startRO>endRO) )
        {
            dst = src;
            GWARN_STREAM("copyAlongROE1(...) : (startRO>=RO) || (endRO>=RO) || (startRO>endRO) ... ");
            return true;
        }

        if ( (startE1>=E1) || (endE1>=E1) || (startE1>endE1) )
        {
            dst = src;
            GWARN_STREAM("copyAlongROE1(...) : (startE1>=E1) || (endE1>=E1) || (startE1>endE1) ... ");
            return true;
        }

        size_t N = dst.get_number_of_elements()/(RO*E1);
        const T* pSrc = src.begin();
        T* pDst = dst.begin();

        long long n;

        #pragma omp parallel for default(none) private(n) shared(N, pSrc, pDst, RO, E1, startRO, endRO, startE1, endE1)
        for ( n=0; n<(long long)N; n++ )
        {
            for ( size_t e1=startE1; e1<=endE1; e1++ )
            {
                size_t offset = n*RO*E1+e1*RO+startRO;
                memcpy(pDst+offset, pSrc+offset, sizeof(T)*(endRO-startRO+1));
            }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::copyAlongROE1(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::
copyAlongROE1E2(const hoNDArray<T>& src, hoNDArray<T>& dst, size_t startRO, size_t endRO, size_t startE1, size_t endE1, size_t startE2, size_t endE2)
{
    try
    {
        size_t NDim = src.get_number_of_dimensions();
        GADGET_CHECK_RETURN_FALSE(NDim>=2);

        size_t RO = dst.get_size(0);
        size_t E1 = dst.get_size(1);
        size_t E2 = dst.get_size(2);

        size_t RO_src = src.get_size(0);
        size_t E1_src = src.get_size(1);
        size_t E2_src = src.get_size(2);

        GADGET_CHECK_RETURN_FALSE(RO==RO_src);
        GADGET_CHECK_RETURN_FALSE(E1==E1_src);
        GADGET_CHECK_RETURN_FALSE(E2==E2_src);
        GADGET_CHECK_RETURN_FALSE(src.get_number_of_elements()==dst.get_number_of_elements());

        if ( (startRO>=RO) || (endRO>=RO) || (startRO>endRO) )
        {
            dst = src;
            GWARN_STREAM("copyAlongROE1E2(...) : (startRO>=RO) || (endRO>=RO) || (startRO>endRO) ... ");
            return true;
        }

        if ( (startE1>=E1) || (endE1>=E1) || (startE1>endE1) )
        {
            dst = src;
            GWARN_STREAM("copyAlongROE1E2(...) : (startE1>=E1) || (endE1>=E1) || (startE1>endE1) ... ");
            return true;
        }

        if ( (startE2>=E2) || (endE2>=E2) || (startE2>endE2) )
        {
            dst = src;
            GWARN_STREAM("copyAlongROE1E2(...) : (startE2>=E2) || (endE2>=E2) || (startE2>endE2) ... ");
            return true;
        }

        size_t N = dst.get_number_of_elements()/(RO*E1*E2);
        const T* pSrc = src.begin();
        T* pDst = dst.begin();

        long long n;

        #pragma omp parallel for default(none) private(n) shared(N, pSrc, pDst, RO, E1, E2, startRO, endRO, startE1, endE1, startE2, endE2)
        for ( n=0; n<(long long)N; n++ )
        {
            for ( size_t e2=startE2; e2<=endE2; e2++ )
            {
                for ( size_t e1=startE1; e1<=endE1; e1++ )
                {
                    size_t offset = n*RO*E1*E2+e2*E1*RO+e1*RO+startRO;
                    memcpy(pDst+offset, pSrc+offset, sizeof(T)*(endRO-startRO+1));
                }
            }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::copyAlongROE1E2(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::
copyAlongROE1TransitionBand(const hoNDArray<T>& src, hoNDArray<T>& dst, size_t startRO, size_t endRO, 
        size_t startE1, size_t endE1, size_t transBandRO, size_t transBandE1)
{
    try
    {
        size_t NDim = src.get_number_of_dimensions();
        GADGET_CHECK_RETURN_FALSE(NDim>=2);

        size_t RO = dst.get_size(0);
        size_t E1 = dst.get_size(1);

        size_t RO_src = src.get_size(0);
        size_t E1_src = src.get_size(1);

        GADGET_CHECK_RETURN_FALSE(RO==RO_src);
        GADGET_CHECK_RETURN_FALSE(E1==E1_src);
        GADGET_CHECK_RETURN_FALSE(src.get_number_of_elements()==dst.get_number_of_elements());

        if ( (startRO>=RO) || (endRO>=RO) || (startRO>endRO) )
        {
            dst = src;
            GWARN_STREAM("copyAlongROE1TransitionBand(...) : (startRO>=RO) || (endRO>=RO) || (startRO>endRO) ... ");
            return true;
        }

        if ( (startE1>=E1) || (endE1>=E1) || (startE1>endE1) )
        {
            dst = src;
            GWARN_STREAM("copyAlongROE1TransitionBand(...) : (startE1>=E1) || (endE1>=E1) || (startE1>endE1) ... ");
            return true;
        }

        while ( transBandRO>1 && startRO+transBandRO > RO/2 )
        {
             transBandRO--;
        }

        while ( transBandRO>1 && endRO-transBandRO < RO/2 )
        {
             transBandRO--;
        }

        while ( transBandE1>1 && startE1+transBandE1 > E1/2 )
        {
             transBandE1--;
        }

        while ( transBandE1>1 && endE1-transBandE1 < E1/2 )
        {
             transBandE1--;
        }

        ISMRMRDKSPACEFILTER filterType = ISMRMRD_FILTER_TAPERED_HANNING;
        bool densityComp = false;

        hoNDArray<T> filter_src_RO, filter_src_E1;

        if ( startRO==0 && endRO==RO-1 )
        {
            Gadgetron::generate_asymmetric_filter(RO, startRO, endRO, filter_src_RO, ISMRMRD_FILTER_NONE, transBandRO, densityComp);
        }
        else
        {
            Gadgetron::generate_asymmetric_filter(RO, startRO, endRO, filter_src_RO, ISMRMRD_FILTER_TAPERED_HANNING, transBandRO, densityComp);
        }

        if ( startE1==0 && endE1==E1-1 )
        {
            Gadgetron::generate_asymmetric_filter(E1, startE1, endE1, filter_src_E1, ISMRMRD_FILTER_NONE, transBandE1, densityComp);
        }
        else
        {
            Gadgetron::generate_asymmetric_filter(E1, startE1, endE1, filter_src_E1, ISMRMRD_FILTER_TAPERED_HANNING, transBandE1, densityComp);
        }

        // in this way, the SNR unit scale property is perserved
        T midValue = filter_src_RO(RO/2);
        T scalFactor = T(1.0)/midValue;
        Gadgetron::scal(scalFactor, filter_src_RO);

        midValue = filter_src_E1(E1/2);
        scalFactor = T(1.0)/midValue;
        Gadgetron::scal(scalFactor, filter_src_E1);

        hoNDArray<T> filter_dst_RO(RO), filter_dst_E1(E1);

        size_t ii;
        for ( ii=0; ii<RO; ii++ )
        {
            filter_dst_RO(ii) = T(1.0) - filter_src_RO(ii);
        }

        for ( ii=0; ii<E1; ii++ )
        {
            filter_dst_E1(ii) = T(1.0) - filter_src_E1(ii);
        }

        hoNDArray<T> srcFiltered(src), dstFiltered(dst);
        if ( startRO==0 && endRO==RO-1 )
        {
            GADGET_CHECK_EXCEPTION_RETURN_FALSE(apply_kspace_filter_E1(src, filter_src_E1, srcFiltered));
            GADGET_CHECK_EXCEPTION_RETURN_FALSE(apply_kspace_filter_E1(dst, filter_dst_E1, dstFiltered));
        }
        else if ( startE1==0 && endE1==E1-1 )
        {
            GADGET_CHECK_EXCEPTION_RETURN_FALSE(apply_kspace_filter_RO(src, filter_src_RO, srcFiltered));
            GADGET_CHECK_EXCEPTION_RETURN_FALSE(apply_kspace_filter_RO(dst, filter_dst_RO, dstFiltered));
        }
        else
        {
            GADGET_CHECK_EXCEPTION_RETURN_FALSE(apply_kspace_filter_ROE1(src, filter_src_RO, filter_src_E1, srcFiltered));

            hoNDArray<T> fxy;
            GADGET_CHECK_EXCEPTION_RETURN_FALSE(compute_2d_filter(filter_src_RO, filter_src_E1, fxy));

            size_t Nxy = RO*E1;
            for ( ii=0; ii<Nxy; ii++ )
            {
                fxy(ii) = T(1.0) - fxy(ii);
            }

            GADGET_CHECK_EXCEPTION_RETURN_FALSE(apply_kspace_filter_ROE1(dst, fxy, dstFiltered));
        }

        Gadgetron::add(srcFiltered, dstFiltered, dst);
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::copyAlongROE1TransitionBand(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::
copyAlongROE1E2TransitionBand(const hoNDArray<T>& src, hoNDArray<T>& dst, size_t startRO, size_t endRO, 
                        size_t startE1, size_t endE1, size_t startE2, size_t endE2, 
                        size_t transBandRO, size_t transBandE1, size_t transBandE2)
{
    try
    {
        size_t NDim = src.get_number_of_dimensions();
        GADGET_CHECK_RETURN_FALSE(NDim>=3);

        size_t RO = dst.get_size(0);
        size_t E1 = dst.get_size(1);
        size_t E2 = dst.get_size(2);

        size_t RO_src = src.get_size(0);
        size_t E1_src = src.get_size(1);
        size_t E2_src = src.get_size(2);

        GADGET_CHECK_RETURN_FALSE(RO==RO_src);
        GADGET_CHECK_RETURN_FALSE(E1==E1_src);
        GADGET_CHECK_RETURN_FALSE(E2==E2_src);
        GADGET_CHECK_RETURN_FALSE(src.get_number_of_elements()==dst.get_number_of_elements());

        if ( (startRO>=RO) || (endRO>=RO) || (startRO>endRO) )
        {
            dst = src;
            GWARN_STREAM("copyAlongROE1TransitionBand(...) : (startRO>=RO) || (endRO>=RO) || (startRO>endRO) ... ");
            return true;
        }

        if ( (startE1>=E1) || (endE1>=E1) || (startE1>endE1) )
        {
            dst = src;
            GWARN_STREAM("copyAlongROE1TransitionBand(...) : (startE1>=E1) || (endE1>=E1) || (startE1>endE1) ... ");
            return true;
        }

        if ( (startE2>=E2) || (endE2>=E2) || (startE2>endE2) )
        {
            dst = src;
            GWARN_STREAM("copyAlongROE1E2TransitionBand(...) : (startE2>=E2) || (endE2>=E2) || (startE2>endE2) ... ");
            return true;
        }

        while ( transBandRO>1 && startRO+transBandRO > RO/2 )
        {
             transBandRO--;
        }

        while ( transBandRO>1 && endRO-transBandRO < RO/2 )
        {
             transBandRO--;
        }

        while ( transBandE1>1 && startE1+transBandE1 > E1/2 )
        {
             transBandE1--;
        }

        while ( transBandE1>1 && endE1-transBandE1 < E1/2 )
        {
             transBandE1--;
        }

        while ( transBandE2>1 && startE2+transBandE2 > E2/2 )
        {
             transBandE2--;
        }

        while ( transBandE2>1 && endE2-transBandE2 < E2/2 )
        {
             transBandE2--;
        }

        ISMRMRDKSPACEFILTER filterType = ISMRMRD_FILTER_TAPERED_HANNING;
        bool densityComp = false;

        hoNDArray<T> filter_src_RO, filter_src_E1, filter_src_E2;

        if ( startRO==0 && endRO==RO-1 )
        {
            Gadgetron::generate_asymmetric_filter(RO, startRO, endRO, filter_src_RO, ISMRMRD_FILTER_NONE, transBandRO, densityComp);
        }
        else
        {
            Gadgetron::generate_asymmetric_filter(RO, startRO, endRO, filter_src_RO, ISMRMRD_FILTER_TAPERED_HANNING, transBandRO, densityComp);
        }

        if ( startE1==0 && endE1==E1-1 )
        {
            Gadgetron::generate_asymmetric_filter(E1, startE1, endE1, filter_src_E1, ISMRMRD_FILTER_NONE, transBandE1, densityComp);
        }
        else
        {
            Gadgetron::generate_asymmetric_filter(E1, startE1, endE1, filter_src_E1, ISMRMRD_FILTER_TAPERED_HANNING, transBandE1, densityComp);
        }

        if ( startE2==0 && endE2==E2-1 )
        {
            Gadgetron::generate_asymmetric_filter(E2, startE2, endE2, filter_src_E2, ISMRMRD_FILTER_NONE, transBandE2, densityComp);
        }
        else
        {
            Gadgetron::generate_asymmetric_filter(E2, startE2, endE2, filter_src_E2, ISMRMRD_FILTER_TAPERED_HANNING, transBandE2, densityComp);
        }

        // in this way, the SNR unit scale property is perserved
        T midValue = filter_src_RO(RO/2);
        T scalFactor = T(1.0)/midValue;
        Gadgetron::scal(scalFactor, filter_src_RO);

        midValue = filter_src_E1(E1/2);
        scalFactor = T(1.0)/midValue;
        Gadgetron::scal(scalFactor, filter_src_E1);

        midValue = filter_src_E2(E2/2);
        scalFactor = T(1.0)/midValue;
        Gadgetron::scal(scalFactor, filter_src_E2);

        hoNDArray<T> filter_dst_RO(RO), filter_dst_E1(E1), filter_dst_E2(E2);

        size_t ii;
        for ( ii=0; ii<RO; ii++ )
        {
            filter_dst_RO(ii) = T(1.0) - filter_src_RO(ii);
        }

        for ( ii=0; ii<E1; ii++ )
        {
            filter_dst_E1(ii) = T(1.0) - filter_src_E1(ii);
        }

        for ( ii=0; ii<E2; ii++ )
        {
            filter_dst_E2(ii) = T(1.0) - filter_src_E2(ii);
        }

        hoNDArray<T> srcFiltered(src), dstFiltered(dst);
        if ( endRO<=RO-1 && startE1==0 && endE1==E1-1 && startE2==0 && endE1==E2-1 )
        {
            GADGET_CHECK_EXCEPTION_RETURN_FALSE(apply_kspace_filter_RO(src, filter_src_E1, srcFiltered));
            GADGET_CHECK_EXCEPTION_RETURN_FALSE(apply_kspace_filter_RO(dst, filter_dst_E1, dstFiltered));
        }
        else if ( startRO==0 && endRO==RO-1 && endE1<=E1-1 && startE2==0 && endE1==E2-1 )
        {
            GADGET_CHECK_EXCEPTION_RETURN_FALSE(apply_kspace_filter_E1(src, filter_src_RO, srcFiltered));
            GADGET_CHECK_EXCEPTION_RETURN_FALSE(apply_kspace_filter_E1(dst, filter_dst_RO, dstFiltered));
        }
        else if ( startRO==0 && endRO==RO-1 && startE1==0 && endE1==E1-1 && endE1<=E2-1 )
        {
            GADGET_CHECK_EXCEPTION_RETURN_FALSE(apply_kspace_filter_E2(src, filter_src_RO, srcFiltered));
            GADGET_CHECK_EXCEPTION_RETURN_FALSE(apply_kspace_filter_E2(dst, filter_dst_RO, dstFiltered));
        }
        else
        {
            GADGET_CHECK_EXCEPTION_RETURN_FALSE(apply_kspace_filter_ROE1E2(src, filter_src_RO, filter_src_E1, filter_src_E2, srcFiltered));

            hoNDArray<T> fxyz;
            GADGET_CHECK_EXCEPTION_RETURN_FALSE(compute_3d_filter(filter_src_RO, filter_src_E1, filter_src_E2, fxyz));

            size_t Nxyz = RO*E1*E2;
            for ( ii=0; ii<Nxyz; ii++ )
            {
                fxyz(ii) = T(1.0) - fxyz(ii);
            }

            GADGET_CHECK_EXCEPTION_RETURN_FALSE(apply_kspace_filter_ROE1E2(dst, fxyz, dstFiltered));
        }

        Gadgetron::add(srcFiltered, dstFiltered, dst);
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::copyAlongROE1E2TransitionBand(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
std::string gtPlusISMRMRDReconUtil<T>::getISMRMRDDimName(const ISMRMRDDIM& dim)
{
    std::ostringstream os;
    switch (dim)
    {
        case DIM_ReadOut:
            os << "DIM_ReadOut";
        break;

        case DIM_Encoding1:
            os << "DIM_Encoding1";
        break;

        case DIM_Channel:
            os << "DIM_Channel";
        break;

        case DIM_Slice:
            os << "DIM_Slice";
        break;

        case DIM_Encoding2:
            os << "DIM_Encoding2";
        break;

        case DIM_Contrast:
            os << "DIM_Contrast";
        break;

        case DIM_Phase:
            os << "DIM_Phase";
        break;

        case DIM_Repetition:
            os << "DIM_Repetition";
        break;

        case DIM_Set:
            os << "DIM_Set";
        break;

        case DIM_Segment:
            os << "DIM_Segment";
        break;

        case DIM_Average:
            os << "DIM_Average";
        break;

        case DIM_other1:
            os << "DIM_other1";
        break;

        case DIM_other2:
            os << "DIM_other2";
        break;

        case DIM_other3:
            os << "DIM_other3";
        break;

        default:
            os << "DIM_NONE";
    }

    std::string dimStr(os.str());
    return dimStr;
}

template <typename T> 
ISMRMRDDIM gtPlusISMRMRDReconUtil<T>::getISMRMRDDimFromName(const std::string& name)
{
    if ( name == "DIM_ReadOut" ) return DIM_ReadOut;
    if ( name == "DIM_Encoding1" ) return DIM_Encoding1;
    if ( name == "DIM_Channel" ) return DIM_Channel;
    if ( name == "DIM_Slice" ) return DIM_Slice;
    if ( name == "DIM_Encoding2" ) return DIM_Encoding2;
    if ( name == "DIM_Contrast" ) return DIM_Contrast;
    if ( name == "DIM_Phase" ) return DIM_Phase;
    if ( name == "DIM_Repetition" ) return DIM_Repetition;
    if ( name == "DIM_Set" ) return DIM_Set;
    if ( name == "DIM_Segment" ) return DIM_Segment;
    if ( name == "DIM_Average" ) return DIM_Average;
    if ( name == "DIM_other1" ) return DIM_other1;
    if ( name == "DIM_other2" ) return DIM_other2;
    if ( name == "DIM_other3" ) return DIM_other3;

    return DIM_NONE;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::getISMRMRDDimIndex(const ISMRMRDDIM& dim, long long& ind)
{
    switch (dim)
    {
        case Gadgetron::DIM_ReadOut:
            ind = 0;
        break;

        case Gadgetron::DIM_Encoding1:
            ind = 1;
        break;

        case Gadgetron::DIM_Channel:
            ind = 2;
        break;

        case Gadgetron::DIM_Slice:
            ind = 3;
        break;

        case Gadgetron::DIM_Encoding2:
            ind = 4;
        break;

        case Gadgetron::DIM_Contrast:
            ind = 5;
        break;

        case Gadgetron::DIM_Phase:
            ind = 6;
        break;

        case Gadgetron::DIM_Repetition:
            ind = 7;
        break;

        case Gadgetron::DIM_Set:
            ind = 8;
        break;

        case Gadgetron::DIM_Segment:
            ind = 9;
        break;

        case Gadgetron::DIM_Average:
            ind = 10;
        break;

        default:
            ind = -1;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::findDimIndex(const std::vector<DimensionRecordType>& dimStartingIndexes, ISMRMRDDIM dim, size_t ind)
{
    size_t N = dimStartingIndexes.size();

    size_t n;
    for ( n=0; n<N; n++ )
    {
        if ( dimStartingIndexes[n].first == dim )
        {
            ind = dimStartingIndexes[n].second;
            return true;
        }
    }

    return false;
}

template <typename T> 
ISMRMRDALGO gtPlusISMRMRDReconUtil<T>::getISMRMRDReconAlgoFromName(const std::string& name)
{
    if ( name == "ISMRMRD_GRAPPA" ) return ISMRMRD_GRAPPA;
    if ( name == "ISMRMRD_SENSE" ) return ISMRMRD_SENSE;
    if ( name == "ISMRMRD_SPIRIT" ) return ISMRMRD_SPIRIT;
    if ( name == "ISMRMRD_SOFTSENSE" ) return ISMRMRD_SOFTSENSE;
    if ( name == "ISMRMRD_L1SOFTSENSE" ) return ISMRMRD_L1SOFTSENSE;
    if ( name == "ISMRMRD_2DTBINNING" ) return ISMRMRD_2DTBINNING;
    if ( name == "ISMRMRD_2DTBINNING_FLOW" ) return ISMRMRD_2DTBINNING_FLOW;
    if ( name == "ISMRMRD_L1SPIRIT_SLEP" ) return ISMRMRD_L1SPIRIT_SLEP;
    if ( name == "ISMRMRD_L1SPIRIT_SLEP_MOTION_COMP" ) return ISMRMRD_L1SPIRIT_SLEP_MOTION_COMP;

    return ISMRMRD_NONE;
}

template <typename T> 
ISMRMRDCOILMAPALGO gtPlusISMRMRDReconUtil<T>::getISMRMRDCoilMapAlgoFromName(const std::string& name)
{
    if ( name == "ISMRMRD_SOUHEIL" ) return ISMRMRD_SOUHEIL;
    if ( name == "ISMRMRD_SOUHEIL_ITER" ) return ISMRMRD_SOUHEIL_ITER;

    return ISMRMRD_SOUHEIL;
}

template <typename T> 
ISMRMRDPFALGO gtPlusISMRMRDReconUtil<T>::getISMRMRDPartialFourierReconAlgoFromName(const std::string& name)
{
    if ( name == "ISMRMRD_PF_HOMODYNE" ) return ISMRMRD_PF_HOMODYNE;
    if ( name == "ISMRMRD_PF_FENGHUANG" ) return ISMRMRD_PF_FENGHUANG;
    if ( name == "ISMRMRD_PF_POCS" ) return ISMRMRD_PF_POCS;
    if ( name == "ISMRMRD_PF_ZEROFILLING_FILTER" ) return ISMRMRD_PF_ZEROFILLING_FILTER;
    if ( name == "ISMRMRD_PF_ZEROFILLING" ) return ISMRMRD_PF_ZEROFILLING;

    return ISMRMRD_PF_NONE;
}

template <typename T> 
std::string gtPlusISMRMRDReconUtil<T>::getNameFromISMRMRDPartialFourierReconAlgo(ISMRMRDPFALGO algo)
{
    if ( algo == ISMRMRD_PF_HOMODYNE ) return std::string("ISMRMRD_PF_HOMODYNE");
    if ( algo == ISMRMRD_PF_FENGHUANG ) return std::string("ISMRMRD_PF_FENGHUANG");
    if ( algo == ISMRMRD_PF_ZEROFILLING_FILTER ) return std::string("ISMRMRD_PF_ZEROFILLING_FILTER");
    if ( algo == ISMRMRD_PF_POCS ) return std::string("ISMRMRD_PF_POCS");
    if ( algo == ISMRMRD_PF_ZEROFILLING ) return std::string("ISMRMRD_PF_ZEROFILLING");

    return std::string("ISMRMRD_PF_NONE");
}

template <typename T> 
ISMRMRDINTERPRETROGATING gtPlusISMRMRDReconUtil<T>::getISMRMRDRetroGatingInterpFromName(const std::string& name)
{
    if ( name == "ISMRMRD_INTERP_RETRO_GATING_LINEAR" ) return ISMRMRD_INTERP_RETRO_GATING_LINEAR;
    if ( name == "ISMRMRD_INTERP_RETRO_GATING_CUBIC" ) return ISMRMRD_INTERP_RETRO_GATING_CUBIC;
    if ( name == "ISMRMRD_INTERP_RETRO_GATING_BSPLINE" ) return ISMRMRD_INTERP_RETRO_GATING_BSPLINE;

    return ISMRMRD_INTERP_RETRO_GATING_LINEAR;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::
extractSubArrayForDim(const hoNDArray<T>& x, hoNDArray<T>& r, ISMRMRDDIM& dim, size_t value, bool lessEqual)
{
    try
    {
        boost::shared_ptr< std::vector<size_t> > dimX = x.get_dimensions();

        long long dimInd;
        GADGET_CHECK_RETURN_FALSE(getISMRMRDDimIndex(dim, dimInd));

        GADGET_CHECK_RETURN_FALSE(value<(*dimX)[dimInd]);

        std::vector<size_t> crop_offset(11, 0);
        crop_offset[0] = 0;
        crop_offset[1] = 0;
        crop_offset[2] = 0;
        crop_offset[3] = 0;
        crop_offset[4] = 0;
        crop_offset[5] = 0;
        crop_offset[6] = 0;
        crop_offset[7] = 0;
        crop_offset[8] = 0;
        crop_offset[9] = 0;
        crop_offset[10] = 0;

        std::vector<size_t> crop_size(11, 0);
        crop_size[0] = (*dimX)[0];
        crop_size[1] = (*dimX)[1];
        crop_size[2] = (*dimX)[2];
        crop_size[3] = (*dimX)[3];
        crop_size[4] = (*dimX)[4];
        crop_size[5] = (*dimX)[5];
        crop_size[6] = (*dimX)[6];
        crop_size[7] = (*dimX)[7];
        crop_size[8] = (*dimX)[8];
        crop_size[9] = (*dimX)[9];
        crop_size[10] = (*dimX)[10];

        if ( lessEqual )
        {
            crop_size[dimInd] = value+1;
        }
        else
        {
            crop_offset[dimInd] = value;
            crop_size[dimInd] = 1;
        }

        GADGET_CHECK_RETURN_FALSE(Gadgetron::cropUpTo11DArray(x, r, crop_offset, crop_size));
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::extractSubArrayForDim(dim, value) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::
extractSubArrayForDim(const hoNDArray<T>& x, hoNDArray<T>& r, ISMRMRDDIM& dim1, size_t value1, ISMRMRDDIM& dim2, size_t value2, bool lessEqual)
{
    try
    {
        boost::shared_ptr< std::vector<size_t> > dimX = x.get_dimensions();

        long long dimInd1, dimInd2;
        GADGET_CHECK_RETURN_FALSE(getISMRMRDDimIndex(dim1, dimInd1));
        GADGET_CHECK_RETURN_FALSE(getISMRMRDDimIndex(dim2, dimInd2));

        GADGET_CHECK_RETURN_FALSE(value1<(*dimX)[dimInd1]);
        GADGET_CHECK_RETURN_FALSE(value2<(*dimX)[dimInd2]);

        std::vector<size_t> crop_offset(11, 0);
        crop_offset[0] = 0;
        crop_offset[1] = 0;
        crop_offset[2] = 0;
        crop_offset[3] = 0;
        crop_offset[4] = 0;
        crop_offset[5] = 0;
        crop_offset[6] = 0;
        crop_offset[7] = 0;
        crop_offset[8] = 0;
        crop_offset[9] = 0;
        crop_offset[10] = 0;

        std::vector<size_t> crop_size(11, 0);
        crop_size[0] = (*dimX)[0];
        crop_size[1] = (*dimX)[1];
        crop_size[2] = (*dimX)[2];
        crop_size[3] = (*dimX)[3];
        crop_size[4] = (*dimX)[4];
        crop_size[5] = (*dimX)[5];
        crop_size[6] = (*dimX)[6];
        crop_size[7] = (*dimX)[7];
        crop_size[8] = (*dimX)[8];
        crop_size[9] = (*dimX)[9];
        crop_size[10] = (*dimX)[10];

        if ( lessEqual )
        {
            crop_size[dimInd1] = value1+1;
            crop_size[dimInd2] = value2+1;
        }
        else
        {
            crop_offset[dimInd1] = value1;
            crop_size[dimInd1] = 1;

            crop_offset[dimInd2] = value2;
            crop_size[dimInd2] = 1;
        }

        GADGET_CHECK_RETURN_FALSE(Gadgetron::cropUpTo11DArray(x, r, crop_offset, crop_size));
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::extractSubArrayForDim(dim1, value1, dim2, value2) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::
extractSubArrayForDim1LessEqualDim2Equal(const hoNDArray<T>& x, hoNDArray<T>& r, ISMRMRDDIM& dim1, size_t value1, ISMRMRDDIM& dim2, size_t value2)
{
    try
    {
        boost::shared_ptr< std::vector<size_t> > dimX = x.get_dimensions();

        long long dimInd1, dimInd2;
        GADGET_CHECK_RETURN_FALSE(getISMRMRDDimIndex(dim1, dimInd1));
        GADGET_CHECK_RETURN_FALSE(getISMRMRDDimIndex(dim2, dimInd2));

        GADGET_CHECK_RETURN_FALSE(value1<(*dimX)[dimInd1]);
        GADGET_CHECK_RETURN_FALSE(value2<(*dimX)[dimInd2]);

        std::vector<size_t> crop_offset(11, 0);
        crop_offset[0] = 0;
        crop_offset[1] = 0;
        crop_offset[2] = 0;
        crop_offset[3] = 0;
        crop_offset[4] = 0;
        crop_offset[5] = 0;
        crop_offset[6] = 0;
        crop_offset[7] = 0;
        crop_offset[8] = 0;
        crop_offset[9] = 0;
        crop_offset[10] = 0;

        std::vector<size_t> crop_size(11, 0);
        crop_size[0] = (*dimX)[0];
        crop_size[1] = (*dimX)[1];
        crop_size[2] = (*dimX)[2];
        crop_size[3] = (*dimX)[3];
        crop_size[4] = (*dimX)[4];
        crop_size[5] = (*dimX)[5];
        crop_size[6] = (*dimX)[6];
        crop_size[7] = (*dimX)[7];
        crop_size[8] = (*dimX)[8];
        crop_size[9] = (*dimX)[9];
        crop_size[10] = (*dimX)[9];

        crop_size[dimInd1] = value1+1;

        crop_offset[dimInd2] = value2;
        crop_size[dimInd2] = 1;

        GADGET_CHECK_RETURN_FALSE(Gadgetron::cropUpTo11DArray(x, r, crop_offset, crop_size));
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::extractSubArrayForDim1LessEqualDim2Equal(dim1, value1, dim2, value2) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::
extractSubArrayForMaxEncodingCounters(const hoNDArray<T>& x, hoNDArray<T>& r, const ISMRMRD::EncodingCounters& maxIdx)
{
    try
    {
        boost::shared_ptr< std::vector<size_t> > dimX = x.get_dimensions();

        std::vector<size_t> crop_offset(11, 0);
        crop_offset[0] = 0;
        crop_offset[1] = 0;
        crop_offset[2] = 0;
        crop_offset[3] = 0;
        crop_offset[4] = 0;
        crop_offset[5] = 0;
        crop_offset[6] = 0;
        crop_offset[7] = 0;
        crop_offset[8] = 0;
        crop_offset[9] = 0;
        crop_offset[10] = 0;

        // [RO E1 Cha Slice E2 Contrast Phase Rep Set Seg Ave]
        std::vector<size_t> crop_size(11, 0);
        crop_size[0] = (*dimX)[0];
        crop_size[1] = (*dimX)[1]; if ( maxIdx.kspace_encode_step_1 < crop_size[1]-1 ) crop_size[1] = maxIdx.kspace_encode_step_1+1;
        crop_size[2] = (*dimX)[2]; 
        crop_size[3] = (*dimX)[3]; if ( maxIdx.slice                < crop_size[3]-1 ) crop_size[3] = maxIdx.slice+1;
        crop_size[4] = (*dimX)[4]; if ( maxIdx.kspace_encode_step_2 < crop_size[4]-1 ) crop_size[4] = maxIdx.kspace_encode_step_2+1;
        crop_size[5] = (*dimX)[5]; if ( maxIdx.contrast             < crop_size[5]-1 ) crop_size[5] = maxIdx.contrast+1;
        crop_size[6] = (*dimX)[6]; if ( maxIdx.phase                < crop_size[6]-1 ) crop_size[6] = maxIdx.phase+1;
        crop_size[7] = (*dimX)[7]; if ( maxIdx.repetition           < crop_size[7]-1 ) crop_size[7] = maxIdx.repetition+1;
        crop_size[8] = (*dimX)[8]; if ( maxIdx.set                  < crop_size[8]-1 ) crop_size[8] = maxIdx.set+1;
        crop_size[9] = (*dimX)[9]; if ( maxIdx.segment              < crop_size[9]-1 ) crop_size[9] = maxIdx.segment+1;
        crop_size[10] = (*dimX)[10]; if ( maxIdx.average            < crop_size[10]-1 ) crop_size[10] = maxIdx.average+1;

        GADGET_CHECK_RETURN_FALSE(Gadgetron::cropUpTo11DArray(x, r, crop_offset, crop_size));
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::extractSubArrayForMaxEncodingCounters(const hoNDArray<T>& x, hoNDArray<T>& r, const ISMRMRD::EncodingCounters& maxIdx) ... ");
        return false;
    }

    return true;
}

template <typename T> 
void gtPlusISMRMRDReconUtil<T>::clearAcquisitionHeaderISMRMRD(ISMRMRD::AcquisitionHeader& acqHeader)
{
    memset(&acqHeader, 0, sizeof(ISMRMRD::AcquisitionHeader));
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::hasIdenticalGeometryISMRMRD(const ISMRMRD::AcquisitionHeader& acqHeader1, const ISMRMRD::AcquisitionHeader& acqHeader2)
{
    long long ii;

    for ( ii=0; ii<ISMRMRD::ISMRMRD_POSITION_LENGTH; ii++ )
    {
        if ( std::abs(acqHeader1.position[ii]-acqHeader2.position[ii]) > GT_IMAGING_GEOMETRY_DELTA ) return false;
        if ( std::abs(acqHeader1.patient_table_position[ii]-acqHeader2.patient_table_position[ii]) > GT_IMAGING_GEOMETRY_DELTA ) return false;
    }

    for ( ii=0; ii<ISMRMRD::ISMRMRD_DIRECTION_LENGTH; ii++ )
    {
        if ( std::abs(acqHeader1.read_dir[ii]-acqHeader2.read_dir[ii]) > GT_IMAGING_GEOMETRY_DELTA ) return false;
        if ( std::abs(acqHeader1.phase_dir[ii]-acqHeader2.phase_dir[ii]) > GT_IMAGING_GEOMETRY_DELTA ) return false;
        if ( std::abs(acqHeader1.slice_dir[ii]-acqHeader2.slice_dir[ii]) > GT_IMAGING_GEOMETRY_DELTA ) return false;
    }

    return true;
}

template <typename T> 
long long gtPlusISMRMRDReconUtil<T>::addPrePostZeros(size_t centre_column, size_t samples)
{
    // 1 : pre zeros
    // 2 : post zeros
    // 0 : no zeros
    if ( 2*centre_column == samples )
    {
        return 0;
    }

    if ( 2*centre_column < samples )
    {
        return 1;
    }

    if ( 2*centre_column > samples )
    {
        return 2;
    }

    return 0;
}

template <typename T> 
void gtPlusISMRMRDReconUtil<T>::findStartEndRO(size_t centre_column, size_t samples, long long& startRO, long long& endRO)
{
    long long zerosFlag = addPrePostZeros(centre_column, samples);

    if ( zerosFlag == 0 )
    {
        startRO = 0;
        endRO = (long long)samples-1;
    }

    if ( zerosFlag == 1 )
    {
        endRO = (long long)2*(samples-centre_column)-1;
        startRO = (long long)endRO-samples+1;
    }

    if ( zerosFlag == 2 )
    {
        startRO = 0;
        endRO = (long long)samples-1;
    }

    return;
}

template <typename T> 
void gtPlusISMRMRDReconUtil<T>::findStartEndROAfterZeroFilling(size_t centre_column, size_t samples_zerofilled, int& startRO, int& endRO)
{
    size_t num = samples_zerofilled/2;

    if ( centre_column == num )
    {
        startRO = 0;
        endRO = (int)samples_zerofilled-1;
    }

    if ( centre_column+num < samples_zerofilled ) // pre zeros
    {
        endRO = (int)samples_zerofilled-1;
        startRO = endRO-(int)(centre_column+num)+1;
    }

    if ( centre_column+num > samples_zerofilled ) // post zeros
    {
        startRO = 0;
        endRO = (int)samples_zerofilled-1;
    }

    return;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::setMetaAttributesFromImageHeaderISMRMRD(const ISMRMRD::ImageHeader& imgHeader, ISMRMRD::MetaContainer& attrib)
{
    try
    {
        unsigned int ii;

        attrib.set(ISMRMRD_IMAGE_version,                 (long)imgHeader.version);
        attrib.set(ISMRMRD_IMAGE_flags,                   (long)imgHeader.flags);
        attrib.set(ISMRMRD_IMAGE_measurement_uid,         (long)imgHeader.measurement_uid);

        // ----------------------------------------

        attrib.set(ISMRMRD_IMAGE_matrix_size, (long)imgHeader.matrix_size[0]);
        attrib.append(ISMRMRD_IMAGE_matrix_size, (long)imgHeader.matrix_size[1]);
        attrib.append(ISMRMRD_IMAGE_matrix_size, (long)imgHeader.matrix_size[2]);

        // ----------------------------------------

        attrib.set(ISMRMRD_IMAGE_field_of_view, (double)imgHeader.field_of_view[0]);
        attrib.append(ISMRMRD_IMAGE_field_of_view, (double)imgHeader.field_of_view[1]);
        attrib.append(ISMRMRD_IMAGE_field_of_view, (double)imgHeader.field_of_view[2]);

        // ----------------------------------------

        attrib.set(ISMRMRD_IMAGE_channels, (long)imgHeader.channels);

        // ----------------------------------------

        attrib.set(ISMRMRD_IMAGE_position, (double)imgHeader.position[0]);
        for ( ii=1; ii<ISMRMRD::ISMRMRD_POSITION_LENGTH; ii++ )
        {
            attrib.append(ISMRMRD_IMAGE_position, (double)imgHeader.position[ii]);
        }

        // ----------------------------------------

        attrib.set(ISMRMRD_IMAGE_read_dir, (double)imgHeader.read_dir[0]);
        for ( ii=1; ii<ISMRMRD::ISMRMRD_DIRECTION_LENGTH; ii++ )
        {
            attrib.append(ISMRMRD_IMAGE_read_dir, (double)imgHeader.read_dir[ii]);
        }

        // ----------------------------------------

        attrib.set(ISMRMRD_IMAGE_phase_dir, (double)imgHeader.phase_dir[0]);
        for ( ii=1; ii<ISMRMRD::ISMRMRD_DIRECTION_LENGTH; ii++ )
        {
            attrib.append(ISMRMRD_IMAGE_phase_dir, (double)imgHeader.phase_dir[ii]);
        }

        // ----------------------------------------

        attrib.set(ISMRMRD_IMAGE_slice_dir, (double)imgHeader.slice_dir[0]);
        for ( ii=1; ii<ISMRMRD::ISMRMRD_DIRECTION_LENGTH; ii++ )
        {
            attrib.append(ISMRMRD_IMAGE_slice_dir, (double)imgHeader.slice_dir[ii]);
        }

        // ----------------------------------------

        attrib.set(ISMRMRD_IMAGE_patient_table_position, (double)imgHeader.patient_table_position[0]);
        for ( ii=1; ii<ISMRMRD::ISMRMRD_POSITION_LENGTH; ii++ )
        {
            attrib.append(ISMRMRD_IMAGE_patient_table_position, (double)imgHeader.patient_table_position[ii]);
        }

        // ----------------------------------------

        attrib.set(ISMRMRD_IMAGE_average,       (long)imgHeader.average);
        attrib.set(ISMRMRD_IMAGE_slice,         (long)imgHeader.slice);
        attrib.set(ISMRMRD_IMAGE_contrast,      (long)imgHeader.contrast);
        attrib.set(ISMRMRD_IMAGE_phase,         (long)imgHeader.phase);
        attrib.set(ISMRMRD_IMAGE_repetition,    (long)imgHeader.repetition);
        attrib.set(ISMRMRD_IMAGE_set,           (long)imgHeader.set);

        // ----------------------------------------

        attrib.set(ISMRMRD_IMAGE_acquisition_time_stamp, (long)imgHeader.acquisition_time_stamp);

        // ----------------------------------------

        attrib.set(ISMRMRD_IMAGE_physiology_time_stamp, (long)imgHeader.physiology_time_stamp[0]);
        for ( ii=1; ii<ISMRMRD::ISMRMRD_PHYS_STAMPS; ii++ )
        {
            attrib.append(ISMRMRD_IMAGE_physiology_time_stamp, (long)imgHeader.physiology_time_stamp[ii]);
        }

        // ----------------------------------------

        attrib.set(ISMRMRD_IMAGE_image_data_type,       (long)imgHeader.data_type);
        attrib.set(ISMRMRD_IMAGE_image_type,            (long)imgHeader.image_type);
        attrib.set(ISMRMRD_IMAGE_image_series_index,    (long)imgHeader.image_series_index);

        // ----------------------------------------

        attrib.set(ISMRMRD_IMAGE_user_int, (long)imgHeader.user_int[0]);
        for ( ii=1; ii<ISMRMRD::ISMRMRD_USER_INTS; ii++ )
        {
            attrib.append(ISMRMRD_IMAGE_user_int, (long)imgHeader.user_int[ii]);
        }

        // ----------------------------------------

        attrib.set(ISMRMRD_IMAGE_user_float, (double)imgHeader.user_float[0]);
        for ( ii=1; ii<ISMRMRD::ISMRMRD_USER_FLOATS; ii++ )
        {
            attrib.append(ISMRMRD_IMAGE_user_float, (double)imgHeader.user_float[ii]);
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::setMetaAttributesFromImageHeaderISMRMRD(const ISMRMRD::ImageHeader& imgHeader, ISMRMRD::MetaContainer& attrib) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::setImageHeaderISMRMRDFromMetaAttributes(const ISMRMRD::MetaContainer& attrib, ISMRMRD::ImageHeader& imgHeader)
{
    try
    {
        unsigned int ii;

        imgHeader.version = (uint16_t)attrib.as_long(ISMRMRD_IMAGE_version, 0);
        imgHeader.flags = (uint64_t)attrib.as_long(ISMRMRD_IMAGE_flags, 0);
        imgHeader.measurement_uid = (uint32_t)attrib.as_long(ISMRMRD_IMAGE_measurement_uid, 0);

        // ----------------------------------------

        imgHeader.matrix_size[0] = (uint16_t)attrib.as_long(ISMRMRD_IMAGE_matrix_size, 0);
        imgHeader.matrix_size[1] = (uint16_t)attrib.as_long(ISMRMRD_IMAGE_matrix_size, 1);
        imgHeader.matrix_size[2] = (uint16_t)attrib.as_long(ISMRMRD_IMAGE_matrix_size, 2);

        // ----------------------------------------

        imgHeader.field_of_view[0] = (float)attrib.as_double(ISMRMRD_IMAGE_field_of_view, 0);
        imgHeader.field_of_view[1] = (float)attrib.as_double(ISMRMRD_IMAGE_field_of_view, 1);
        imgHeader.field_of_view[2] = (float)attrib.as_double(ISMRMRD_IMAGE_field_of_view, 2);

        // ----------------------------------------

        imgHeader.channels = (uint16_t)attrib.as_long(ISMRMRD_IMAGE_channels, 0);;

        // ----------------------------------------

        for ( ii=0; ii<ISMRMRD::ISMRMRD_POSITION_LENGTH; ii++ )
        {
            imgHeader.position[ii] = (float)attrib.as_double(ISMRMRD_IMAGE_position, ii);
        }

        // ----------------------------------------

        for ( ii=0; ii<ISMRMRD::ISMRMRD_DIRECTION_LENGTH; ii++ )
        {
            imgHeader.read_dir[ii] = (float)attrib.as_double(ISMRMRD_IMAGE_read_dir, ii);
        }

        // ----------------------------------------

        for ( ii=0; ii<ISMRMRD::ISMRMRD_DIRECTION_LENGTH; ii++ )
        {
            imgHeader.phase_dir[ii] = (float)attrib.as_double(ISMRMRD_IMAGE_phase_dir, ii);
        }

        // ----------------------------------------

        for ( ii=0; ii<ISMRMRD::ISMRMRD_DIRECTION_LENGTH; ii++ )
        {
            imgHeader.slice_dir[ii] = (float)attrib.as_double(ISMRMRD_IMAGE_slice_dir, ii);
        }

        // ----------------------------------------

        for ( ii=0; ii<ISMRMRD::ISMRMRD_POSITION_LENGTH; ii++ )
        {
            imgHeader.patient_table_position[ii] = (float)attrib.as_double(ISMRMRD_IMAGE_patient_table_position, ii);
        }

        // ----------------------------------------

        imgHeader.average = (uint16_t)attrib.as_long(ISMRMRD_IMAGE_average, 0);
        imgHeader.slice = (uint16_t)attrib.as_long(ISMRMRD_IMAGE_slice, 0);
        imgHeader.contrast = (uint16_t)attrib.as_long(ISMRMRD_IMAGE_contrast, 0);
        imgHeader.phase = (uint16_t)attrib.as_long(ISMRMRD_IMAGE_phase, 0);
        imgHeader.repetition = (uint16_t)attrib.as_long(ISMRMRD_IMAGE_repetition, 0);
        imgHeader.set = (uint16_t)attrib.as_long(ISMRMRD_IMAGE_set, 0);

        // ----------------------------------------

        imgHeader.acquisition_time_stamp = (uint32_t)attrib.as_long(ISMRMRD_IMAGE_acquisition_time_stamp, 0);

        // ----------------------------------------

        for ( ii=0; ii<ISMRMRD::ISMRMRD_PHYS_STAMPS; ii++ )
        {
            imgHeader.physiology_time_stamp[ii] = (uint32_t)attrib.as_long(ISMRMRD_IMAGE_physiology_time_stamp, ii);
        }

        // ----------------------------------------

        imgHeader.data_type = (uint16_t)attrib.as_long(ISMRMRD_IMAGE_image_data_type, 0);
        imgHeader.image_type = (uint16_t)attrib.as_long(ISMRMRD_IMAGE_image_type, 0);
        imgHeader.image_series_index = (uint16_t)attrib.as_long(ISMRMRD_IMAGE_image_series_index, 0);

        // ----------------------------------------

        for ( ii=0; ii<ISMRMRD::ISMRMRD_USER_INTS; ii++ )
        {
            imgHeader.user_int[ii] = (int32_t)attrib.as_long(ISMRMRD_IMAGE_user_int, ii);
        }

        // ----------------------------------------

        for ( ii=0; ii<ISMRMRD::ISMRMRD_USER_FLOATS; ii++ )
        {
            imgHeader.user_float[ii] = (float)attrib.as_double(ISMRMRD_IMAGE_user_float, ii);
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::setImageHeaderISMRMRDFromMetaAttributes(const ISMRMRD::MetaContainer& attrib, ISMRMRD::ImageHeader& imgHeader) ... ");
        return false;
    }

    return true;
}

//#ifdef USE_CUDA
//
//template <typename T> 
//bool gtPlusISMRMRDReconUtil<T>::
//cudaJobSplitter(const std::vector<unsigned int>& jobIDs, size_t jobSize, size_t minimalMemoryForValidDevice, 
//                std::vector< std::pair<unsigned int, std::vector<std::vector<unsigned int> > > >& jobSchedule)
//{
//    try
//    {
//        unsigned int numOfJobs = jobIDs.size();
//        if ( numOfJobs == 0 )
//        {
//            GWARN_STREAM("numOfJobs == 0");
//            return true;
//        }
//
//        // find valid device
//        int numOfDevices(0);
//        GADGET_CHECK_RETURN_FALSE(cudaGetDeviceCount( &numOfDevices )==cudaSuccess);
//
//        if ( numOfDevices == 0 )
//        {
//            GWARN_STREAM("numOfDevices == 0");
//            return true;
//        }
//
//        std::vector<unsigned int> validDevices;
//        int d;
//        for ( d=0; d<numOfDevices; d++ )
//        {
//            size_t totalMem = cudaDeviceManager::Instance()->total_global_mem(d);
//            if ( totalMem >= minimalMemoryForValidDevice )
//            {
//                validDevices.push_back(d);
//            }
//        }
//
//        if ( validDevices.empty() )
//        {
//            GERROR_STREAM("No valid device can be found : " << minimalMemoryForValidDevice);
//            return false;
//        }
//
//        std::vector<unsigned int> maxJobN(validDevices.size());
//        for ( d=0; d<validDevices.size(); d++ )
//        {
//            size_t totalMem = cudaDeviceManager::Instance()->total_global_mem(validDevices[d]);
//            maxJobN[d] = totalMem/jobSize;
//        }
//
//        jobSchedule.clear();
//
//        size_t job = 0;
//        unsigned int validDevice = 0;
//        while ( job < numOfJobs )
//        {
//            size_t start = job;
//            size_t end = job + maxJobN[validDevice] - 1;
//
//            if ( end >= numOfJobs ) end = numOfJobs - 1;
//
//            unsigned int deviceID = validDevices[validDevice];
//
//            unsigned int loc;
//            for ( loc=0; loc<jobSchedule.size(); loc++ )
//            {
//                if ( jobSchedule[loc].first == deviceID ) break;
//            }
//
//            if ( loc < jobSchedule.size() )
//            {
//                // insert a new job package
//                std::vector<unsigned int> jobPackage;
//                for ( unsigned int jj=start; jj<=end; jj++ )
//                {
//                    jobPackage.push_back(jobIDs[jj]);
//                }
//
//                jobSchedule[loc].second.push_back(jobPackage);
//            }
//            else
//            {
//                // create a new entry
//                std::pair<unsigned int, std::vector<std::vector<unsigned int> > > jobItem;
//                jobItem.first = deviceID;
//
//                std::vector<unsigned int> jobPackage;
//                for ( unsigned int jj=start; jj<=end; jj++ )
//                {
//                    jobPackage.push_back(jobIDs[jj]);
//                }
//                jobItem.second.push_back(jobPackage);
//
//                jobSchedule.push_back(jobItem);
//            }
//
//            job = end+1;
//            validDevice++;
//
//            if ( validDevice >= validDevices.size() )
//            {
//                validDevice = 0;
//            }
//        }
//    }
//    catch(...)
//    {
//        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::cudaJobSplitter(...) ... ");
//        return false;
//    }
//
//    return true;
//}
//
//template <typename T> 
//bool gtPlusISMRMRDReconUtil<T>::
//cudaJobSplitter(unsigned int numOfJobs, size_t jobSize, size_t minimalMemoryForValidDevice, 
//            std::vector< std::pair<unsigned int, std::vector<std::vector<unsigned int> > > >& jobSchedule)
//{
//    if ( numOfJobs == 0 )
//    {
//        GWARN_STREAM("numOfJobs == 0");
//        return true;
//    }
//
//    std::vector<unsigned int> jobIDs(numOfJobs, 0);
//    unsigned int ii;
//    for ( ii=0; ii<numOfJobs; ii++ ) jobIDs[ii] = ii;
//    return cudaJobSplitter(jobIDs, jobSize, minimalMemoryForValidDevice, jobSchedule);
//}
//
//#endif // USE_CUDA

template <typename T> 
void gtPlusISMRMRDReconUtil<T>::
compareAgainstGroundTruthArray(const std::string& gt_filename, const hoNDArray<T>& x, typename realType<T>::Type& normDiff, typename realType<T>::Type& maxNormDiff)
{
    hoNDArray<T> gt;

    gtPlusIOAnalyze gt_io;
    gt_io.importArray(gt, gt_filename);

    compareAgainstGroundTruthArray(gt, x, normDiff, maxNormDiff);
}

template <typename T> 
void gtPlusISMRMRDReconUtil<T>::
compareAgainstGroundTruthArray(const hoNDArray<T>& gt, const hoNDArray<T>& x, typename realType<T>::Type& normDiff, typename realType<T>::Type& maxNormDiff)
{
    hoNDArray<T> diff(x);
    Gadgetron::subtract(gt, x, diff);

    typename realType<T>::Type v;
    Gadgetron::norm2(diff, v);
    normDiff = v;

    T maxV;
    size_t ind;
    Gadgetron::maxAbsolute(diff, maxV, ind);
    maxNormDiff = std::abs(maxV);
}

// ========================================================================================== //

template <typename T> 
gtPlusISMRMRDReconUtilComplex<T>::gtPlusISMRMRDReconUtilComplex() {}

template <typename T> 
gtPlusISMRMRDReconUtilComplex<T>::~gtPlusISMRMRDReconUtilComplex() {}

template <typename T> 
void gtPlusISMRMRDReconUtilComplex<T>::printInfo(std::ostream& os)
{
    using namespace std;
    os << "-------------- GTPlus ISMRMRD Recon Util Complex -------------" << endl;
    os << "Implementation of recon utilities for ISMRMRD complex data type" << endl;
    os << "--------------------------------------------------------------" << endl;
}

// ------------------------------------------------------------------------
// noise prewhitening
// ------------------------------------------------------------------------
template <typename T> 
bool gtPlusISMRMRDReconUtilComplex<T>::
computeNoisePrewhiteningMatrix(const hoNDArray<T>& noise, double noiseBandWidth, double receiverBWRatio, double ADCSamplingTimeinSecond, hoMatrix<T>& prewhiteningMatrix)
{
    try
    {
        size_t RO = noise.get_size(0);
        size_t E1 = noise.get_size(1);
        size_t CHA = noise.get_size(2);

        GADGET_CHECK_RETURN_FALSE(prewhiteningMatrix.createMatrix(CHA, CHA));
        Gadgetron::clear(prewhiteningMatrix);

        typedef typename realType<T>::Type ValueType;

        // noise sampling time in second
        ValueType noiseSamplingTimeinSecond = (ValueType)(1.0/(noiseBandWidth*RO));

        // scaling factor
        ValueType scaling = (ValueType)(noiseSamplingTimeinSecond/ADCSamplingTimeinSecond/receiverBWRatio);
        scaling /= (RO*E1-1);

        // compute the noise covariance matrix
        hoMatrix<T> R(RO*E1, CHA, const_cast<T*>(noise.begin()));

        // R'*R --> CHA by CHA covariance matrix
        Gadgetron::gemm(prewhiteningMatrix, R, true, R, false);
        Gadgetron::scal(scaling, prewhiteningMatrix);

        // 0.5*(R+R')
        hoMatrix<T> RH(prewhiteningMatrix);
        conjugatetrans(prewhiteningMatrix, RH);
        Gadgetron::add(prewhiteningMatrix, RH, prewhiteningMatrix);
        Gadgetron::scal( (ValueType)0.5, prewhiteningMatrix);

        Gadgetron::potrf(prewhiteningMatrix, 'U');
        Gadgetron::trtri(prewhiteningMatrix, 'U');
        Gadgetron::scal( (value_type)(std::sqrt((double)2.0)), prewhiteningMatrix);
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtilComplex<T>::computeNoisePrewhiteningMatrix(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtilComplex<T>::
performNoisePrewhitening(hoNDArray<T>& data, const hoMatrix<T>& prewhiteningMatrix)
{
    try
    {
        size_t RO = data.get_size(0);
        size_t E1 = data.get_size(1);
        size_t CHA = data.get_size(2);

        GADGET_CHECK_RETURN_FALSE(prewhiteningMatrix.rows()==CHA);
        GADGET_CHECK_RETURN_FALSE(prewhiteningMatrix.cols()==CHA);

        size_t N = data.get_number_of_elements()/(RO*E1*CHA);

        long long n;
        #pragma omp parallel default(none) private(n) shared(RO, E1, CHA, N, data, prewhiteningMatrix)
        {
            hoMatrix<T> tmp(RO*E1, CHA);
            Gadgetron::clear(tmp);

            #pragma omp for
            for ( n=0; n<(long long)N; n++ )
            {
                hoMatrix<T> D(RO*E1, CHA, data.begin()+n*RO*E1*CHA);
                Gadgetron::gemm(tmp, D, false, prewhiteningMatrix, false);
                memcpy(data.begin()+n*RO*E1*CHA, tmp.begin(), sizeof(T)*RO*E1*CHA);
            }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtilComplex<T>::performNoisePrewhitening(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtilComplex<T>::
zpadResize2D(const hoNDArray<T>& data, size_t sizeX, size_t sizeY, hoNDArray<T>& dataResized)
{
    try
    {
        size_t RO = data.get_size(0);
        size_t E1 = data.get_size(1);

        GADGET_CHECK_RETURN_FALSE(sizeX>=RO);
        GADGET_CHECK_RETURN_FALSE(sizeY>=E1);

        if ( RO==sizeX && E1==sizeY )
        {
            dataResized = data;
            return true;
        }

        if ( dataResized.get_size(0)!=sizeX || dataResized.get_size(1)!=sizeY )
        {
            dataResized.create(sizeX, sizeY);
        }

        Gadgetron::clear(&dataResized);

        hoNDArray<T> kspace(data);
        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft2c(data, kspace);
        GADGET_CHECK_RETURN_FALSE(zpadResize2DOnKSpace(kspace, sizeX, sizeY, dataResized));
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtilComplex<T>::zpadResize2D(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtilComplex<T>::
zpadResize2DOnKSpace(const hoNDArray<T>& kspace, size_t sizeX, size_t sizeY, hoNDArray<T>& dataResized)
{
    try
    {
        size_t RO = kspace.get_size(0);
        size_t E1 = kspace.get_size(1);

        GADGET_CHECK_RETURN_FALSE(sizeX>=RO);
        GADGET_CHECK_RETURN_FALSE(sizeY>=E1);

        if ( RO==sizeX && E1==sizeY )
        {
            dataResized = kspace;
            Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft2c(dataResized);
            return true;
        }

        if ( dataResized.get_size(0)!=sizeX || dataResized.get_size(1)!=sizeY )
        {
            dataResized.create(sizeX, sizeY);
        }

        Gadgetron::clear(&dataResized);

        // GADGET_CHECK_RETURN_FALSE(this->zeropad2D(kspace, sizeX, sizeY, dataResized));
        // GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::zeropad2D(kspace, sizeX, sizeY, dataResized));
        GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::pad(sizeX, sizeY, const_cast<hoNDArray<T>*>(&kspace), &dataResized));
        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft2c(dataResized);

        typename realType<T>::Type scaling = (typename realType<T>::Type)(std::sqrt((double)sizeX*sizeY)/std::sqrt((double)RO*E1));
        Gadgetron::scal(scaling, dataResized);
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtilComplex<T>::zpadResize2DOnKSpace(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtilComplex<T>::
zpadResize3D(const hoNDArray<T>& data, size_t sizeX, size_t sizeY, size_t sizeZ, hoNDArray<T>& dataResized)
{
    try
    {
        size_t RO = data.get_size(0);
        size_t E1 = data.get_size(1);
        size_t E2 = data.get_size(2);

        GADGET_CHECK_RETURN_FALSE(sizeX>=RO);
        GADGET_CHECK_RETURN_FALSE(sizeY>=E1);
        GADGET_CHECK_RETURN_FALSE(sizeZ>=E2);

        if ( RO==sizeX && E1==sizeY && E2==sizeZ )
        {
            dataResized = data;
            return true;
        }

        if ( dataResized.get_size(0)!=sizeX || dataResized.get_size(1)!=sizeY || dataResized.get_size(2)!=sizeZ )
        {
            dataResized.create(sizeX, sizeY, sizeZ);
        }

        Gadgetron::clear(&dataResized);

        hoNDArray<T> kspace(data);
        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft3c(data, kspace);
        GADGET_CHECK_RETURN_FALSE(zpadResize3DOnKSpace(kspace, sizeX, sizeY, sizeZ, dataResized));
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtilComplex<T>::zpadResize3D(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtilComplex<T>::
zpadResize3DOnKSpace(const hoNDArray<T>& kspace, size_t sizeX, size_t sizeY, size_t sizeZ, hoNDArray<T>& dataResized)
{
    try
    {
        size_t RO = kspace.get_size(0);
        size_t E1 = kspace.get_size(1);
        size_t E2 = kspace.get_size(2);

        GADGET_CHECK_RETURN_FALSE(sizeX>=RO);
        GADGET_CHECK_RETURN_FALSE(sizeY>=E1);
        GADGET_CHECK_RETURN_FALSE(sizeZ>=E2);

        if ( RO==sizeX && E1==sizeY && E2==sizeZ )
        {
            dataResized = kspace;
            Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft3c(dataResized);
            return true;
        }

        if ( dataResized.get_size(0)!=sizeX || dataResized.get_size(1)!=sizeY || dataResized.get_size(2)!=sizeZ )
        {
            dataResized.create(sizeX, sizeY, sizeZ);
        }

        Gadgetron::clear(&dataResized);

        // GADGET_CHECK_RETURN_FALSE(this->zeropad3D(kspace, sizeX, sizeY, sizeZ, dataResized));
        // GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::zeropad3D(kspace, sizeX, sizeY, sizeZ, dataResized));
        GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::pad(sizeX, sizeY, sizeZ, const_cast<hoNDArray<T>*>(&kspace), &dataResized));
        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft3c(dataResized);

        typename realType<T>::Type scaling = (typename realType<T>::Type)(std::sqrt((double)sizeX*sizeY*sizeZ)/std::sqrt((double)RO*E1*E2));
        Gadgetron::scal(scaling, dataResized);
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtilComplex<T>::zpadResize3DOnKSpace(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtilComplex<T>::
zpadResize2DFilter(const hoNDArray<T>& data, size_t sizeX, size_t sizeY, const hoNDArray<T>& filter2D, hoNDArray<T>& dataResized)
{
    try
    {
        size_t RO = data.get_size(0);
        size_t E1 = data.get_size(1);

        GADGET_CHECK_RETURN_FALSE(sizeX>=RO);
        GADGET_CHECK_RETURN_FALSE(sizeY>=E1);

        GADGET_CHECK_RETURN_FALSE(filter2D.get_size(0)==sizeX);
        GADGET_CHECK_RETURN_FALSE(filter2D.get_size(1)==sizeY);

        if ( RO==sizeX && E1==sizeY )
        {
            dataResized = data;
            return true;
        }

        if ( dataResized.get_size(0)!=sizeX || dataResized.get_size(1)!=sizeY )
        {
            dataResized.create(sizeX, sizeY);
        }

        Gadgetron::clear(&dataResized);

        hoNDArray<T> kspace(data);
        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft2c(data, kspace);
        // GADGET_CHECK_RETURN_FALSE(this->zeropad2D(kspace, sizeX, sizeY, dataResized));
        // GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::zeropad2D(kspace, sizeX, sizeY, dataResized));
        GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::pad(sizeX, sizeY, const_cast<hoNDArray<T>*>(&kspace), &dataResized));
        GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::multiply(dataResized, filter2D, dataResized));
        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft2c(dataResized);

        typename realType<T>::Type scaling = (typename realType<T>::Type)(std::sqrt((double)sizeX*sizeY)/std::sqrt((double)RO*E1));
        Gadgetron::scal(scaling, dataResized);
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtilComplex<T>::zpadResize2DFilter(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtilComplex<T>::
zpadResize3DFilter(const hoNDArray<T>& data, size_t sizeX, size_t sizeY, size_t sizeZ, const hoNDArray<T>& filter3D, hoNDArray<T>& dataResized)
{
    try
    {
        size_t RO = data.get_size(0);
        size_t E1 = data.get_size(1);
        size_t E2 = data.get_size(2);

        GADGET_CHECK_RETURN_FALSE(sizeX>=RO);
        GADGET_CHECK_RETURN_FALSE(sizeY>=E1);
        GADGET_CHECK_RETURN_FALSE(sizeZ>=E2);

        GADGET_CHECK_RETURN_FALSE(filter3D.get_size(0)==sizeX);
        GADGET_CHECK_RETURN_FALSE(filter3D.get_size(1)==sizeY);
        GADGET_CHECK_RETURN_FALSE(filter3D.get_size(2)==sizeZ);

        if ( RO==sizeX && E1==sizeY && E2==sizeZ )
        {
            dataResized = data;
            return true;
        }

        if ( dataResized.get_size(0)!=sizeX || dataResized.get_size(1)!=sizeY || dataResized.get_size(2)!=sizeZ )
        {
            dataResized.create(sizeX, sizeY, sizeZ);
        }

        Gadgetron::clear(&dataResized);

        hoNDArray<T> kspace(data);
        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft3c(data, kspace);
        // GADGET_CHECK_RETURN_FALSE(this->zeropad3D(kspace, sizeX, sizeY, sizeZ, dataResized));
        // GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::zeropad3D(kspace, sizeX, sizeY, sizeZ, dataResized));
        GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::pad(sizeX, sizeY, sizeZ, &kspace, &dataResized));
        GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::multiply(dataResized, filter3D, dataResized));
        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft3c(dataResized);

        typename realType<T>::Type scaling = (typename realType<T>::Type)(std::sqrt((double)sizeX*sizeY*sizeZ)/std::sqrt((double)RO*E1*E2));
        Gadgetron::scal(scaling, dataResized);
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtilComplex<T>::zpadResize3DFilter(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtilComplex<T>::
kspacefilterROImage(hoNDArray<T>& data, const hoNDArray<T>& fRO)
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(data.get_size(0)==fRO.get_number_of_elements());
        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft1c(data);
        GADGET_CHECK_EXCEPTION_RETURN_FALSE(apply_kspace_filter_RO(data, fRO));
        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft1c(data);
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtilComplex<T>::kspacefilterROImage(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtilComplex<T>::
kspacefilterROImage(const hoNDArray<T>& data, const hoNDArray<T>& fRO, hoNDArray<T>& dataFiltered)
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(data.get_size(0)==fRO.get_number_of_elements());
        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft1c(data, dataFiltered);
        GADGET_CHECK_EXCEPTION_RETURN_FALSE(apply_kspace_filter_RO(dataFiltered, fRO));
        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft1c(dataFiltered);
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtilComplex<T>::kspacefilterROImage(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtilComplex<T>::
kspacefilterE1Image(const hoNDArray<T>& data, const hoNDArray<T>& fE1, hoNDArray<T>& dataFiltered)
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(data.get_size(1)==fE1.get_number_of_elements());

        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft2c(data, dataFiltered);
        GADGET_CHECK_EXCEPTION_RETURN_FALSE(apply_kspace_filter_RO(dataFiltered, fE1, dataFiltered));
        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft2c(dataFiltered);
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtilComplex<T>::kspacefilterE1Image(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtilComplex<T>::
kspacefilterE2Image(const hoNDArray<T>& data, const hoNDArray<T>& fE2, hoNDArray<T>& dataFiltered)
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(data.get_size(2)==fE2.get_number_of_elements());

        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft3c(data, dataFiltered);
        GADGET_CHECK_EXCEPTION_RETURN_FALSE(apply_kspace_filter_RO(dataFiltered, fE2, dataFiltered));
        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft3c(dataFiltered);
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtilComplex<T>::kspacefilterE2Image(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtilComplex<T>::
kspacefilterE1E2Image(const hoNDArray<T>& data, const hoNDArray<T>& fE1, const hoNDArray<T>& fE2, hoNDArray<T>& dataFiltered)
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(data.get_size(1)==fE1.get_number_of_elements());
        GADGET_CHECK_RETURN_FALSE(data.get_size(2)==fE2.get_number_of_elements());

        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft3c(data, dataFiltered);
        GADGET_CHECK_EXCEPTION_RETURN_FALSE(apply_kspace_filter_E1E2(dataFiltered, fE1, fE2, dataFiltered));
        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft3c(dataFiltered);
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtilComplex<T>::kspacefilterE1E2Image(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtilComplex<T>::
kspacefilterROE1E2Image(const hoNDArray<T>& data, const hoNDArray<T>& fRO, const hoNDArray<T>& fE1, const hoNDArray<T>& fE2, hoNDArray<T>& dataFiltered)
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(data.get_size(0)==fRO.get_number_of_elements());
        GADGET_CHECK_RETURN_FALSE(data.get_size(1)==fE1.get_number_of_elements());
        GADGET_CHECK_RETURN_FALSE(data.get_size(2)==fE2.get_number_of_elements());

        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft3c(data, dataFiltered);
        GADGET_CHECK_EXCEPTION_RETURN_FALSE(apply_kspace_filter_ROE1E2(dataFiltered, fRO, fE1, fE2, dataFiltered));
        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft3c(dataFiltered);
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtilComplex<T>::kspacefilterROE1E2Image(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtilComplex<T>::
coilMap2DNIH(const hoNDArray<T>& data, hoNDArray<T>& coilMap, ISMRMRDCOILMAPALGO algo, size_t ks, size_t power, size_t iterNum, typename realType<T>::Type thres, bool useGPU)
{
    try
    {
        typedef typename realType<T>::Type value_type;

        long long RO = data.get_size(0);
        long long E1 = data.get_size(1);
        long long CHA = data.get_size(2);

        size_t N = data.get_number_of_elements()/(RO*E1*CHA);
        size_t num = RO*E1*CHA;

        if ( !data.dimensions_equal(&coilMap) )
        {
            coilMap = data;
        }

        if ( ks%2 != 1 )
        {
            ks++;
        }

        long long n;

        hoNDArray<T> dataCurr(RO, E1, 1, CHA, N, const_cast<T*>(data.begin()));
        hoNDArray<T> coilMapCurr(RO, E1, 1, CHA, N, coilMap.begin());
        if (algo == ISMRMRD_SOUHEIL_ITER)
        {
            GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::coil_map_Inati_Iter(dataCurr, coilMapCurr, ks, ks, iterNum, thres));
        }
        else
        {
            GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::coil_map_Inati(dataCurr, coilMapCurr, ks, power));
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtilComplex<T>::coilMap2DNIH(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtilComplex<T>::
coilMap3DNIH(const hoNDArray<T>& data, hoNDArray<T>& coilMap, ISMRMRDCOILMAPALGO algo, size_t ks, size_t power, size_t iterNum, typename realType<T>::Type thres, bool true3D)
{
    try
    {
        typedef typename realType<T>::Type value_type;

        size_t RO = data.get_size(0);
        size_t E1 = data.get_size(1);
        size_t E2 = data.get_size(2);
        size_t CHA = data.get_size(3);

        size_t N = data.get_number_of_elements()/(RO*E1*E2*CHA);

        if ( !data.dimensions_equal(&coilMap) )
        {
            coilMap = data;
        }

        Gadgetron::clear(coilMap);

        if ( ks%2 != 1 )
        {
            ks++;
        }

        hoNDArray<T> dataCurr(RO, E1, E2, CHA, N, const_cast<T*>(data.begin()));
        hoNDArray<T> coilMapCurr(RO, E1, E2, CHA, N, coilMap.begin());
        if (algo == ISMRMRD_SOUHEIL_ITER)
        {
            GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::coil_map_Inati_Iter(dataCurr, coilMapCurr, ks, ks, iterNum, thres));
        }
        else
        {
            GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::coil_map_Inati(dataCurr, coilMapCurr, ks, power));
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtilComplex<T>::coilMap3DNIH(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtilComplex<T>::
sumOfSquare(const hoNDArray<T>& data, hoNDArray<T>& sos)
{
    try
    {
        size_t NDim = data.get_number_of_dimensions();
        GADGET_CHECK_RETURN_FALSE(NDim>=3);

        hoNDArray<T> tmp(data);
        Gadgetron::multiplyConj(data, data, tmp);

        if ( NDim == 3 )
        {
            // GADGET_CHECK_RETURN_FALSE(Gadgetron::sumOverLastDimension(tmp, sos));
            GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::sum_over_dimension(tmp, sos, 2));

            std::vector<size_t> dim(2);
            dim[0] = sos.get_size(0);
            dim[1] = sos.get_size(1);

            sos.reshape(dim);
        }
        else if ( NDim == 4 )
        {
            // GADGET_CHECK_RETURN_FALSE(Gadgetron::sumOverSecondLastDimension(tmp, sos));

            GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::sum_over_dimension(tmp, sos, 2));

            std::vector<size_t> dim(3);
            dim[0] = sos.get_size(0);
            dim[1] = sos.get_size(1);
            dim[2] = sos.get_size(3);

            sos.reshape(dim);
        }
        else
        {
            GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::sum_over_dimension(tmp, sos, 2));
        }

        Gadgetron::sqrt(sos, sos);
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtilComplex<T>::sumOfSquare(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtilComplex<T>::
coilCombine(const hoNDArray<T>& data, const hoNDArray<T>& coilMap, hoNDArray<T>& combined)
{
    try
    {
        Gadgetron::coil_combine(data, coilMap, 2, combined);
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtilComplex<T>::coilCombine(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtilComplex<T>::
coilCombine3D(const hoNDArray<T>& data, const hoNDArray<T>& coilMap, hoNDArray<T>& combined)
{
    try
    {
        Gadgetron::coil_combine(data, coilMap, 3, combined);
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtilComplex<T>::coilCombine3D(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtilComplex<T>::
conjugateSymmetry2D(const hoNDArray<T>& kspace, hoNDArray<T>& kspaceConj)
{
    try
    {
        if ( !kspaceConj.dimensions_equal(&kspace) )
        {
            kspaceConj.create(kspace.get_dimensions());
        }

        long long RO = kspace.get_size(0);
        long long E1 = kspace.get_size(1);
        long long num = kspace.get_number_of_elements()/(RO*E1);

        long long centerRO = RO/2;
        long long centerE1 = E1/2;

        long long ii;

        #pragma omp parallel for default(none) private(ii) shared(RO, E1, num, centerRO, centerE1, kspace, kspaceConj)
        for ( ii=0; ii<num; ii++ )
        {
            ho2DArray<T> src(RO, E1, const_cast<T*>(kspace.begin()+ii*RO*E1));
            ho2DArray<T> dst(RO, E1, const_cast<T*>(kspaceConj.begin()+ii*RO*E1));

            long long ro, e1;
            long long cro, ce1;

            for ( e1=0; e1<E1; e1++ )
            {
                ce1 = 2*centerE1-e1;
                if ( ce1 > E1-1 )
                {
                    ce1 -= E1;
                }

                for ( ro=0; ro<RO; ro++ )
                {
                    cro = 2*centerRO-ro;
                    if ( cro > RO-1 )
                    {
                        cro -= RO;
                    }

                    dst(ro, e1) = std::conj(src(cro, ce1));
                }
            }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtilComplex<T>::conjugateSymmetry2D(const hoNDArray<T>& kspace, hoNDArray<T>& kspaceConj) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtilComplex<T>::
conjugateSymmetry3D(const hoNDArray<T>& kspace, hoNDArray<T>& kspaceConj)
{
    try
    {
        if ( !kspaceConj.dimensions_equal(&kspace) )
        {
            kspaceConj.create(kspace.get_dimensions());
        }

        long long RO = kspace.get_size(0);
        long long E1 = kspace.get_size(1);
        long long E2 = kspace.get_size(2);
        long long num = kspace.get_number_of_elements()/(RO*E1*E2);

        long long centerRO = RO/2;
        long long centerE1 = E1/2;
        long long centerE2 = E2/2;

        long long ii;

        #pragma omp parallel for default(none) private(ii) shared(RO, E1, E2, num, centerRO, centerE1, centerE2, kspace, kspaceConj)
        for ( ii=0; ii<num; ii++ )
        {
            ho3DArray<T> src(RO, E1, E2, const_cast<T*>(kspace.begin()+ii*RO*E1*E2));
            ho3DArray<T> dst(RO, E1, E2, const_cast<T*>(kspaceConj.begin()+ii*RO*E1*E2));

            long long ro, e1, e2;
            long long cro, ce1, ce2;

            for ( e2=0; e2<E2; e2++ )
            {
                ce2 = 2*centerE2-e2;
                if ( ce2 > E2-1 )
                {
                    ce2 -= E2;
                }

                for ( e1=0; e1<E1; e1++ )
                {
                    ce1 = 2*centerE1-e1;
                    if ( ce1 > E1-1 )
                    {
                        ce1 -= E1;
                    }

                    for ( ro=0; ro<RO; ro++ )
                    {
                        cro = 2*centerRO-ro;
                        if ( cro > RO-1 )
                        {
                            cro -= RO;
                        }

                        dst(ro, e1, e2) = std::conj(src(cro, ce1, ce2));
                    }
                }
            }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtilComplex<T>::conjugateSymmetry3D(const hoNDArray<T>& kspace, hoNDArray<T>& kspaceConj) ... ");
        return false;
    }

    return true;
}

}}
