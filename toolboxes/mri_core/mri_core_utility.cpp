
/** \file   mri_core_utility.cpp
    \brief  Implementation useful utility functionalities for 2D and 3D MRI parallel imaging
    \author Hui Xue
*/

#include "mri_core_utility.h"

#include "hoNDFFT.h"
#include "hoNDArray_linalg.h"
#include "hoNDArray_utils.h"
#include "hoNDArray_elemwise.h"
#include "hoNDArray_reductions.h"

namespace Gadgetron
{

// ------------------------------------------------------------------------
// zero-padding resize
// ------------------------------------------------------------------------
void zpadRange(size_t srcSize, size_t dstSize, size_t& start, size_t& end)
{
    try
    {
        if ( srcSize >= dstSize )
        {
            start = 0;
            end = srcSize-1;
            return;
        }

        //unsigned srcCenterInd = srcSize/2;
        //unsigned dstCenterInd = dstSize/2;

        start = (dstSize/2) - (srcSize/2);
        end = srcSize + start -1;

        //start = std::floor((double)dstSize/2.0)+1+std::ceil(-1.0 * (double)srcSize/2.0)-1;
        //end = std::floor((double)dstSize/2.0)+std::ceil(srcSize/2.0)-1;
    }
    catch(...)
    {
        GADGET_THROW("Errors in zpadRange(...) ... ");
    }
}

// ------------------------------------------------------------------------

template <typename T> 
void zeropad2D(const hoNDArray<T>& data, size_t sizeX, size_t sizeY, hoNDArray<T>& dataPadded)
{
    try
    {
        size_t RO = data.get_size(0);
        size_t E1 = data.get_size(1);

        GADGET_CHECK_THROW(sizeX>=RO);
        GADGET_CHECK_THROW(sizeY>=E1);

        if ( RO==sizeX && E1==sizeY )
        {
            dataPadded = data;
            return;
        }

        size_t sRO, eRO, sE1, eE1;
        zpadRange(RO, sizeX, sRO, eRO);
        zpadRange(E1, sizeY, sE1, eE1);

        boost::shared_ptr< std::vector<size_t> > dimPadded = data.get_dimensions();
        (*dimPadded)[0] = sizeX;
        (*dimPadded)[1] = sizeY;
        dataPadded.create(dimPadded);
        Gadgetron::clear(&dataPadded);

        size_t num = data.get_number_of_elements()/(RO*E1);

        long long n;

        const T* pData = data.begin();
        T* pDataPadded = dataPadded.begin();

        #pragma omp parallel for default(none) private(n) shared(num, sE1, eE1, sRO, RO, E1, pData, pDataPadded, sizeX, sizeY)
        for ( n=0; n<(long long)num; n++ )
        {
            for ( size_t y=sE1; y<=eE1; y++ )
            {
                memcpy(pDataPadded+n*sizeX*sizeY+y*sizeX+sRO, pData+n*RO*E1+(y-sE1)*RO, sizeof(T)*RO);
            }
        }
    }
    catch(...)
    {
        GADGET_THROW("Errors in zeropad2D(...) ... ");
    }
}

template EXPORTMRICORE void zeropad2D(const hoNDArray<float>& data, size_t sizeX, size_t sizeY, hoNDArray<float>& dataPadded);
template EXPORTMRICORE void zeropad2D(const hoNDArray<double>& data, size_t sizeX, size_t sizeY, hoNDArray<double>& dataPadded);
template EXPORTMRICORE void zeropad2D(const hoNDArray< std::complex<float> >& data, size_t sizeX, size_t sizeY, hoNDArray< std::complex<float> >& dataPadded);
template EXPORTMRICORE void zeropad2D(const hoNDArray< std::complex<double> >& data, size_t sizeX, size_t sizeY, hoNDArray< std::complex<double> >& dataPadded);

// ------------------------------------------------------------------------

template <typename T> 
void zeropad3D(const hoNDArray<T>& data, size_t sizeX, size_t sizeY, size_t sizeZ, hoNDArray<T>& dataPadded)
{
    try
    {
        size_t RO = data.get_size(0);
        size_t E1 = data.get_size(1);
        size_t E2 = data.get_size(2);

        GADGET_CHECK_THROW(sizeX>=RO);
        GADGET_CHECK_THROW(sizeY>=E1);
        GADGET_CHECK_THROW(sizeZ>=E2);

        if ( RO==sizeX && E1==sizeY && E2==sizeZ )
        {
            dataPadded = data;
            return;
        }

        size_t sRO, eRO, sE1, eE1, sE2, eE2;
        zpadRange(RO, sizeX, sRO, eRO);
        zpadRange(E1, sizeY, sE1, eE1);
        zpadRange(E2, sizeZ, sE2, eE2);

        boost::shared_ptr< std::vector<size_t> > dimPadded = data.get_dimensions();
        (*dimPadded)[0] = sizeX;
        (*dimPadded)[1] = sizeY;
        (*dimPadded)[2] = sizeZ;
        dataPadded.create(dimPadded);
        Gadgetron::clear(&dataPadded);

        size_t num = data.get_number_of_elements()/(RO*E1*E2);

        long long n;

        const T* pData = data.begin();
        T* pDataPadded = dataPadded.begin();

        #pragma omp parallel for default(none) private(n) shared(num, sE2, eE2, sE1, eE1, sRO, RO, E1, E2, pData, pDataPadded, sizeX, sizeY, sizeZ)
        for ( n=0; n<(long long)num; n++ )
        {
            T* pDst = pDataPadded+n*sizeX*sizeY*sizeZ;
            T* pSrc = const_cast<T*>(pData)+n*RO*E1*E2;

            long long z;
            // #pragma omp parallel for default(none) private(z) shared(pDst, pSrc, sE2, eE2, sE1, eE1, sRO, RO, E1, E2, sizeX, sizeY, sizeZ) num_threads(2)
            for ( z=(long long)sE2; z<=(long long)eE2; z++ )
            {
                long long o1 = z*sizeX*sizeY + sRO;
                long long o2 = (z-sE2)*RO*E1;
                for ( size_t y=sE1; y<=eE1; y++ )
                {
                    memcpy(pDst+o1+y*sizeX, pSrc+o2+(y-sE1)*RO, sizeof(T)*RO);
                }
            }
        }
    }
    catch(...)
    {
        GADGET_THROW("Errors in zeropad3D(...) ... ");
    }
}

template EXPORTMRICORE void zeropad3D(const hoNDArray<float>& data, size_t sizeX, size_t sizeY, size_t sizeZ, hoNDArray<float>& dataPadded);
template EXPORTMRICORE void zeropad3D(const hoNDArray<double>& data, size_t sizeX, size_t sizeY, size_t sizeZ, hoNDArray<double>& dataPadded);
template EXPORTMRICORE void zeropad3D(const hoNDArray< std::complex<float> >& data, size_t sizeX, size_t sizeY, size_t sizeZ, hoNDArray< std::complex<float> >& dataPadded);
template EXPORTMRICORE void zeropad3D(const hoNDArray< std::complex<double> >& data, size_t sizeX, size_t sizeY, size_t sizeZ, hoNDArray< std::complex<double> >& dataPadded);

// ------------------------------------------------------------------------

template <typename T> 
void zeropad3DNoPresetZeros(const hoNDArray<T>& data, size_t sizeX, size_t sizeY, size_t sizeZ, hoNDArray<T>& dataPadded)
{
    try
    {
        size_t RO = data.get_size(0);
        size_t E1 = data.get_size(1);
        size_t E2 = data.get_size(2);
        size_t srcCHA = data.get_size(3);
        size_t dstCHA = data.get_size(4);

        GADGET_CHECK_THROW(sizeX>=RO);
        GADGET_CHECK_THROW(sizeY>=E1);
        GADGET_CHECK_THROW(sizeZ>=E2);

        if ( RO==sizeX && E1==sizeY && E2==sizeZ )
        {
            dataPadded = data;
            return;
        }

        size_t sRO, eRO, sE1, eE1, sE2, eE2;
        zpadRange(RO, sizeX, sRO, eRO);
        zpadRange(E1, sizeY, sE1, eE1);
        zpadRange(E2, sizeZ, sE2, eE2);

        GADGET_CHECK_THROW(dataPadded.get_size(0)==sizeX);
        GADGET_CHECK_THROW(dataPadded.get_size(1)==sizeY);
        GADGET_CHECK_THROW(dataPadded.get_size(2)==sizeZ);
        GADGET_CHECK_THROW(dataPadded.get_size(3)==srcCHA);
        GADGET_CHECK_THROW(dataPadded.get_size(4)==dstCHA);

        size_t num = data.get_number_of_elements()/(RO*E1*E2);

        long long n;

        const T* pData = data.begin();
        T* pDataPadded = dataPadded.begin();

        #pragma omp parallel for default(none) private(n) shared(num, sE2, eE2, sE1, eE1, sRO, RO, E1, E2, pData, pDataPadded, sizeX, sizeY, sizeZ)
        for ( n=0; n<(long long)num; n++ )
        {
            T* pDst = pDataPadded+n*sizeX*sizeY*sizeZ;
            T* pSrc = const_cast<T*>(pData)+n*RO*E1*E2;

            long long z;
            //#pragma omp parallel for default(none) private(z) shared(pDst, pSrc, sE2, eE2, sE1, eE1, sRO, RO, E1, E2, sizeX, sizeY, sizeZ) num_threads(2)
            for ( z=(long long)sE2; z<=(long long)eE2; z++ )
            {
                long long o1 = z*sizeX*sizeY + sRO;
                long long o2 = (z-sE2)*RO*E1;
                for ( size_t y=sE1; y<=eE1; y++ )
                {
                    memcpy(pDst+o1+y*sizeX, pSrc+o2+(y-sE1)*RO, sizeof(T)*RO);
                }
            }
        }
    }
    catch(...)
    {
        GADGET_THROW("Errors in zeropad3DNoPresetZeros(...) ... ");
    }

    return;
}

template EXPORTMRICORE void zeropad3DNoPresetZeros(const hoNDArray<float>& data, size_t sizeX, size_t sizeY, size_t sizeZ, hoNDArray<float>& dataPadded);
template EXPORTMRICORE void zeropad3DNoPresetZeros(const hoNDArray<double>& data, size_t sizeX, size_t sizeY, size_t sizeZ, hoNDArray<double>& dataPadded);
template EXPORTMRICORE void zeropad3DNoPresetZeros(const hoNDArray< std::complex<float> >& data, size_t sizeX, size_t sizeY, size_t sizeZ, hoNDArray< std::complex<float> >& dataPadded);
template EXPORTMRICORE void zeropad3DNoPresetZeros(const hoNDArray< std::complex<double> >& data, size_t sizeX, size_t sizeY, size_t sizeZ, hoNDArray< std::complex<double> >& dataPadded);

// ------------------------------------------------------------------------

template <typename T> 
void cutpad2D(const hoNDArray<T>& data, size_t sizeX, size_t sizeY, hoNDArray<T>& dataCut)
{
    try
    {
        size_t RO = data.get_size(0);
        size_t E1 = data.get_size(1);

        GADGET_CHECK_THROW(sizeX<=RO);
        GADGET_CHECK_THROW(sizeY<=E1);

        if ( RO==sizeX && E1==sizeY )
        {
            dataCut = data;
            return;
        }

        size_t sRO, eRO, sE1, eE1;
        zpadRange(sizeX, RO, sRO, eRO);
        zpadRange(sizeY, E1, sE1, eE1);

        boost::shared_ptr< std::vector<size_t> > dim = data.get_dimensions();
        (*dim)[0] = sizeX;
        (*dim)[1] = sizeY;
        dataCut.create(dim);

        size_t num = data.get_number_of_elements()/(RO*E1);

        long long n;

        const T* pData = data.begin();
        T* pDataCut = dataCut.begin();

        #pragma omp parallel for default(none) private(n) shared(num, sE1, eE1, sRO, RO, E1, pData, pDataCut, sizeX, sizeY)
        for ( n=0; n<(long long)num; n++ )
        {
            for ( size_t y=sE1; y<=eE1; y++ )
            {
                memcpy(pDataCut+n*sizeX*sizeY+(y-sE1)*sizeX, pData+n*RO*E1+y*RO+sRO, sizeof(T)*sizeX);
            }
        }
    }
    catch(...)
    {
        GADGET_THROW("Errors in cutpad2D(...) ... ");
    }
}

template EXPORTMRICORE void cutpad2D(const hoNDArray<float>& data, size_t sizeX, size_t sizeY, hoNDArray<float>& dataCut);
template EXPORTMRICORE void cutpad2D(const hoNDArray<double>& data, size_t sizeX, size_t sizeY, hoNDArray<double>& dataCut);
template EXPORTMRICORE void cutpad2D(const hoNDArray< std::complex<float> >& data, size_t sizeX, size_t sizeY, hoNDArray< std::complex<float> >& dataCut);
template EXPORTMRICORE void cutpad2D(const hoNDArray< std::complex<double> >& data, size_t sizeX, size_t sizeY, hoNDArray< std::complex<double> >& dataCut);

// ------------------------------------------------------------------------

template <typename T> 
void cutpad3D(const hoNDArray<T>& data, size_t sizeX, size_t sizeY, size_t sizeZ, hoNDArray<T>& dataCut)
{
    try
    {
        size_t RO = data.get_size(0);
        size_t E1 = data.get_size(1);
        size_t E2 = data.get_size(2);

        GADGET_CHECK_THROW(sizeX<=RO);
        GADGET_CHECK_THROW(sizeY<=E1);
        GADGET_CHECK_THROW(sizeZ<=E2);

        if ( RO==sizeX && E1==sizeY && E2==sizeZ )
        {
            dataCut = data;
            return;
        }

        size_t sRO, eRO, sE1, eE1, sE2, eE2;
        zpadRange(sizeX, RO, sRO, eRO);
        zpadRange(sizeY, E1, sE1, eE1);
        zpadRange(sizeZ, E2, sE2, eE2);

        boost::shared_ptr< std::vector<size_t> > dim = data.get_dimensions();
        (*dim)[0] = sizeX;
        (*dim)[1] = sizeY;
        (*dim)[2] = sizeZ;
        dataCut.create(dim);

        size_t num = data.get_number_of_elements()/(RO*E1*E2);

        long long n;

        const T* pData = data.begin();
        T* pDataCut = dataCut.begin();

        #pragma omp parallel for default(none) private(n) shared(num, sE2, eE2, sE1, eE1, sRO, RO, E1, E2, pData, pDataCut, sizeX, sizeY, sizeZ)
        for ( n=0; n<(long long)num; n++ )
        {
            for ( size_t z=sE2; z<=eE2; z++ )
            {
                for ( size_t y=sE1; y<=eE1; y++ )
                {
                    memcpy(pDataCut+n*sizeX*sizeY*sizeZ+(z-sE2)*sizeX*sizeY+(y-sE1)*sizeX, pData+n*RO*E1*E2+z*RO*E1+y*RO+sRO, sizeof(T)*sizeX);
                }
            }
        }
    }
    catch(...)
    {
        GADGET_THROW("Errors in cutpad3D(...) ... ");
    }
}

template EXPORTMRICORE void cutpad3D(const hoNDArray<float>& data, size_t sizeX, size_t sizeY, size_t sizeZ, hoNDArray<float>& dataCut);
template EXPORTMRICORE void cutpad3D(const hoNDArray<double>& data, size_t sizeX, size_t sizeY, size_t sizeZ, hoNDArray<double>& dataCut);
template EXPORTMRICORE void cutpad3D(const hoNDArray< std::complex<float> >& data, size_t sizeX, size_t sizeY, size_t sizeZ, hoNDArray< std::complex<float> >& dataCut);
template EXPORTMRICORE void cutpad3D(const hoNDArray< std::complex<double> >& data, size_t sizeX, size_t sizeY, size_t sizeZ, hoNDArray< std::complex<double> >& dataCut);

// ------------------------------------------------------------------------

template<typename T> 
bool sumOver1stDimension(const hoNDArray<T>& x, hoNDArray<T>& r)
{
    try
    {
        size_t RO = x.get_size(0);
        size_t num = x.get_number_of_elements()/(RO);

        boost::shared_ptr< std::vector<size_t> > dim = x.get_dimensions();

        std::vector<size_t> dimAve(*dim);
        dimAve[0] = 1;
        r.create(&dimAve);

        const T* pX = x.begin();
        T* pR = r.begin();

        int n;
        #pragma omp parallel for default(none) private(n) shared(RO, num, pX, pR)
        for ( n=0; n<(int)num; n++ )
        {
            T xsum = pX[n*RO];
            for (size_t ro=1; ro<RO; ro++ )
            {
                xsum += pX[n*RO+ro];
            }

            pR[n] = xsum;
        }
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in sumOver1stDimension(...) ... ");
        return false;
    }

    return true;
}

template<typename T> 
bool sumOver2ndDimension(const hoNDArray<T>& x, hoNDArray<T>& r)
{
    try
    {
        size_t NDim = x.get_number_of_dimensions();

        if ( NDim < 2 )
        {
            r = x;
            return true;
        }

        size_t RO = x.get_size(0);
        size_t E1 = x.get_size(1);

        size_t num = x.get_number_of_elements()/(RO*E1);

        boost::shared_ptr< std::vector<size_t> > dim = x.get_dimensions();

        std::vector<size_t> dimAve(*dim);
        dimAve[1] = 1;
        r.create(&dimAve);

        int n;
        #ifdef GCC_OLD_FLAG
            #pragma omp parallel for default(none) private(n) shared(RO, E1, num)
        #else
            #pragma omp parallel for default(none) private(n) shared(RO, E1, num, x, r)
        #endif
        for ( n=0; n<(int)num; n++ )
        {
            hoNDArray<T> xsum(RO, const_cast<T*>(r.begin()+n*RO));
            memcpy(xsum.begin(), x.begin()+n*RO*E1, xsum.get_number_of_bytes());

            for (size_t e1=1; e1<E1; e1++ )
            {
                hoNDArray<T> x1D(RO, const_cast<T*>(x.begin()+n*RO*E1+e1*RO));
                Gadgetron::add(x1D, xsum, xsum);
            }
        }
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in sumOver2ndDimension(...) ... ");
        return false;
    }

    return true;
}

template<typename T> 
bool sumOver3rdDimension(const hoNDArray<T>& x, hoNDArray<T>& r)
{
    try
    {
        size_t NDim = x.get_number_of_dimensions();

        if ( NDim < 3 )
        {
            r = x;
            return true;
        }

        size_t RO = x.get_size(0);
        size_t E1 = x.get_size(1);
        size_t CHA = x.get_size(2);

        size_t num = x.get_number_of_elements()/(RO*E1*CHA);

        boost::shared_ptr< std::vector<size_t> > dim = x.get_dimensions();

        std::vector<size_t> dimAve(*dim);
        dimAve[2] = 1;
        r.create(&dimAve);

        int n;
        #ifdef GCC_OLD_FLAG
            #pragma omp parallel default(none) private(n) shared(RO, E1, CHA, num) if (num>1)
        #else
            #pragma omp parallel default(none) private(n) shared(RO, E1, CHA, num, x, r) if (num>1)
        #endif
        {
            hoNDArray<T> xsum;
            hoNDArray<T> x2D;

            #pragma omp for
            for ( n=0; n<(int)num; n++ )
            {
                xsum.create(RO, E1, const_cast<T*>(r.begin()+n*RO*E1));
                memcpy(xsum.begin(), x.begin()+n*RO*E1*CHA, xsum.get_number_of_bytes());

                for (size_t cha=1; cha<CHA; cha++ )
                {
                    x2D.create(RO, E1, const_cast<T*>(x.begin()+n*RO*E1*CHA+cha*RO*E1));
                    Gadgetron::add(x2D, xsum, xsum);
                }
            }
        }
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in sumOver3rdDimension(...) ... ");
        return false;
    }

    return true;
}

template<typename T> bool sumOver4thDimension(const hoNDArray<T>& x, hoNDArray<T>& r)
{
    try
    {
        size_t NDim = x.get_number_of_dimensions();

        if ( NDim < 4 )
        {
            r = x;
            return true;
        }

        size_t RO = x.get_size(0);
        size_t E1 = x.get_size(1);
        size_t CHA = x.get_size(2);
        size_t N = x.get_size(3);

        size_t num = x.get_number_of_elements()/(RO*E1*CHA*N);

        boost::shared_ptr< std::vector<size_t> > dim = x.get_dimensions();

        std::vector<size_t> dimAve(*dim);
        dimAve[3] = 1;
        r.create(&dimAve);

        int n;
        #ifdef GCC_OLD_FLAG
            #pragma omp parallel for default(none) private(n) shared(RO, E1, CHA, N, num)
        #else
            #pragma omp parallel for default(none) private(n) shared(RO, E1, CHA, N, num, x, r)
        #endif
        for ( n=0; n<(int)num; n++ )
        {
            hoNDArray<T> xsum(RO, E1, CHA, const_cast<T*>(r.begin()+n*RO*E1*CHA));
            memcpy(xsum.begin(), x.begin()+n*RO*E1*CHA*N, xsum.get_number_of_bytes());

            for (size_t nn=1; nn<N; nn++ )
            {
                hoNDArray<T> x3D(RO, E1, CHA, const_cast<T*>(x.begin()+n*RO*E1*CHA*N+nn*RO*E1*CHA));
                Gadgetron::add(x3D, xsum, xsum);
            }
        }
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in sumOver4thDimension(const hoNDArray<T>& x, hoNDArray<T>& r) ... ");
        return false;
    }

    return true;
}

template<typename T> bool sumOver5thDimension(const hoNDArray<T>& x, hoNDArray<T>& r)
{
    try
    {
        size_t NDim = x.get_number_of_dimensions();

        if ( NDim < 5 )
        {
            r = x;
            return true;
        }

        size_t RO = x.get_size(0);
        size_t E1 = x.get_size(1);
        size_t CHA = x.get_size(2);
        size_t N = x.get_size(3);
        size_t S = x.get_size(4);

        size_t num = x.get_number_of_elements()/(RO*E1*CHA*N*S);

        boost::shared_ptr< std::vector<size_t> > dim = x.get_dimensions();

        std::vector<size_t> dimAve(*dim);
        dimAve[4] = 1;
        r.create(&dimAve);

        int n;
        #ifdef GCC_OLD_FLAG
            #pragma omp parallel for default(none) private(n) shared(RO, E1, CHA, N, S, num) if (num > 4)
        #else
            #pragma omp parallel for default(none) private(n) shared(RO, E1, CHA, N, S, num, x, r) if (num > 4)
        #endif
        for ( n=0; n<(int)num; n++ )
        {
            hoNDArray<T> xsum(RO, E1, CHA, N, const_cast<T*>(r.begin()+n*RO*E1*CHA*N));
            memcpy(xsum.begin(), x.begin()+n*RO*E1*CHA*N*S, xsum.get_number_of_bytes());

            for (size_t s=1; s<S; s++ )
            {
                hoNDArray<T> x4D(RO, E1, CHA, N, const_cast<T*>(x.begin()+n*RO*E1*CHA*N*S+s*RO*E1*CHA*N));
                Gadgetron::add(x4D, xsum, xsum);
            }
        }
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in sumOver5thDimension(const hoNDArray<T>& x, hoNDArray<T>& r) ... ");
        return false;
    }

    return true;
}

template EXPORTMRICORE bool sumOver1stDimension(const hoNDArray<float>& x, hoNDArray<float>& r);
template EXPORTMRICORE bool sumOver1stDimension(const hoNDArray<double>& x, hoNDArray<double>& r);
template EXPORTMRICORE bool sumOver1stDimension(const hoNDArray< std::complex<float> >& x, hoNDArray< std::complex<float> >& r);
template EXPORTMRICORE bool sumOver1stDimension(const hoNDArray< std::complex<double> >& x, hoNDArray< std::complex<double> >& r);

template EXPORTMRICORE bool sumOver2ndDimension(const hoNDArray<float>& x, hoNDArray<float>& r);
template EXPORTMRICORE bool sumOver2ndDimension(const hoNDArray<double>& x, hoNDArray<double>& r);
template EXPORTMRICORE bool sumOver2ndDimension(const hoNDArray< std::complex<float> >& x, hoNDArray< std::complex<float> >& r);
template EXPORTMRICORE bool sumOver2ndDimension(const hoNDArray< std::complex<double> >& x, hoNDArray< std::complex<double> >& r);

template EXPORTMRICORE bool sumOver3rdDimension(const hoNDArray<float>& x, hoNDArray<float>& r);
template EXPORTMRICORE bool sumOver3rdDimension(const hoNDArray<double>& x, hoNDArray<double>& r);
template EXPORTMRICORE bool sumOver3rdDimension(const hoNDArray< std::complex<float> >& x, hoNDArray< std::complex<float> >& r);
template EXPORTMRICORE bool sumOver3rdDimension(const hoNDArray< std::complex<double> >& x, hoNDArray< std::complex<double> >& r);

template EXPORTMRICORE bool sumOver4thDimension(const hoNDArray<float>& x, hoNDArray<float>& r);
template EXPORTMRICORE bool sumOver4thDimension(const hoNDArray<double>& x, hoNDArray<double>& r);
template EXPORTMRICORE bool sumOver4thDimension(const hoNDArray< std::complex<float> >& x, hoNDArray< std::complex<float> >& r);
template EXPORTMRICORE bool sumOver4thDimension(const hoNDArray< std::complex<double> >& x, hoNDArray< std::complex<double> >& r);

template EXPORTMRICORE bool sumOver5thDimension(const hoNDArray<float>& x, hoNDArray<float>& r);
template EXPORTMRICORE bool sumOver5thDimension(const hoNDArray<double>& x, hoNDArray<double>& r);
template EXPORTMRICORE bool sumOver5thDimension(const hoNDArray< std::complex<float> >& x, hoNDArray< std::complex<float> >& r);
template EXPORTMRICORE bool sumOver5thDimension(const hoNDArray< std::complex<double> >& x, hoNDArray< std::complex<double> >& r);

// ------------------------------------------------------------------------

template<typename T> 
bool permuteE2To3rdDimension(const hoNDArray<T>& x, hoNDArray<T>& r)
{
    try
    {
        boost::shared_ptr< std::vector<size_t> > dimX = x.get_dimensions();

        size_t NDim = dimX->size();

        if ( NDim <= 5 )
        {
            r = x;
            return true;
        }

        size_t RO = x.get_size(0);
        size_t E1 = x.get_size(1);
        size_t CHA = x.get_size(2);
        size_t SLC = x.get_size(3);
        size_t E2 = x.get_size(4);

        std::vector<size_t> dimR(*dimX);
        dimR[2] = E2;
        dimR[3] = CHA;
        dimR[4] = SLC;

        r.create(&dimR);

        size_t N2D = RO*E1;
        size_t N5D = RO*E1*CHA*E2*SLC;

        size_t N = x.get_number_of_elements()/N5D;

        const T* pX = x.begin();
        T* pR = r.begin();

        size_t n;
        for ( n=0; n<N; n++ )
        {
            int e2;
            #pragma omp parallel for default(none) private(e2) shared(N5D, N2D, pX, pR, CHA, SLC, E2, n)
            for ( e2=0; e2<E2; e2++ )
            {
                for ( size_t slc=0; slc<SLC; slc++ )
                {
                    for ( size_t cha=0; cha<CHA; cha++ )
                    {
                        memcpy(pR+n*N5D+slc*CHA*E2*N2D+cha*E2*N2D+e2*N2D, pX+n*N5D+e2*SLC*CHA*N2D+slc*CHA*N2D+cha*N2D, sizeof(T)*N2D);
                    }
                }
            }
        }
    }
    catch (...)
    {
        GADGET_ERROR_MSG("Errors in permuteE2To3rdDimension(const hoNDArray<T>& x, hoNDArray<T>& r) ... ");
        return false;
    }
    return true;
}

template<typename T> 
bool permuteE2To5thDimension(const hoNDArray<T>& x, hoNDArray<T>& r)
{
    try
    {
        boost::shared_ptr< std::vector<size_t> > dimX = x.get_dimensions();

        size_t NDim = dimX->size();

        if ( NDim < 5 )
        {
            r = x;
            return true;
        }

        size_t RO = x.get_size(0);
        size_t E1 = x.get_size(1);
        size_t E2 = x.get_size(2);
        size_t CHA = x.get_size(3);
        size_t SLC = x.get_size(4);

        std::vector<size_t> dimR(*dimX);
        dimR[2] = CHA;
        dimR[3] = SLC;
        dimR[4] = E2;

        r.create(&dimR);

        size_t N2D = RO*E1;
        size_t N5D = RO*E1*CHA*E2*SLC;

        size_t N = x.get_number_of_elements()/N5D;

        const T* pX = x.begin();
        T* pR = r.begin();

        size_t n;
        for ( n=0; n<N; n++ )
        {
            int e2;
            #pragma omp parallel for default(none) private(e2) shared(N5D, N2D, pX, pR, CHA, SLC, E2, n)
            for ( e2=0; e2<E2; e2++ )
            {
                for ( size_t slc=0; slc<SLC; slc++ )
                {
                    for ( size_t cha=0; cha<CHA; cha++ )
                    {
                        memcpy(pR+n*N5D+e2*SLC*CHA*N2D+slc*CHA*N2D+cha*N2D, pX+n*N5D+slc*CHA*E2*N2D+cha*E2*N2D+e2*N2D, sizeof(T)*N2D);
                    }
                }
            }
        }
    }
    catch (...)
    {
        GADGET_ERROR_MSG("Errors in permuteE2To5thDimension(const hoNDArray<T>& x, hoNDArray<T>& r) ... ");
        return false;
    }
    return true;
}

template<typename T> 
bool permuteROTo3rdDimensionFor3DRecon(const hoNDArray<T>& x, hoNDArray<T>& r)
{
    try
    {
        boost::shared_ptr< std::vector<size_t> > dimX = x.get_dimensions();

        size_t NDim = dimX->size();

        if ( NDim < 3 )
        {
            r = x;
            return true;
        }

        size_t RO = x.get_size(0);
        size_t E1 = x.get_size(1);
        size_t E2 = x.get_size(2);

        std::vector<size_t> dimR(*dimX);
        dimR[0] = E1;
        dimR[1] = E2;
        dimR[2] = RO;

        r.create(&dimR);

        size_t N3D = RO*E1*E2;

        size_t N = x.get_number_of_elements()/N3D;

        const T* pX = x.begin();
        T* pR = r.begin();

        long long n;

        #pragma omp parallel for default(none) private(n) shared(RO, E1, E2, N, pR, N3D, pX)
        for ( n=0; n<(long long)N; n++ )
        {
            T* pRn = pR + n*N3D;
            T* pXn = const_cast<T*>(pX) + n*N3D;

            for ( size_t e2=0; e2<E2; e2++ )
            {
                for ( size_t e1=0; e1<E1; e1++ )
                {
                    for ( size_t ro=0; ro<RO; ro++ )
                    {
                        pRn[e1+e2*E1+ro*E1*E2] = pXn[ro+e1*RO+e2*RO*E1];
                    }
                }
            }
        }
    }
    catch (...)
    {
        GADGET_ERROR_MSG("Errors in permuteROTo3rdDimensionFor3DRecon(const hoNDArray<T>& x, hoNDArray<T>& r) ... ");
        return false;
    }
    return true;
}

template<typename T> 
bool permuteROTo4thDimensionFor3DRecon(const hoNDArray<T>& x, hoNDArray<T>& r)
{
    try
    {
        boost::shared_ptr< std::vector<size_t> > dimX = x.get_dimensions();

        size_t NDim = dimX->size();

        if ( NDim < 4 )
        {
            r = x;
            return true;
        }

        size_t RO = x.get_size(0);
        size_t E1 = x.get_size(1);
        size_t E2 = x.get_size(2);
        size_t CHA = x.get_size(3);

        std::vector<size_t> dimR(*dimX);
        dimR[0] = E1;
        dimR[1] = E2;
        dimR[2] = CHA;
        dimR[3] = RO;

        r.create(&dimR);

        size_t N4D = RO*E1*E2*CHA;

        size_t N = x.get_number_of_elements()/N4D;

        const T* pX = x.begin();
        T* pR = r.begin();

        long long n;
        for ( n=0; n<(long long)N; n++ )
        {
            T* pRn = pR + n*N4D;
            T* pXn = const_cast<T*>(pX) + n*N4D;

            long long cha;

            #pragma omp parallel for default(none) private(cha) shared(RO, E1, E2, CHA, pXn, pRn)
            for ( cha=0; cha<(long long)CHA; cha++ )
            {
                for ( size_t e2=0; e2<E2; e2++ )
                {
                    for ( size_t e1=0; e1<E1; e1++ )
                    {
                        for ( size_t ro=0; ro<RO; ro++ )
                        {
                            pRn[e1+e2*E1+cha*E1*E2+ro*E1*E2*CHA] = pXn[ro+e1*RO+e2*RO*E1+cha*RO*E1*E2];
                        }
                    }
                }
            }
        }
    }
    catch (...)
    {
        GADGET_ERROR_MSG("Errors in permuteROTo4thDimensionFor3DRecon(const hoNDArray<T>& x, hoNDArray<T>& r) ... ");
        return false;
    }
    return true;
}

template<typename T> 
bool permuteROTo1stDimensionFor3DRecon(const hoNDArray<T>& x, hoNDArray<T>& r)
{
    try
    {
        boost::shared_ptr< std::vector<size_t> > dimX = x.get_dimensions();

        size_t NDim = dimX->size();

        if ( NDim < 4 )
        {
            r = x;
            return true;
        }

        size_t E1 = x.get_size(0);
        size_t E2 = x.get_size(1);
        size_t CHA = x.get_size(2);
        size_t RO = x.get_size(3);

        std::vector<size_t> dimR(*dimX);
        dimR[0] = RO;
        dimR[1] = E1;
        dimR[2] = E2;
        dimR[3] = CHA;

        r.create(&dimR);

        size_t N4D = RO*E1*E2*CHA;

        size_t N = x.get_number_of_elements()/N4D;

        const T* pX = x.begin();
        T* pR = r.begin();

        long long n;
        for ( n=0; n<(long long)N; n++ )
        {
            T* pRn = pR + n*N4D;
            T* pXn = const_cast<T*>(pX) + n*N4D;

            long long cha;

            #pragma omp parallel for default(none) private(cha) shared(RO, E1, E2, CHA, pXn, pRn)
            for ( cha=0; cha<(long long)CHA; cha++ )
            {
                for ( size_t e2=0; e2<E2; e2++ )
                {
                    for ( size_t e1=0; e1<E1; e1++ )
                    {
                        size_t indRn = e1*RO+e2*RO*E1+cha*RO*E1*E2;
                        size_t indXn = e1+e2*E1+cha*E1*E2;
                        for ( size_t ro=0; ro<RO; ro++ )
                        {
                            pRn[ro+indRn] = pXn[indXn+ro*E1*E2*CHA];
                        }
                    }
                }
            }
        }
    }
    catch (...)
    {
        GADGET_ERROR_MSG("Errors in permuteROTo1stDimensionFor3DRecon(const hoNDArray<T>& x, hoNDArray<T>& r) ... ");
        return false;
    }
    return true;
}

template<typename T> 
bool permute3rdDimensionTo1stDimension(const hoNDArray<T>& x, hoNDArray<T>& r)
{
    try
    {
        boost::shared_ptr< std::vector<size_t> > dimX = x.get_dimensions();

        size_t NDim = dimX->size();

        if ( NDim < 3 )
        {
            r = x;
            return true;
        }

        size_t RO = x.get_size(0);
        size_t E1 = x.get_size(1);
        size_t E2 = x.get_size(2);

        std::vector<size_t> dimR(*dimX);
        dimR[0] = E2;
        dimR[1] = RO;
        dimR[2] = E1;

        r.create(&dimR);

        size_t N3D = RO*E1*E2;

        size_t N = x.get_number_of_elements()/N3D;

        const T* pX = x.begin();
        T* pR = r.begin();

        long long n, e2;
        for ( n=0; n<(long long)N; n++ )
        {
            T* pRn = pR + n*N3D;
            T* pXn = const_cast<T*>(pX) + n*N3D;

            #pragma omp parallel for default(none) private(e2) shared(RO, E1, E2, pXn, pRn)
            for ( e2=0; e2<(long long)E2; e2++ )
            {
                for ( size_t e1=0; e1<E1; e1++ )
                {
                    size_t indRn = e2+e1*E2*RO;
                    size_t indXn = e1*RO+e2*RO*E1;
                    for ( size_t ro=0; ro<RO; ro++ )
                    {
                        pRn[ro*E2+indRn] = pXn[ro+indXn];
                    }
                }
            }
        }
    }
    catch (...)
    {
        GADGET_ERROR_MSG("Errors in permute3rdDimensionTo1stDimension(const hoNDArray<T>& x, hoNDArray<T>& r) ... ");
        return false;
    }
    return true;
}

template<typename T> 
bool permuteROTo5thDimensionFor3DRecon(const hoNDArray<T>& x, hoNDArray<T>& r)
{
    try
    {
        boost::shared_ptr< std::vector<size_t> > dimX = x.get_dimensions();

        size_t NDim = dimX->size();

        if ( NDim < 5 )
        {
            r = x;
            return true;
        }

        size_t RO = x.get_size(0);
        size_t E1 = x.get_size(1);
        size_t E2 = x.get_size(2);
        size_t srcCHA = x.get_size(3);
        size_t dstCHA = x.get_size(4);

        std::vector<size_t> dimR(*dimX);
        dimR[0] = E1;
        dimR[1] = E2;
        dimR[2] = srcCHA;
        dimR[3] = dstCHA;
        dimR[4] = RO;

        r.create(&dimR);

        size_t N5D = RO*E1*E2*srcCHA*dstCHA;

        size_t N = x.get_number_of_elements()/N5D;

        const T* pX = x.begin();
        T* pR = r.begin();

        long long n;
        for ( n=0; n<(long long)N; n++ )
        {
            T* pRn = pR + n*N5D;
            T* pXn = const_cast<T*>(pX) + n*N5D;

            long long dcha;

            #pragma omp parallel for default(none) private(dcha) shared(RO, E1, E2, srcCHA, dstCHA, pXn, pRn)
            for ( dcha=0; dcha<(long long)dstCHA; dcha++ )
            {
                for ( size_t scha=0; scha<(int)srcCHA; scha++ )
                {
                    for ( size_t e2=0; e2<E2; e2++ )
                    {
                        for ( size_t e1=0; e1<E1; e1++ )
                        {
                            size_t indRn = e1+e2*E1+scha*E1*E2+dcha*E1*E2*srcCHA;
                            size_t indXn = e1*RO+e2*RO*E1+scha*RO*E1*E2+dcha*RO*E1*E2*srcCHA;
                            for ( size_t ro=0; ro<RO; ro++ )
                            {
                                pRn[indRn+ro*E1*E2*srcCHA*dstCHA] = pXn[ro+indXn];
                            }
                        }
                    }
                }
            }
        }
    }
    catch (...)
    {
        GADGET_ERROR_MSG("Errors in permuteROTo5thDimensionFor3DRecon(const hoNDArray<T>& x, hoNDArray<T>& r) ... ");
        return false;
    }
    return true;
}

template EXPORTMRICORE bool permuteE2To3rdDimension(const hoNDArray<float>& x, hoNDArray<float>& r);
template EXPORTMRICORE bool permuteE2To3rdDimension(const hoNDArray<double>& x, hoNDArray<double>& r);
template EXPORTMRICORE bool permuteE2To3rdDimension(const hoNDArray< std::complex<float> >& x, hoNDArray< std::complex<float> >& r);
template EXPORTMRICORE bool permuteE2To3rdDimension(const hoNDArray< std::complex<double> >& x, hoNDArray< std::complex<double> >& r);

template EXPORTMRICORE bool permuteE2To5thDimension(const hoNDArray<float>& x, hoNDArray<float>& r);
template EXPORTMRICORE bool permuteE2To5thDimension(const hoNDArray<double>& x, hoNDArray<double>& r);
template EXPORTMRICORE bool permuteE2To5thDimension(const hoNDArray< std::complex<float> >& x, hoNDArray< std::complex<float> >& r);
template EXPORTMRICORE bool permuteE2To5thDimension(const hoNDArray< std::complex<double> >& x, hoNDArray< std::complex<double> >& r);

template EXPORTMRICORE bool permuteROTo3rdDimensionFor3DRecon(const hoNDArray<float>& x, hoNDArray<float>& r);
template EXPORTMRICORE bool permuteROTo3rdDimensionFor3DRecon(const hoNDArray<double>& x, hoNDArray<double>& r);
template EXPORTMRICORE bool permuteROTo3rdDimensionFor3DRecon(const hoNDArray< std::complex<float> >& x, hoNDArray< std::complex<float> >& r);
template EXPORTMRICORE bool permuteROTo3rdDimensionFor3DRecon(const hoNDArray< std::complex<double> >& x, hoNDArray< std::complex<double> >& r);

template EXPORTMRICORE bool permuteROTo4thDimensionFor3DRecon(const hoNDArray<float>& x, hoNDArray<float>& r);
template EXPORTMRICORE bool permuteROTo4thDimensionFor3DRecon(const hoNDArray<double>& x, hoNDArray<double>& r);
template EXPORTMRICORE bool permuteROTo4thDimensionFor3DRecon(const hoNDArray< std::complex<float> >& x, hoNDArray< std::complex<float> >& r);
template EXPORTMRICORE bool permuteROTo4thDimensionFor3DRecon(const hoNDArray< std::complex<double> >& x, hoNDArray< std::complex<double> >& r);

template EXPORTMRICORE bool permuteROTo1stDimensionFor3DRecon(const hoNDArray<float>& x, hoNDArray<float>& r);
template EXPORTMRICORE bool permuteROTo1stDimensionFor3DRecon(const hoNDArray<double>& x, hoNDArray<double>& r);
template EXPORTMRICORE bool permuteROTo1stDimensionFor3DRecon(const hoNDArray< std::complex<float> >& x, hoNDArray< std::complex<float> >& r);
template EXPORTMRICORE bool permuteROTo1stDimensionFor3DRecon(const hoNDArray< std::complex<double> >& x, hoNDArray< std::complex<double> >& r);

template EXPORTMRICORE bool permute3rdDimensionTo1stDimension(const hoNDArray<float>& x, hoNDArray<float>& r);
template EXPORTMRICORE bool permute3rdDimensionTo1stDimension(const hoNDArray<double>& x, hoNDArray<double>& r);
template EXPORTMRICORE bool permute3rdDimensionTo1stDimension(const hoNDArray< std::complex<float> >& x, hoNDArray< std::complex<float> >& r);
template EXPORTMRICORE bool permute3rdDimensionTo1stDimension(const hoNDArray< std::complex<double> >& x, hoNDArray< std::complex<double> >& r);

template EXPORTMRICORE bool permuteROTo5thDimensionFor3DRecon(const hoNDArray<float>& x, hoNDArray<float>& r);
template EXPORTMRICORE bool permuteROTo5thDimensionFor3DRecon(const hoNDArray<double>& x, hoNDArray<double>& r);
template EXPORTMRICORE bool permuteROTo5thDimensionFor3DRecon(const hoNDArray< std::complex<float> >& x, hoNDArray< std::complex<float> >& r);
template EXPORTMRICORE bool permuteROTo5thDimensionFor3DRecon(const hoNDArray< std::complex<double> >& x, hoNDArray< std::complex<double> >& r);

}
