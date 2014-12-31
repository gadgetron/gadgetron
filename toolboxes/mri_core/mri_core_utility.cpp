
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
void zeropad2D(const hoNDArray<T>& data, size_t sizeX, size_t sizeY, hoNDArray<T>& dataPadded, bool presetZeros)
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

        std::vector<size_t> dimPadded;
        data.get_dimensions(dimPadded);
        dimPadded[0] = sizeX;
        dimPadded[1] = sizeY;

        if ( !dataPadded.dimensions_equal(&dimPadded) )
        {
            dataPadded.create(dimPadded);
        }

        if ( presetZeros )
        {
            Gadgetron::clear(&dataPadded);
        }

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

template EXPORTMRICORE void zeropad2D(const hoNDArray<float>& data, size_t sizeX, size_t sizeY, hoNDArray<float>& dataPadded, bool presetZeros);
template EXPORTMRICORE void zeropad2D(const hoNDArray<double>& data, size_t sizeX, size_t sizeY, hoNDArray<double>& dataPadded, bool presetZeros);
template EXPORTMRICORE void zeropad2D(const hoNDArray< std::complex<float> >& data, size_t sizeX, size_t sizeY, hoNDArray< std::complex<float> >& dataPadded, bool presetZeros);
template EXPORTMRICORE void zeropad2D(const hoNDArray< std::complex<double> >& data, size_t sizeX, size_t sizeY, hoNDArray< std::complex<double> >& dataPadded, bool presetZeros);

// ------------------------------------------------------------------------

template <typename T> 
void zeropad3D(const hoNDArray<T>& data, size_t sizeX, size_t sizeY, size_t sizeZ, hoNDArray<T>& dataPadded, bool presetZeros)
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

        std::vector<size_t> dimPadded;
        data.get_dimensions(dimPadded);
        dimPadded[0] = sizeX;
        dimPadded[1] = sizeY;
        dimPadded[2] = sizeZ;

        if ( !dataPadded.dimensions_equal(&dimPadded) )
        {
            dataPadded.create(dimPadded);
        }

        if ( presetZeros )
        {
            Gadgetron::clear(&dataPadded);
        }

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

template EXPORTMRICORE void zeropad3D(const hoNDArray<float>& data, size_t sizeX, size_t sizeY, size_t sizeZ, hoNDArray<float>& dataPadded, bool presetZeros);
template EXPORTMRICORE void zeropad3D(const hoNDArray<double>& data, size_t sizeX, size_t sizeY, size_t sizeZ, hoNDArray<double>& dataPadded, bool presetZeros);
template EXPORTMRICORE void zeropad3D(const hoNDArray< std::complex<float> >& data, size_t sizeX, size_t sizeY, size_t sizeZ, hoNDArray< std::complex<float> >& dataPadded, bool presetZeros);
template EXPORTMRICORE void zeropad3D(const hoNDArray< std::complex<double> >& data, size_t sizeX, size_t sizeY, size_t sizeZ, hoNDArray< std::complex<double> >& dataPadded, bool presetZeros);

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
bool sumOverLastDimension(const hoNDArray<T>& x, hoNDArray<T>& r)
{
    try
    {
        boost::shared_ptr< std::vector<size_t> > dim = x.get_dimensions();
        size_t NDim = dim->size();

        std::vector<size_t> dimR(NDim-1);

        size_t d;
        for ( d=0; d<NDim-1; d++ )
        {
            dimR[d] = (*dim)[d];
        }

        if ( !r.dimensions_equal(&dimR) )
        {
            r.create(&dimR);
        }

        // Gadgetron::clear(&r);

        if ( x.get_size(NDim-1) <= 1 )
        {
            memcpy(r.begin(), x.begin(), x.get_number_of_bytes());
            return true;
        }

        size_t lastDim = x.get_size(NDim-1);
        size_t NR = r.get_number_of_elements();
        T* pA = const_cast<T*>(x.begin());
        T* pR = r.begin();

        memcpy(pR, pA, sizeof(T)*NR);

        // sum over the last dim
        hoNDArray<T> tmp;
        for ( d=1; d<lastDim; d++ )
        {
            tmp.create(&dimR, pA+d*NR);
            add(tmp, r, r);
        }
    }
    catch (...)
    {
        GADGET_ERROR_MSG("Errors in sumOverLastDimension(const hoNDArray<T>& x, hoNDArray<T>& r) ... ");
        return false;
    }
    return true;
}

template EXPORTMRICORE bool sumOverLastDimension(const hoNDArray<float>& x, hoNDArray<float>& r);
template EXPORTMRICORE bool sumOverLastDimension(const hoNDArray<double>& x, hoNDArray<double>& r);
template EXPORTMRICORE bool sumOverLastDimension(const hoNDArray< std::complex<float> >& x, hoNDArray< std::complex<float> >& r);
template EXPORTMRICORE bool sumOverLastDimension(const hoNDArray< std::complex<double> >& x, hoNDArray< std::complex<double> >& r);

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

inline void multiplyCplx(size_t N, const  std::complex<float> * x, const  std::complex<float> * y,  std::complex<float> * r)
{
    long long n;
    #pragma omp parallel for default(none) private(n) shared(N, x, y, r) if (N>64*1024)
    for (n = 0; n < (long long)N; n++)
    {
        const std::complex<float>& a1 = x[n];
        const std::complex<float>& b1 = y[n];
        const float a = a1.real();
        const float b = a1.imag();
        const float c = b1.real();
        const float d = b1.imag();

        reinterpret_cast<float(&)[2]>(r[n])[0] = a*c-b*d;
        reinterpret_cast<float(&)[2]>(r[n])[1] = a*d+b*c;
    }
}

inline void multiplyCplx(size_t N, const  std::complex<double> * x, const  std::complex<double> * y,  std::complex<double> * r)
{
    long long n;
    #pragma omp parallel for default(none) private(n) shared(N, x, y, r) if (N>64*1024)
    for (n = 0; n < (long long)N; n++)
    {
        const std::complex<double>& a1 = x[n];
        const std::complex<double>& b1 = y[n];
        const double a = a1.real();
        const double b = a1.imag();
        const double c = b1.real();
        const double d = b1.imag();

        reinterpret_cast<double(&)[2]>(r[n])[0] = a*c-b*d;
        reinterpret_cast<double(&)[2]>(r[n])[1] = a*d+b*c;
    }
}

template <typename T> 
bool multipleMultiply(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& r)
{
    GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()<=y.get_number_of_elements());
    if ( r.get_number_of_elements()!=y.get_number_of_elements())
    {
        r = y;
    }

    long long Nx = x.get_number_of_elements();
    long long N = y.get_number_of_elements() / Nx;

    const T* pX = x.begin();

    long long n;

    if ( typeid(T)==typeid(float) )
    {
        #pragma omp parallel for default(none) private(n) shared(pX, y, r, Nx, N)
        for ( n=0; n<N; n++ )
        {
            const T* pY = y.begin()+n*Nx;
            T* pR = r.begin() + n*Nx;

            long long ii;
            for ( ii=0; ii<Nx; ii++ )
            {
                pR[ii] = pX[ii] * pY[ii];
            }
        }
    }
    else if ( typeid(T)==typeid(double) )
    {
        #pragma omp parallel for default(none) private(n) shared(pX, y, r, Nx, N)
        for ( n=0; n<N; n++ )
        {
            const T* pY = y.begin()+n*Nx;
            T* pR = r.begin() + n*Nx;

            long long ii;
            for ( ii=0; ii<Nx; ii++ )
            {
                pR[ii] = pX[ii] * pY[ii];
            }
        }
    }
    else if ( typeid(T)==typeid( std::complex<float> ) )
    {
        #ifdef GCC_OLD_FLAG
            #pragma omp parallel for default(none) private(n) shared(Nx, N)
        #else
            #pragma omp parallel for default(none) private(n) shared(x, y, r, Nx, N)
        #endif
        for ( n=0; n<N; n++ )
        {
            multiplyCplx(x.get_number_of_elements(), (const std::complex<float>*)(x.begin()), (const std::complex<float>*)(y.begin()+n*Nx), (std::complex<float>*)(r.begin()+n*Nx));
        }
    }
    else if ( typeid(T)==typeid( std::complex<double> ) )
    {
        #ifdef GCC_OLD_FLAG
            #pragma omp parallel for default(none) private(n) shared(Nx, N)
        #else
            #pragma omp parallel for default(none) private(n) shared(x, y, r, Nx, N)
        #endif
        for ( n=0; n<N; n++ )
        {
            multiplyCplx(x.get_number_of_elements(), (const std::complex<double>*)(x.begin()), (const std::complex<double>*)(y.begin()+n*Nx), (std::complex<double>*)(r.begin()+n*Nx));
        }
    }
    else
    {
        GADGET_ERROR_MSG("multipleMultiply : unsupported type " << typeid(T).name());
        return false;
    }

    return true;
}

template EXPORTMRICORE bool multipleMultiply(const hoNDArray<float>& x, const hoNDArray<float>& y, hoNDArray<float>& r);
template EXPORTMRICORE bool multipleMultiply(const hoNDArray<double>& x, const hoNDArray<double>& y, hoNDArray<double>& r);
template EXPORTMRICORE bool multipleMultiply(const hoNDArray< std::complex<float> >& x, const hoNDArray< std::complex<float> >& y, hoNDArray< std::complex<float> >& r);
template EXPORTMRICORE bool multipleMultiply(const hoNDArray< std::complex<double> >& x, const hoNDArray< std::complex<double> >& y, hoNDArray< std::complex<double> >& r);

}
