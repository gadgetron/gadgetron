
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
}
