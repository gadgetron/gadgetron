
/** \file   mri_core_utility.cpp
    \brief  Implementation useful utility functionalities for 2D and 3D MRI parallel imaging
    \author Hui Xue
*/

#include "mri_core_utility.h"
#include "hoNDArray_elemwise.h"

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

template<typename T> 
void cutpad(const hoNDArray<T>& data, const std::vector<size_t>& cutSize, hoNDArray<T>& dataCut)
{
    try
    {
        size_t NDim = data.get_number_of_dimensions();
        size_t NDimCut = cutSize.size();

        GADGET_CHECK_THROW(NDim >= NDimCut);

        bool needToCut = false;

        size_t n;
        for (n = 0; n < NDimCut; n++)
        {
            GADGET_CHECK_THROW(cutSize[n] <= data.get_size(n));

            if (cutSize[n]<data.get_size(n))
            {
                needToCut = true;
            }
        }

        if (!needToCut)
        {
            dataCut = data;
            return;
        }

        std::vector<size_t> dims;
        data.get_dimensions(dims);

        std::vector<size_t> dimsCut;
        dimsCut = dims;

        for (n = 0; n < NDimCut; n++)
        {
            dimsCut[n] = cutSize[n];
        }

        if (!dataCut.dimensions_equal(&dimsCut))
        {
            dataCut.create(dimsCut);
        }

        std::vector<size_t> start(NDimCut), end(NDimCut);
        for (n = 0; n < NDimCut; n++)
        {
            GADGET_CATCH_THROW(zpadRange(dimsCut[n], dims[n], start[n], end[n]));
        }

        size_t len = dimsCut[0];
        size_t num = dataCut.get_number_of_elements() / len;

        long long k;

        const T* pData = data.begin();
        T* pCutData = dataCut.begin();

        #pragma omp parallel default(none) private(k, n) shared(pData, pCutData, num, NDimCut, len, data, start, dataCut)
        {
            std::vector<size_t> ind;

            #pragma omp for 
            for (k = 0; k < (long long)num; k++)
            {
                T* pCutDataCur = pCutData + k*len;
                ind = dataCut.calculate_index(k*len);

                for (n = 0; n < NDimCut; n++)
                {
                    ind[n] += start[n];
                }

                size_t offset = data.calculate_offset(ind);
                const T* pDataCur = pData + offset;

                memcpy(pCutDataCur, pDataCur, sizeof(T)*len);
            }
        }
    }
    catch (...)
    {
        GADGET_THROW("Errors in cutpad(...) ... ");
    }
}

// ------------------------------------------------------------------------

template <typename T>
void cutpad2D(const hoNDArray<T>& data, size_t sizeX, size_t sizeY, hoNDArray<T>& dataCut)
{
    std::vector<size_t> cutSize(2);
    cutSize[0] = sizeX;
    cutSize[1] = sizeY;

    return cutpad(data, cutSize, dataCut);
}

template EXPORTMRICORE void cutpad2D(const hoNDArray<float>& data, size_t sizeX, size_t sizeY, hoNDArray<float>& dataCut);
template EXPORTMRICORE void cutpad2D(const hoNDArray<double>& data, size_t sizeX, size_t sizeY, hoNDArray<double>& dataCut);
template EXPORTMRICORE void cutpad2D(const hoNDArray< std::complex<float> >& data, size_t sizeX, size_t sizeY, hoNDArray< std::complex<float> >& dataCut);
template EXPORTMRICORE void cutpad2D(const hoNDArray< std::complex<double> >& data, size_t sizeX, size_t sizeY, hoNDArray< std::complex<double> >& dataCut);

// ------------------------------------------------------------------------

template <typename T>
void cutpad3D(const hoNDArray<T>& data, size_t sizeX, size_t sizeY, size_t sizeZ, hoNDArray<T>& dataCut)
{
    std::vector<size_t> cutSize(3);
    cutSize[0] = sizeX;
    cutSize[1] = sizeY;
    cutSize[2] = sizeZ;

    return cutpad(data, cutSize, dataCut);
}

template EXPORTMRICORE void cutpad3D(const hoNDArray<float>& data, size_t sizeX, size_t sizeY, size_t sizeZ, hoNDArray<float>& dataCut);
template EXPORTMRICORE void cutpad3D(const hoNDArray<double>& data, size_t sizeX, size_t sizeY, size_t sizeZ, hoNDArray<double>& dataCut);
template EXPORTMRICORE void cutpad3D(const hoNDArray< std::complex<float> >& data, size_t sizeX, size_t sizeY, size_t sizeZ, hoNDArray< std::complex<float> >& dataCut);
template EXPORTMRICORE void cutpad3D(const hoNDArray< std::complex<double> >& data, size_t sizeX, size_t sizeY, size_t sizeZ, hoNDArray< std::complex<double> >& dataCut);
}
