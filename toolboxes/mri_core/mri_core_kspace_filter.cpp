
/** \file   mri_core_kspace_filter.cpp
    \brief  Implementation kspace filter functionalities for 2D and 3D MRI parallel imaging
    \author Hui Xue
*/

#include "mri_core_kspace_filter.h"
#include "hoNDArray_elemwise.h"

#ifdef __clang__
    #define unary_function  __unary_function
#endif

#include <boost/algorithm/string.hpp>

#ifdef M_PI
    #undef M_PI
#endif // M_PI
#define M_PI 3.14159265358979323846

namespace Gadgetron
{

ISMRMRDKSPACEFILTER get_kspace_filter_type(const std::string& name)
{
    std::string name_lower(name);
    boost::algorithm::to_lower(name_lower);

    if (name_lower == "gaussian")
    {
        return ISMRMRD_FILTER_GAUSSIAN;
    }
    else if (name_lower == "hanning")
    {
        return ISMRMRD_FILTER_HANNING;
    }
    else if (name_lower == "taperedhanning")
    {
        return ISMRMRD_FILTER_TAPERED_HANNING;
    }
    else if (name_lower == "none")
    {
        return ISMRMRD_FILTER_NONE;
    }

    GERROR_STREAM("Unrecognized kspace filter name : " << name);

    return ISMRMRD_FILTER_NONE;
}

std::string get_kspace_filter_name(ISMRMRDKSPACEFILTER v)
{
    std::string name;

    switch (v)
    {
        case ISMRMRD_FILTER_GAUSSIAN:
        {
            name = "Gaussian";
            break;
        }

        case ISMRMRD_FILTER_HANNING:
        {
            name = "Hanning";
            break;
        }

        case ISMRMRD_FILTER_TAPERED_HANNING:
        {
            name = "TaperedHanning";
            break;
        }

        case ISMRMRD_FILTER_NONE:
        {
            name = "none";
            break;
        }

        default:
        {
            GERROR_STREAM("Unrecognized kspace filter type : " << v);
            name = "none"; 
        }
    }

    return name;
}

template<typename T>
void generate_symmetric_filter(size_t len, hoNDArray<T>& filter, ISMRMRDKSPACEFILTER filterType, double sigma, size_t width)
{
    try
    {
        if (len == 0) return;

        filter.create(len);

        if (width == 0 || width >= len) width = 1;

        size_t ii;
        if (filterType == ISMRMRD_FILTER_GAUSSIAN)
        {
            double r = -1.0*sigma*sigma / 2;

            if (len % 2 == 0)
            {
                // to make sure the zero points match and boundary of filters are symmetric
                double stepSize = 2.0 / (len - 2);
                std::vector<double> x(len - 1);

                for (ii = 0; ii<len - 1; ii++)
                {
                    x[ii] = -1 + ii*stepSize;
                }

                for (ii = 0; ii<len - 1; ii++)
                {
                    filter(ii + 1) = T(std::exp(r*(x[ii] * x[ii])));
                }

                filter(0) = T(0);
            }
            else
            {
                double stepSize = 2.0 / (len - 1);
                std::vector<double> x(len);

                for (ii = 0; ii<len; ii++)
                {
                    x[ii] = -1 + ii*stepSize;
                }

                for (ii = 0; ii<len; ii++)
                {
                    filter(ii) = T(std::exp(r*(x[ii] * x[ii])));
                }
            }
        }
        else if (filterType == ISMRMRD_FILTER_TAPERED_HANNING)
        {
            hoNDArray<T> w(width);

            for (ii = 1; ii <= width; ii++)
            {
                w(ii - 1) = T((0.5 * (1 - std::cos(2.0*M_PI*ii / (2 * width + 1)))));
            }

            // make sure the center of the filter will end up being 1:
            Gadgetron::fill(filter, T(1.0));
            
            if (len % 2 == 0)
            {
                for (ii = 1; ii <= width; ii++)
                {
                    filter(ii) = w(ii - 1);
                    filter(len - ii) = filter(ii);
                }

                filter(0) = T(0);
            }
            else
            {
                for (ii = 1; ii <= width; ii++)
                {
                    filter(ii - 1) = w(ii - 1);
                    filter(len - ii) = filter(ii - 1);
                }
            }
        }
        else if (filterType == ISMRMRD_FILTER_HANNING)
        {
            if (len % 2 == 0)
            {
                size_t N = len - 1;
                double halfLen = (double)((N + 1) / 2);
                for (ii = 1; ii <= halfLen; ii++)
                {
                    filter(ii) = T((0.5 * (1 - std::cos(2.0*M_PI*ii / (N + 1)))));
                }

                for (ii = (size_t)halfLen; ii<N; ii++)
                {
                    filter(ii + 1) = filter(N - ii);
                }

                filter(0) = T(0);
            }
            else
            {
                double halfLen = (double)((len + 1) / 2);
                for (ii = 1; ii <= (size_t)halfLen; ii++)
                {
                    filter(ii - 1) = T((0.5 * (1 - std::cos(2.0*M_PI*ii / (len + 1)))));
                }

                for (ii = (size_t)halfLen; ii<len; ii++)
                {
                    filter(ii) = filter(len - 1 - ii);
                }
            }
        }
        else if (filterType == ISMRMRD_FILTER_NONE)
        {
            Gadgetron::fill(filter, T(1.0));
        }
        else
        {
            GADGET_THROW("generate_symmetric_filter, unrecognized fiter type ... ");
        }

        T sos = 0.0f;
        for (ii = 0; ii<len; ii++)
        {
            sos += filter(ii)*filter(ii);
        }

        T r = T(1.0 / std::sqrt(std::abs(sos) / (len)));
        for (ii = 0; ii<len; ii++)
        {
            filter(ii) *= r;
        }
    }
    catch (...)
    {
        GADGET_THROW("Errors in generate_symmetric_filter(...) ... ");
    }
}

template EXPORTMRICORE void generate_symmetric_filter(size_t len, hoNDArray<float>& filter, ISMRMRDKSPACEFILTER filterType, double sigma, size_t width);
template EXPORTMRICORE void generate_symmetric_filter(size_t len, hoNDArray<double>& filter, ISMRMRDKSPACEFILTER filterType, double sigma, size_t width);
template EXPORTMRICORE void generate_symmetric_filter(size_t len, hoNDArray< std::complex<float> >& filter, ISMRMRDKSPACEFILTER filterType, double sigma, size_t width);
template EXPORTMRICORE void generate_symmetric_filter(size_t len, hoNDArray< std::complex<double> >& filter, ISMRMRDKSPACEFILTER filterType, double sigma, size_t width);

// ------------------------------------------------------------------------

template<typename T>
void generate_asymmetric_filter(size_t len, size_t start, size_t end, hoNDArray<T>& filter, ISMRMRDKSPACEFILTER filterType, size_t width, bool densityComp)
{
    try
    {
        if (len == 0) return;

        if (start > len - 1) start = 0;
        if (end > len - 1) end = len - 1;

        if (start > end)
        {
            start = 0;
            end = len - 1;
        }

        filter.create(len);
        Gadgetron::clear(filter);

        size_t ii;
        for (ii = start; ii <= end; ii++)
        {
            filter(ii) = T(1.0);
        }

        if (width == 0 || width >= len) width = 1;

        hoNDArray<T> w(width);

        if (filterType == ISMRMRD_FILTER_TAPERED_HANNING)
        {
            for (ii = 1; ii <= width; ii++)
            {
                w(ii - 1) = T((0.5 * (1 - std::cos(2.0*M_PI*ii / (2 * width + 1)))));
            }
        }
        else if (filterType == ISMRMRD_FILTER_NONE)
        {
            Gadgetron::fill(w, T(1.0));
        }
        else
        {
            GADGET_THROW("generate_symmetric_filter, unrecognized fiter type ... ");
        }

        if (densityComp)
        {
            size_t startSym(0), endSym(len - 1);
            find_symmetric_sampled_region(start, end, len / 2, startSym, endSym);

            if (start == 0 && end == len - 1)
            {
                for (ii = 1; ii <= width; ii++)
                {
                    filter(ii - 1) = w(ii - 1);
                    filter(len - ii) = filter(ii - 1);
                }
            }

            if (start == 0 && end<len - 1)
            {
                for (ii = 0; ii<startSym; ii++)
                {
                    filter(ii) = 2.0;
                }

                for (ii = 1; ii <= width; ii++)
                {
                    filter(ii - 1 + startSym) = T(1.0) + w(width - ii);
                    filter(end - ii + 1) = w(ii - 1);
                }
            }

            if (start>0 && end == len - 1)
            {
                for (ii = endSym + 1; ii<len; ii++)
                {
                    filter(ii) = 2.0;
                }

                for (ii = 1; ii <= width; ii++)
                {
                    filter(endSym - ii + 1) = T(1.0) + w(width - ii);
                    filter(start + ii - 1) = w(ii - 1);
                }
            }

            if (start>0 && end<len - 1)
            {
                if (start == startSym && end == endSym)
                {
                    for (ii = 1; ii <= width; ii++)
                    {
                        filter(start + ii - 1) = w(ii - 1);
                        filter(end - ii + 1) = w(ii - 1);
                    }
                }
                else if (start == startSym && end>endSym)
                {
                    for (ii = endSym + 1; ii <= end; ii++)
                    {
                        filter(ii) = 2.0;
                    }

                    for (ii = 1; ii <= width; ii++)
                    {
                        filter(end - ii + 1) = T(1.0) + w(ii - 1);
                        filter(endSym - ii + 1) = w(width - ii);
                        filter(start + ii - 1) = w(ii - 1);
                    }
                }
                else if (start<startSym && end == endSym)
                {
                    for (ii = start; ii<startSym; ii++)
                    {
                        filter(ii) = 2.0;
                    }

                    for (ii = 1; ii <= width; ii++)
                    {
                        filter(ii - 1 + start) = T(1.0) + w(ii - 1);
                        filter(ii - 1 + startSym) = w(width - ii);
                        filter(end - ii + 1) = w(ii - 1);
                    }
                }
                else
                {
                    for (ii = 1; ii <= width; ii++)
                    {
                        filter(start + ii - 1) = w(ii - 1);
                        filter(end - ii + 1) = w(ii - 1);
                    }
                }
            }
        }
        else
        {
            if (start == 0 && end == len - 1)
            {
                for (ii = 1; ii <= width; ii++)
                {
                    filter(ii - 1) = w(ii - 1);
                    filter(len - ii) = filter(ii - 1);
                }
            }

            if (start == 0 && end<len - 1)
            {
                for (ii = 1; ii <= width; ii++)
                {
                    filter(end - ii + 1) = w(ii - 1);
                }
            }

            if (start>0 && end == len - 1)
            {
                for (ii = 1; ii <= width; ii++)
                {
                    filter(start + ii - 1) = w(ii - 1);
                }
            }

            if (start>0 && end<len - 1)
            {
                for (ii = 1; ii <= width; ii++)
                {
                    filter(start + ii - 1) = w(ii - 1);
                    filter(end - ii + 1) = w(ii - 1);
                }
            }
        }

        T sos = 0.0f;
        for (ii = 0; ii<len; ii++)
        {
            sos += filter(ii)*filter(ii);
        }

        T r = (T)(1.0 / std::sqrt(std::abs(sos) / (end - start + 1))); // SNR unit filter
        for (ii = 0; ii<len; ii++)
        {
            filter(ii) *= r;
        }
    }
    catch (...)
    {
        GADGET_THROW("Errors in generate_asymmetric_filter(...) ... ");
    }
}

template EXPORTMRICORE void generate_asymmetric_filter(size_t len, size_t start, size_t end, hoNDArray<float>& filter, ISMRMRDKSPACEFILTER filterType, size_t width, bool densityComp);
template EXPORTMRICORE void generate_asymmetric_filter(size_t len, size_t start, size_t end, hoNDArray<double>& filter, ISMRMRDKSPACEFILTER filterType, size_t width, bool densityComp);
template EXPORTMRICORE void generate_asymmetric_filter(size_t len, size_t start, size_t end, hoNDArray< std::complex<float> >& filter, ISMRMRDKSPACEFILTER filterType, size_t width, bool densityComp);
template EXPORTMRICORE void generate_asymmetric_filter(size_t len, size_t start, size_t end, hoNDArray< std::complex<double> >& filter, ISMRMRDKSPACEFILTER filterType, size_t width, bool densityComp);

// ------------------------------------------------------------------------

template<typename T>
void generate_symmetric_filter_ref(size_t len, size_t start, size_t end, hoNDArray<T>& filter)
{
    try
    {
        GADGET_CHECK_THROW(len >= 2);
        GADGET_CHECK_THROW(start >= 0 && end <= len - 1 && start <= end);

        if (start == 0 && end == len - 1)
        {
            generate_symmetric_filter(len, filter, ISMRMRD_FILTER_HANNING);
            return;
        }

        size_t centerInd = len / 2;

        size_t lenFilter(0); // make a symmetric filter with zero at the center
        size_t lenFilterEnd = 2 * (end - centerInd) + 1;
        size_t lenFilterStart = 2 * (centerInd - start) + 1;

        if (start == 0 && end<len - 1)
        {
            lenFilter = lenFilterEnd;
        }
        else if (start>0 && end == len - 1)
        {
            lenFilter = lenFilterStart;
        }
        else if (start>0 && end<len - 1)
        {
            lenFilter = ((lenFilterStart<lenFilterEnd) ? lenFilterStart : lenFilterEnd);
        }
        else
        {
            GERROR_STREAM("generate_symmetric_filter_ref, invalid inputs : start - end - len ... " << start << " " << end << " " << len);
            GADGET_THROW("generate_symmetric_filter_ref, invalid inputs ... ");
        }

        GADGET_CHECK_THROW(lenFilter>0);

        hoNDArray<T> filterSym(lenFilter);
        generate_symmetric_filter(lenFilter, filterSym, ISMRMRD_FILTER_HANNING);

        filter.create(len);
        Gadgetron::clear(&filter);

        if (start == 0 && end<len - 1)
        {
            memcpy(filter.begin() + end - lenFilter + 1, filterSym.begin(), filterSym.get_number_of_bytes());
            return;
        }
        else if (start>0 && end == len - 1)
        {
            memcpy(filter.begin() + start, filterSym.begin(), filterSym.get_number_of_bytes());
            return;
        }
        else if (start>0 && end<len - 1)
        {
            if (lenFilter == lenFilterStart)
            {
                memcpy(filter.begin() + start, filterSym.begin(), filterSym.get_number_of_bytes());
            }
            else
            {
                memcpy(filter.begin() + end - lenFilter + 1, filterSym.begin(), filterSym.get_number_of_bytes());
            }

            return;
        }
        else
        {
            GERROR_STREAM("Invalid inputs : start - end - len : " << start << " " << end << " " << len);
            GADGET_THROW("generate_symmetric_filter_ref, invalid inputs : start - end - len");
        }
    }
    catch (...)
    {
        GADGET_THROW("Errors in generate_symmetric_filter_ref(...) ... ");
    }
}

template EXPORTMRICORE void generate_symmetric_filter_ref(size_t len, size_t start, size_t end, hoNDArray<float>& filter);
template EXPORTMRICORE void generate_symmetric_filter_ref(size_t len, size_t start, size_t end, hoNDArray<double>& filter);
template EXPORTMRICORE void generate_symmetric_filter_ref(size_t len, size_t start, size_t end, hoNDArray< std::complex<float> >& filter);
template EXPORTMRICORE void generate_symmetric_filter_ref(size_t len, size_t start, size_t end, hoNDArray< std::complex<double> >& filter);

// ------------------------------------------------------------------------

template <typename T> 
void compute_2d_filter(const hoNDArray<T>& fx, const hoNDArray<T>& fy, hoNDArray<T>& fxy)
{
    try
    {
        size_t RO = fx.get_size(0);
        size_t E1 = fy.get_size(0);

        fxy.create(RO, E1);
        T* pFxy = fxy.begin();

        size_t x, y;

        for (y = 0; y<E1; y++)
        {
            for (x = 0; x<RO; x++)
            {
                pFxy[y*RO + x] = fx(x) * fy(y);
            }
        }
    }
    catch (...)
    {
        GADGET_THROW("Errors in compute_2d_filter(...) ... ");
    }
}

template EXPORTMRICORE void compute_2d_filter(const hoNDArray<float>& fx, const hoNDArray<float>& fy, hoNDArray<float>& fxy);
template EXPORTMRICORE void compute_2d_filter(const hoNDArray<double>& fx, const hoNDArray<double>& fy, hoNDArray<double>& fxy);
template EXPORTMRICORE void compute_2d_filter(const hoNDArray< std::complex<float> >& fx, const hoNDArray< std::complex<float> >& fy, hoNDArray< std::complex<float> >& fxy);
template EXPORTMRICORE void compute_2d_filter(const hoNDArray< std::complex<double> >& fx, const hoNDArray< std::complex<double> >& fy, hoNDArray< std::complex<double> >& fxy);

// ------------------------------------------------------------------------

void compute_2d_filter(const hoNDArray<float>& fx, const hoNDArray<float>& fy, hoNDArray< std::complex<float> >& fxy)
{
    try
    {
        size_t RO = fx.get_size(0);
        size_t E1 = fy.get_size(0);

        fxy.create(RO, E1);
        std::complex<float> * pFxy = fxy.begin();

        size_t x, y;

        for (y = 0; y<E1; y++)
        {
            for (x = 0; x<RO; x++)
            {
                pFxy[y*RO + x] = std::complex<float>(fx(x) * fy(y));
            }
        }
    }
    catch (...)
    {
        GADGET_THROW("Errors in compute_2d_filter(float) ... ");
    }
}

// ------------------------------------------------------------------------

void compute_2d_filter(const hoNDArray<double>& fx, const hoNDArray<double>& fy, hoNDArray< std::complex<double> >& fxy)
{
    try
    {
        size_t RO = fx.get_size(0);
        size_t E1 = fy.get_size(0);

        fxy.create(RO, E1);
        std::complex<double> * pFxy = fxy.begin();

        size_t x, y;

        for (y = 0; y<E1; y++)
        {
            for (x = 0; x<RO; x++)
            {
                pFxy[y*RO + x] = std::complex<double>(fx(x) * fy(y));
            }
        }
    }
    catch (...)
    {
        GADGET_THROW("Errors in compute_2d_filter(double) ... ");
    }
}

// ------------------------------------------------------------------------

template <typename T> 
void compute_3d_filter(const hoNDArray<T>& fx, const hoNDArray<T>& fy, const hoNDArray<T>& fz, hoNDArray<T>& fxyz)
{
    try
    {
        size_t RO = fx.get_size(0);
        size_t E1 = fy.get_size(0);
        size_t E2 = fz.get_size(0);

        fxyz.create(RO, E1, E2);
        T* pFxyz = fxyz.begin();

        const T* px = fx.begin();
        const T* py = fy.begin();
        const T* pz = fz.begin();

        size_t x, y, z;

        T vz, vy, vx;

        for (z = 0; z<E2; z++)
        {
            vz = pz[z];
            for (y = 0; y<E1; y++)
            {
                vy = py[y];
                for (x = 0; x<RO; x++)
                {
                    vx = px[x];
                    pFxyz[x+y*RO+z*RO*E1] = (vx*vz*vy);
                }
            }
        }
    }
    catch (...)
    {
        GADGET_THROW("Errors in compute_3d_filter(...) ... ");
    }
}

template EXPORTMRICORE void compute_3d_filter(const hoNDArray<float>& fx, const hoNDArray<float>& fy, const hoNDArray<float>& fz, hoNDArray<float>& fxyz);
template EXPORTMRICORE void compute_3d_filter(const hoNDArray<double>& fx, const hoNDArray<double>& fy, const hoNDArray<double>& fz, hoNDArray<double>& fxyz);
template EXPORTMRICORE void compute_3d_filter(const hoNDArray< std::complex<float> >& fx, const hoNDArray< std::complex<float> >& fy, const hoNDArray< std::complex<float> >& fz, hoNDArray< std::complex<float> >& fxyz);
template EXPORTMRICORE void compute_3d_filter(const hoNDArray< std::complex<double> >& fx, const hoNDArray< std::complex<double> >& fy, const hoNDArray< std::complex<double> >& fz, hoNDArray< std::complex<double> >& fxyz);

// ------------------------------------------------------------------------

void compute_3d_filter(const hoNDArray<float>& fx, const hoNDArray<float>& fy, const hoNDArray<float>& fz, hoNDArray< std::complex<float> >& fxyz)
{
    try
    {
        size_t RO = fx.get_size(0);
        size_t E1 = fy.get_size(0);
        size_t E2 = fz.get_size(0);

        fxyz.create(RO, E1, E2);
        std::complex<float> * pFxyz = fxyz.begin();

        size_t x, y, z;

        for (z = 0; z<E2; z++)
        {
            for (y = 0; y<E1; y++)
            {
                for (x = 0; x<RO; x++)
                {
                    pFxyz[z*RO*E1 + y*RO + x] = std::complex<float>(fx(x)*fy(y)*fz(z));
                }
            }
        }
    }
    catch (...)
    {
        GADGET_THROW("Errors in compute_3d_filter(float) ... ");
    }
}

// ------------------------------------------------------------------------

void compute_3d_filter(const hoNDArray<double>& fx, const hoNDArray<double>& fy, const hoNDArray<double>& fz, hoNDArray< std::complex<double> >& fxyz)
{
    try
    {
        size_t RO = fx.get_size(0);
        size_t E1 = fy.get_size(0);
        size_t E2 = fz.get_size(0);

        fxyz.create(RO, E1, E2);
        std::complex<double> * pFxyz = fxyz.begin();

        size_t x, y, z;

        for (z = 0; z<E2; z++)
        {
            for (y = 0; y<E1; y++)
            {
                for (x = 0; x<RO; x++)
                {
                    pFxyz[z*RO*E1 + y*RO + x] = std::complex<double>(fx(x)*fy(y)*fz(z));
                }
            }
        }
    }
    catch (...)
    {
        GADGET_THROW("Errors in compute_3d_filter(double) ... ");
    }
}

// ------------------------------------------------------------------------

template <typename T>
void apply_kspace_filter_RO(hoNDArray<T>& data, const hoNDArray<T>& fRO)
{
    try
    {
        GADGET_CHECK_THROW(data.get_size(0) == fRO.get_number_of_elements());
        Gadgetron::multiply(data, fRO, data);
    }
    catch (...)
    {
        GADGET_THROW("Errors in apply_kspace_filter_RO(...) ... ");
    }
}

template EXPORTMRICORE void apply_kspace_filter_RO(hoNDArray<float>& data, const hoNDArray<float>& fRO);
template EXPORTMRICORE void apply_kspace_filter_RO(hoNDArray<double>& data, const hoNDArray<double>& fRO);
template EXPORTMRICORE void apply_kspace_filter_RO(hoNDArray< std::complex<float> >& data, const hoNDArray< std::complex<float> >& fRO);
template EXPORTMRICORE void apply_kspace_filter_RO(hoNDArray< std::complex<double> >& data, const hoNDArray< std::complex<double> >& fRO);

template <typename T> 
void apply_kspace_filter_RO(const hoNDArray<T>& data, const hoNDArray<T>& fRO, hoNDArray<T>& dataFiltered)
{
    try
    {
        GADGET_CHECK_THROW(data.get_size(0) == fRO.get_number_of_elements());
        Gadgetron::multiply(data, fRO, dataFiltered);
    }
    catch (...)
    {
        GADGET_THROW("Errors in apply_kspace_filter_RO(...) ... ");
    }
}

template EXPORTMRICORE void apply_kspace_filter_RO(const hoNDArray<float>& data, const hoNDArray<float>& fRO, hoNDArray<float>& dataFiltered);
template EXPORTMRICORE void apply_kspace_filter_RO(const hoNDArray<double>& data, const hoNDArray<double>& fRO, hoNDArray<double>& dataFiltered);
template EXPORTMRICORE void apply_kspace_filter_RO(const hoNDArray< std::complex<float> >& data, const hoNDArray< std::complex<float> >& fRO, hoNDArray< std::complex<float> >& dataFiltered);
template EXPORTMRICORE void apply_kspace_filter_RO(const hoNDArray< std::complex<double> >& data, const hoNDArray< std::complex<double> >& fRO, hoNDArray< std::complex<double> >& dataFiltered);

// ------------------------------------------------------------------------

template <typename T> 
void apply_kspace_filter_E1(const hoNDArray<T>& data, const hoNDArray<T>& fE1, hoNDArray<T>& dataFiltered)
{
    try
    {
        GADGET_CHECK_THROW(data.get_size(1) == fE1.get_number_of_elements());

        hoNDArray<T> fRO(data.get_size(0));
        fRO.fill(T(1.0));

        hoNDArray<T> fxy;
        compute_2d_filter(fRO, fE1, fxy);

        Gadgetron::multiply(data, fxy, dataFiltered);
    }
    catch (...)
    {
        GADGET_THROW("Errors in apply_kspace_filter_E1(...) ... ");
    }
}

template EXPORTMRICORE void apply_kspace_filter_E1(const hoNDArray<float>& data, const hoNDArray<float>& fE1, hoNDArray<float>& dataFiltered);
template EXPORTMRICORE void apply_kspace_filter_E1(const hoNDArray<double>& data, const hoNDArray<double>& fE1, hoNDArray<double>& dataFiltered);
template EXPORTMRICORE void apply_kspace_filter_E1(const hoNDArray< std::complex<float> >& data, const hoNDArray< std::complex<float> >& fE1, hoNDArray< std::complex<float> >& dataFiltered);
template EXPORTMRICORE void apply_kspace_filter_E1(const hoNDArray< std::complex<double> >& data, const hoNDArray< std::complex<double> >& fE1, hoNDArray< std::complex<double> >& dataFiltered);

// ------------------------------------------------------------------------

template <typename T> 
void apply_kspace_filter_ROE1(const hoNDArray<T>& data, const hoNDArray<T>& fROE1, hoNDArray<T>& dataFiltered)
{
    try
    {
        GADGET_CHECK_THROW(data.get_size(0) == fROE1.get_size(0));
        GADGET_CHECK_THROW(data.get_size(1) == fROE1.get_size(1));

        Gadgetron::multiply(data, fROE1, dataFiltered);
    }
    catch (...)
    {
        GADGET_THROW("Errors in apply_kspace_filter_ROE1(...) ... ");
    }
}

template EXPORTMRICORE void apply_kspace_filter_ROE1(const hoNDArray<float>& data, const hoNDArray<float>& fROE1, hoNDArray<float>& dataFiltered);
template EXPORTMRICORE void apply_kspace_filter_ROE1(const hoNDArray<double>& data, const hoNDArray<double>& fROE1, hoNDArray<double>& dataFiltered);
template EXPORTMRICORE void apply_kspace_filter_ROE1(const hoNDArray< std::complex<float> >& data, const hoNDArray< std::complex<float> >& fROE1, hoNDArray< std::complex<float> >& dataFiltered);
template EXPORTMRICORE void apply_kspace_filter_ROE1(const hoNDArray< std::complex<double> >& data, const hoNDArray< std::complex<double> >& fROE1, hoNDArray< std::complex<double> >& dataFiltered);

// ------------------------------------------------------------------------

template <typename T> 
void apply_kspace_filter_ROE1(const hoNDArray<T>& data, const hoNDArray<T>& fRO, const hoNDArray<T>& fE1, hoNDArray<T>& dataFiltered)
{
    try
    {
        GADGET_CHECK_THROW(data.get_size(0) == fRO.get_size(0));
        GADGET_CHECK_THROW(data.get_size(1) == fE1.get_size(0));

        hoNDArray<T> fROE1;
        compute_2d_filter(fRO, fE1, fROE1);

        Gadgetron::multiply(data, fROE1, dataFiltered);
    }
    catch (...)
    {
        GADGET_THROW("Errors in apply_kspace_filter_ROE1(...) ... ");
    }
}

template EXPORTMRICORE void apply_kspace_filter_ROE1(const hoNDArray<float>& data, const hoNDArray<float>& fRO, const hoNDArray<float>& fE1, hoNDArray<float>& dataFiltered);
template EXPORTMRICORE void apply_kspace_filter_ROE1(const hoNDArray<double>& data, const hoNDArray<double>& fRO, const hoNDArray<double>& fE1, hoNDArray<double>& dataFiltered);
template EXPORTMRICORE void apply_kspace_filter_ROE1(const hoNDArray< std::complex<float> >& data, const hoNDArray< std::complex<float> >& fRO, const hoNDArray< std::complex<float> >& fE1, hoNDArray< std::complex<float> >& dataFiltered);
template EXPORTMRICORE void apply_kspace_filter_ROE1(const hoNDArray< std::complex<double> >& data, const hoNDArray< std::complex<double> >& fRO, const hoNDArray< std::complex<double> >& fE1, hoNDArray< std::complex<double> >& dataFiltered);

// ------------------------------------------------------------------------

template <typename T> 
void apply_kspace_filter_E2(const hoNDArray<T>& data, const hoNDArray<T>& fE2, hoNDArray<T>& dataFiltered)
{
    try
    {
        GADGET_CHECK_THROW(data.get_size(2) == fE2.get_number_of_elements());

        hoNDArray<T> fRO(data.get_size(0));
        fRO.fill(T(1.0));

        hoNDArray<T> fE1(data.get_size(1));
        fE1.fill(T(1.0));

        hoNDArray<T> fxyz;
        compute_3d_filter(fRO, fE1, fE2, fxyz);

        Gadgetron::multiply(data, fxyz, dataFiltered);
    }
    catch (...)
    {
        GADGET_THROW("Errors in apply_kspace_filter_E2(...) ... ");
    }
}

template EXPORTMRICORE void apply_kspace_filter_E2(const hoNDArray<float>& data, const hoNDArray<float>& fE2, hoNDArray<float>& dataFiltered);
template EXPORTMRICORE void apply_kspace_filter_E2(const hoNDArray<double>& data, const hoNDArray<double>& fE2, hoNDArray<double>& dataFiltered);
template EXPORTMRICORE void apply_kspace_filter_E2(const hoNDArray< std::complex<float> >& data, const hoNDArray< std::complex<float> >& fE2, hoNDArray< std::complex<float> >& dataFiltered);
template EXPORTMRICORE void apply_kspace_filter_E2(const hoNDArray< std::complex<double> >& data, const hoNDArray< std::complex<double> >& fE2, hoNDArray< std::complex<double> >& dataFiltered);

// ------------------------------------------------------------------------

template <typename T> 
void apply_kspace_filter_ROE2(const hoNDArray<T>& data, const hoNDArray<T>& fRO, const hoNDArray<T>& fE2, hoNDArray<T>& dataFiltered)
{
    try
    {
        GADGET_CHECK_THROW(data.get_size(0) == fRO.get_number_of_elements());
        GADGET_CHECK_THROW(data.get_size(2) == fE2.get_number_of_elements());

        hoNDArray<T> fE1(data.get_size(1));
        fE1.fill(T(1.0));

        hoNDArray<T> fxyz;
        compute_3d_filter(fRO, fE1, fE2, fxyz);

        Gadgetron::multiply(data, fxyz, dataFiltered);
    }
    catch (...)
    {
        GADGET_THROW("Errors in apply_kspace_filter_ROE2(...) ... ");
    }
}

template EXPORTMRICORE void apply_kspace_filter_ROE2(const hoNDArray<float>& data, const hoNDArray<float>& fRO, const hoNDArray<float>& fE2, hoNDArray<float>& dataFiltered);
template EXPORTMRICORE void apply_kspace_filter_ROE2(const hoNDArray<double>& data, const hoNDArray<double>& fRO, const hoNDArray<double>& fE2, hoNDArray<double>& dataFiltered);
template EXPORTMRICORE void apply_kspace_filter_ROE2(const hoNDArray< std::complex<float> >& data, const hoNDArray< std::complex<float> >& fRO, const hoNDArray< std::complex<float> >& fE2, hoNDArray< std::complex<float> >& dataFiltered);
template EXPORTMRICORE void apply_kspace_filter_ROE2(const hoNDArray< std::complex<double> >& data, const hoNDArray< std::complex<double> >& fRO, const hoNDArray< std::complex<double> >& fE2, hoNDArray< std::complex<double> >& dataFiltered);

// ------------------------------------------------------------------------

template <typename T> 
void apply_kspace_filter_E1E2(const hoNDArray<T>& data, const hoNDArray<T>& fE1, const hoNDArray<T>& fE2, hoNDArray<T>& dataFiltered)
{
    try
    {
        GADGET_CHECK_THROW(data.get_size(1) == fE1.get_number_of_elements());
        GADGET_CHECK_THROW(data.get_size(2) == fE2.get_number_of_elements());

        hoNDArray<T> fRO(data.get_size(0));
        fRO.fill(T(1.0));

        hoNDArray<T> fxyz;
        compute_3d_filter(fRO, fE1, fE2, fxyz);

        Gadgetron::multiply(data, fxyz, dataFiltered);
    }
    catch (...)
    {
        GADGET_THROW("Errors in apply_kspace_filter_E1E2(...) ... ");
    }
}

template EXPORTMRICORE void apply_kspace_filter_E1E2(const hoNDArray<float>& data, const hoNDArray<float>& fE1, const hoNDArray<float>& fE2, hoNDArray<float>& dataFiltered);
template EXPORTMRICORE void apply_kspace_filter_E1E2(const hoNDArray<double>& data, const hoNDArray<double>& fE1, const hoNDArray<double>& fE2, hoNDArray<double>& dataFiltered);
template EXPORTMRICORE void apply_kspace_filter_E1E2(const hoNDArray< std::complex<float> >& data, const hoNDArray< std::complex<float> >& fE1, const hoNDArray< std::complex<float> >& fE2, hoNDArray< std::complex<float> >& dataFiltered);
template EXPORTMRICORE void apply_kspace_filter_E1E2(const hoNDArray< std::complex<double> >& data, const hoNDArray< std::complex<double> >& fE1, const hoNDArray< std::complex<double> >& fE2, hoNDArray< std::complex<double> >& dataFiltered);

// ------------------------------------------------------------------------

template <typename T> 
void apply_kspace_filter_ROE1E2(const hoNDArray<T>& data, const hoNDArray<T>& fROE1E2, hoNDArray<T>& dataFiltered)
{
    try
    {
        GADGET_CHECK_THROW(data.get_size(0) == fROE1E2.get_size(0));
        GADGET_CHECK_THROW(data.get_size(1) == fROE1E2.get_size(1));
        GADGET_CHECK_THROW(data.get_size(2) == fROE1E2.get_size(2));

        Gadgetron::multiply(data, fROE1E2, dataFiltered);
    }
    catch (...)
    {
        GADGET_THROW("Errors in apply_kspace_filter_ROE1E2(...) ... ");
    }
}

template EXPORTMRICORE void apply_kspace_filter_ROE1E2(const hoNDArray<float>& data, const hoNDArray<float>& fROE1E2, hoNDArray<float>& dataFiltered);
template EXPORTMRICORE void apply_kspace_filter_ROE1E2(const hoNDArray<double>& data, const hoNDArray<double>& fROE1E2, hoNDArray<double>& dataFiltered);
template EXPORTMRICORE void apply_kspace_filter_ROE1E2(const hoNDArray< std::complex<float> >& data, const hoNDArray< std::complex<float> >& fROE1E2, hoNDArray< std::complex<float> >& dataFiltered);
template EXPORTMRICORE void apply_kspace_filter_ROE1E2(const hoNDArray< std::complex<double> >& data, const hoNDArray< std::complex<double> >& fROE1E2, hoNDArray< std::complex<double> >& dataFiltered);

// ------------------------------------------------------------------------

template <typename T> 
void apply_kspace_filter_ROE1E2(const hoNDArray<T>& data, const hoNDArray<T>& fRO, const hoNDArray<T>& fE1, const hoNDArray<T>& fE2, hoNDArray<T>& dataFiltered)
{
    try
    {
        GADGET_CHECK_THROW(data.get_size(0) == fRO.get_number_of_elements());
        GADGET_CHECK_THROW(data.get_size(1) == fE1.get_number_of_elements());
        GADGET_CHECK_THROW(data.get_size(2) == fE2.get_number_of_elements());

        hoNDArray<T> fxyz;
        compute_3d_filter(fRO, fE1, fE2, fxyz);

        Gadgetron::multiply(data, fxyz, dataFiltered);
    }
    catch (...)
    {
        GADGET_THROW("Errors in apply_kspace_filter_ROE1E2(...) ... ");
    }
}

template EXPORTMRICORE void apply_kspace_filter_ROE1E2(const hoNDArray<float>& data, const hoNDArray<float>& fRO, const hoNDArray<float>& fE1, const hoNDArray<float>& fE2, hoNDArray<float>& dataFiltered);
template EXPORTMRICORE void apply_kspace_filter_ROE1E2(const hoNDArray<double>& data, const hoNDArray<double>& fRO, const hoNDArray<double>& fE1, const hoNDArray<double>& fE2, hoNDArray<double>& dataFiltered);
template EXPORTMRICORE void apply_kspace_filter_ROE1E2(const hoNDArray< std::complex<float> >& data, const hoNDArray< std::complex<float> >& fRO, const hoNDArray< std::complex<float> >& fE1, const hoNDArray< std::complex<float> >& fE2, hoNDArray< std::complex<float> >& dataFiltered);
template EXPORTMRICORE void apply_kspace_filter_ROE1E2(const hoNDArray< std::complex<double> >& data, const hoNDArray< std::complex<double> >& fRO, const hoNDArray< std::complex<double> >& fE1, const hoNDArray< std::complex<double> >& fE2, hoNDArray< std::complex<double> >& dataFiltered);

// ------------------------------------------------------------------------

void find_symmetric_sampled_region(size_t start, size_t end, size_t center, size_t& startSym, size_t& endSym)
{
    GADGET_CHECK_THROW(end >= start);
    GADGET_CHECK_THROW(center >= start);
    GADGET_CHECK_THROW(end >= center);

    size_t halfSizeStart = center - start;
    size_t halfSizeEnd = end - center;

    if (halfSizeStart > halfSizeEnd)
    {
        startSym = center - halfSizeEnd;
        endSym = center + halfSizeEnd;
    }
    else
    {
        startSym = center - halfSizeStart;
        endSym = center + halfSizeStart;
    }
}

// ------------------------------------------------------------------------

template<typename T>
void compute_filter_SNR_unit_scale_factor(const hoNDArray<T>& filter, T& scalFactor)
{
    size_t ii, len;

    len = filter.get_number_of_elements();
    if (len == 0)
    {
        scalFactor = T(1.0);
        return;
    }

    T sos(0.0);
    for (ii = 0; ii<len; ii++)
    {
        sos += filter(ii)*filter(ii);
    }

    scalFactor = (T)(1.0 / std::sqrt(std::abs(sos) / len));
}

template EXPORTMRICORE void compute_filter_SNR_unit_scale_factor(const hoNDArray<float>& filter, float& scalFactor);
template EXPORTMRICORE void compute_filter_SNR_unit_scale_factor(const hoNDArray<double>& filter, double& scalFactor);

template<> EXPORTMRICORE
void compute_filter_SNR_unit_scale_factor(const hoNDArray< std::complex<float> >& filter, std::complex<float> & scalFactor)
{
    size_t ii, len;

    len = filter.get_number_of_elements();
    if (len == 0)
    {
        scalFactor = std::complex<float>(1.0);
        return;
    }

    std::complex<float> sos(0.0);
    for (ii = 0; ii<len; ii++)
    {
        sos += filter(ii)*std::conj(filter(ii));
    }

    scalFactor = (std::complex<float>)(1.0 / std::sqrt(std::abs(sos) / len));
}

template<> EXPORTMRICORE
void compute_filter_SNR_unit_scale_factor(const hoNDArray< std::complex<double> >& filter, std::complex<double> & scalFactor)
{
    size_t ii, len;

    len = filter.get_number_of_elements();
    if (len == 0)
    {
        scalFactor = std::complex<double>(1.0);
        return;
    }

    std::complex<double> sos(0.0);
    for (ii = 0; ii<len; ii++)
    {
        sos += filter(ii)*std::conj(filter(ii));
    }

    scalFactor = (std::complex<double>)(1.0 / std::sqrt(std::abs(sos) / len));
}

// ------------------------------------------------------------------------

}
