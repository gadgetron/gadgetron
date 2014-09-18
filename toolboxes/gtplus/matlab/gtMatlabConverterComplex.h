/********************************************************************
    created:    2013/10/03
    created:    3:10:2013   14:06
    author:     Hui Xue

    purpose:    Gadgetron complex data structure to matlab conversion
*********************************************************************/

#pragma once

#include "gtMatlabImage.h"

namespace Gadgetron
{

template <typename T> 
class gtMatlabConverterComplex
{
public:

    gtMatlabConverterComplex() {}
    virtual ~gtMatlabConverterComplex() {}

    virtual bool hoNDArray2Matlab(const hoNDArray<T>& a, mxArray*& aMx);
    virtual bool Matlab2hoNDArray(const mxArray* aMx, hoNDArray<T>& a);

    template <unsigned int D> 
    bool hoNDImage2Matlab(const hoNDImage<T, D>& a, mxArray*& aMx, mxArray*& aHeader)
    {
        std::vector<size_t> dim(D);
        a.get_dimensions(dim);

        hoNDArray<T> buf(dim, const_cast<T*>(a.get_data_ptr()), false);
        GADGET_CHECK_RETURN_FALSE(hoNDArray2Matlab(buf, aMx));

        gtMatlabImageHeader<T, D> header(a);
        GADGET_CHECK_RETURN_FALSE(header.toMatlab(aHeader));

        return true;
    }

    template <unsigned int D> 
    bool Matlab2hoNDImage(const mxArray* aMx, const mxArray* aHeader, hoNDImage<T, D>& a)
    {
        hoNDArray<T> buf;
        GADGET_CHECK_RETURN_FALSE(Matlab2hoNDArray(aMx, buf));
        GADGET_CHECK_RETURN_FALSE(buf.get_number_of_dimensions()==D);

        a.from_NDArray(buf);

        gtMatlabImageHeader<T, D> header;
        GADGET_CHECK_RETURN_FALSE(header.fromMatlab(aHeader));

        unsigned int ii;
        for ( ii=0; ii<D; ii++ )
        {
            a.set_pixel_size(ii, header.pixelSize_[ii]);
            a.set_origin(ii, header.origin_[ii]);
            a.set_axis(ii, header.axis_[ii]);
        }

        return true;
    }

    virtual void printInfo(std::ostream& os) const;

protected:
};

template <typename T> 
bool gtMatlabConverterComplex<T>::
hoNDArray2Matlab(const hoNDArray<T>& a, mxArray*& aMx)
{
    try
    {
        boost::shared_ptr< std::vector<size_t> > dim = a.get_dimensions();

        int ndim = dim->size();
        mwSize* dims = new mwSize[ndim];

        size_t ii;
        for ( ii=0; ii<ndim; ii++ )
        {
            dims[ii] = static_cast<mwSize>( (*dim)[ii] );
        }

        size_t N = a.get_number_of_elements();
        const T* pA = a.begin();

        if ( typeid(T) == typeid(std::complex<float>) )
        {
            aMx = mxCreateNumericArray(ndim, dims, mxSINGLE_CLASS, mxCOMPLEX);
            float* pr = static_cast<float*>(mxGetData(aMx));
            float* pi = static_cast<float*>(mxGetImagData(aMx));

            for ( ii=0; ii<N; ii++ )
            {
                pr[ii] = static_cast<float>(pA[ii].real());
                pi[ii] = static_cast<float>(pA[ii].imag());
            }
        }
        else if ( typeid(T) == typeid(std::complex<double>) )
        {
            aMx = mxCreateNumericArray(ndim, dims, mxDOUBLE_CLASS, mxCOMPLEX);
            double* pr = static_cast<double*>(mxGetData(aMx));
            double* pi = static_cast<double*>(mxGetImagData(aMx));

            for ( ii=0; ii<N; ii++ )
            {
                pr[ii] = static_cast<double>(pA[ii].real());
                pi[ii] = static_cast<double>(pA[ii].imag());
            }
        }
    }
    catch(...)
    {
        mexErrMsgTxt("Errors happened in gtMatlabConverterComplex::hoNDArray2Matlab(const hoNDArray<T>& a, mxArray*& aMx) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtMatlabConverterComplex<T>::
Matlab2hoNDArray(const mxArray* aMx, hoNDArray<T>& a)
{
    try
    {
        mwSize ndim = mxGetNumberOfDimensions(aMx);
        const mwSize* dims = mxGetDimensions(aMx);

        std::vector<size_t> dim(ndim);

        size_t ii;
        for ( ii=0; ii<ndim; ii++ )
        {
            dim[ii] = static_cast<size_t>(dims[ii]);
        }

        a.create(&dim);
        size_t N = a.get_number_of_elements();
        T* pA = a.begin();

        if ( mxIsComplex(aMx) && mxIsDouble(aMx) )
        {
            double* pr = static_cast<double*>(mxGetData(aMx));
            double* pi = static_cast<double*>(mxGetImagData(aMx));

            for ( ii=0; ii<N; ii++ )
            {
                pA[ii] = T(pr[ii], pi[ii]);
            }
        }
        else if ( mxIsComplex(aMx) && mxIsSingle(aMx) )
        {
            float* pr = static_cast<float*>(mxGetData(aMx));
            float* pi = static_cast<float*>(mxGetImagData(aMx));

            for ( ii=0; ii<N; ii++ )
            {
                pA[ii] = T(pr[ii], pi[ii]);
            }
        }
    }
    catch(...)
    {
        mexErrMsgTxt("Errors happened in gtMatlabConverterComplex::Matlab2hoNDArray(const mxArray* aMx, hoNDArray<T>& a) ... ");
        return false;
    }

    return true;
}

template <typename T> 
void gtMatlabConverterComplex<T>::printInfo(std::ostream& os) const
{
    using namespace std;
    os << "--------------------------------------------------" << endl;
    os << "Gadgetron matlab Converter for complex arrays ..." << endl;
    os << "--------------------------------------------------" << endl;
}

}
