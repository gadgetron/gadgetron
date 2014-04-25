/********************************************************************
    created:    2013/10/03
    created:    3:10:2013   14:06
    author:     Hui Xue

    purpose:    Gadgetron data structure to matlab conversion
*********************************************************************/

#pragma once

#include "gtMatlabImage.h"

namespace Gadgetron
{

template <typename T> 
class gtMatlabConverter
{
public:

    typedef gtMatlabConverter<T> Self;

    gtMatlabConverter() {}
    virtual ~gtMatlabConverter() {}

    virtual bool hoNDArray2Matlab(const hoNDArray<T>& a, mxArray*& aMx);
    virtual bool Matlab2hoNDArray(const mxArray* aMx, hoNDArray<T>& a);

    virtual bool Vec2Matlab(const std::vector<T>& vec, mxArray*& aMx);
    virtual bool Matlab2Vec(const mxArray* aMx, std::vector<T>& vec);

    virtual bool Str2Matlab(const std::string& str, mxArray*& aMx);
    virtual bool Matlab2Str(const mxArray* aMx, std::string& str);

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
        mwSize ndim = mxGetNumberOfDimensions(aMx);
        if ( ndim != D ) return false;

        hoNDArray<T> buf;
        GADGET_CHECK_RETURN_FALSE(Matlab2hoNDArray(aMx, buf));
        GADGET_CHECK_RETURN_FALSE(buf.get_number_of_dimensions()<=D);

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
};

template <typename T> 
bool gtMatlabConverter<T>::
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

        if ( typeid(T) == typeid(float) )
        {
            aMx = mxCreateNumericArray(ndim, dims, mxSINGLE_CLASS, mxREAL);
            float* ptr = static_cast<float*>(mxGetData(aMx));

            for ( ii=0; ii<N; ii++ )
            {
                ptr[ii] = pA[ii];
            }
        }
        else
        {
            aMx = mxCreateNumericArray(ndim, dims, mxDOUBLE_CLASS, mxREAL);
            double* ptr = static_cast<double*>(mxGetData(aMx));

            for ( ii=0; ii<N; ii++ )
            {
                ptr[ii] = pA[ii];
            }
        }
    }
    catch(...)
    {
        mexErrMsgTxt("Errors happened in gtMatlabConverter::hoNDArray2Matlab(const hoNDArray<T>& a, mxArray*& aMx) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtMatlabConverter<T>::
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

        if ( mxIsSingle(aMx) )
        {
            float* ptr = static_cast<float*>(mxGetData(aMx));
            for ( ii=0; ii<N; ii++ )
            {
                pA[ii] = static_cast<T>(ptr[ii]);
            }
        }
        else
        {
            double* ptr = static_cast<double*>(mxGetData(aMx));
            for ( ii=0; ii<N; ii++ )
            {
                pA[ii] = static_cast<T>(ptr[ii]);
            }
        }
    }
    catch(...)
    {
        mexErrMsgTxt("Errors happened in gtMatlabConverter::Matlab2hoNDArray(const mxArray* aMx, hoNDArray<T>& a) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtMatlabConverter<T>::
Vec2Matlab(const std::vector<T>& vec, mxArray*& aMx)
{
    try
    {
        aMx = mxCreateNumericMatrix(vec.size(), 1, mxDOUBLE_CLASS, mxREAL);
        double* ptr = static_cast<double*>(mxGetData(aMx));
        for ( size_t ii=0; ii<vec.size(); ii++ )
        {
            ptr[ii] = static_cast<double>(vec[ii]);
        }
    }
    catch(...)
    {
        mexErrMsgTxt("Errors happened in gtMatlabConverter::Vec2Matlab(const std::vector<T>& vec, mxArray*& aMx) ... ");
        return false;
    }

    return true; 
}

template <typename T> 
bool gtMatlabConverter<T>::
Matlab2Vec(const mxArray* aMx, std::vector<T>& vec)
{
    try
    {
        mwSize N = mxGetNumberOfElements(aMx);
        vec.resize(N);

        if ( mxIsSingle(aMx) )
        {
            float* ptr = static_cast<float*>(mxGetData(aMx));
            for ( size_t ii=0; ii<N; ii++ )
            {
                vec[ii] = static_cast<T>(ptr[ii]);
            }
        }
        else
        {
            double* ptr = static_cast<double*>(mxGetData(aMx));
            for ( size_t ii=0; ii<N; ii++ )
            {
                vec[ii] = static_cast<T>(ptr[ii]);
            }
        }
    }
    catch(...)
    {
        mexErrMsgTxt("Errors happened in gtMatlabConverter::Matlab2Vec(const mxArray* aMx, std::vector<T>& vec) ... ");
        return false;
    }

    return true; 
}

template <typename T> 
bool gtMatlabConverter<T>::
Str2Matlab(const std::string& str, mxArray*& aMx)
{
    aMx = mxCreateString(str.c_str());
    return (aMx != NULL);
}

template <typename T> 
bool gtMatlabConverter<T>::
Matlab2Str(const mxArray* aMx, std::string& str)
{
    int N = mxGetNumberOfElements(aMx) + 1;

    std::vector<char> buf(N, '\0');
    if (mxGetString(aMx, &buf[0], N) != 0)
    {
        return false;
    }
    str = std::string(&buf[0]);

    return true;
}

template <typename T> 
void gtMatlabConverter<T>::printInfo(std::ostream& os) const
{
    using namespace std;
    os << "--------------------------------------------------" << endl;
    os << "Gadgetron matlab Converter ..." << endl;
    os << "--------------------------------------------------" << endl;
}

}
