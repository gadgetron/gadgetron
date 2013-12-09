/********************************************************************
    created:    2013/10/03
    created:    3:10:2013   14:06
    author:     Hui Xue

    purpose:    Gadgetron complex data structure to matlab conversion
*********************************************************************/

#pragma once

#include <matrix.h>
#include <mat.h>
#include <mexGT.h>
#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <strstream>

#include "hoNDArray.h"

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

    virtual void printInfo(std::ostream& os) const;

protected:
};

template <typename T> 
bool gtMatlabConverterComplex<T>::
hoNDArray2Matlab(const hoNDArray<T>& a, mxArray*& aMx)
{
    try
    {
        boost::shared_ptr< std::vector<unsigned long long> > dim = a.get_dimensions();

        int ndim = dim->size();
        mwSize* dims = new mwSize[ndim];

        unsigned long long ii;
        for ( ii=0; ii<ndim; ii++ )
        {
            dims[ii] = static_cast<mwSize>( *dim[ii] );
        }

        unsigned long long N = a.get_number_of_elements();
        const T* pA = a.begin();

        if ( typeid(T) == typeid(std::complex<float) )
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
        else if ( typeid(T) == typeid(std::complex<double) )
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

        std::vector<unsigned long long> dim(ndim);

        unsigned long long ii;
        for ( ii=0; ii<ndim; ii++ )
        {
            dim[ii] = static_cast<unsigned long long>(dims[ii]);
        }

        a.create(&dim);
        unsigned long long N = a.get_number_of_elements();
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
