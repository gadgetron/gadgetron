/********************************************************************
    created:    2014/02/25
    author:     Hui Xue

    purpose:    Gadgetron data structure for ND image matlab conversion
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
#include "hoNDImage.h"
#include "gtMatlab.h"

namespace Gadgetron
{

template <typename T, unsigned int D>
class gtMatlabImageHeader
{
public:

    typedef hoNDImage<T, D> ImageType;

    typedef typename ImageType::value_type value_type;
    typedef typename ImageType::coord_type coord_type;
    typedef typename ImageType::a_axis_type a_axis_type;
    typedef typename ImageType::axis_type axis_type;

    coord_type pixelSize_[D];
    coord_type origin_[D];
    hoNDPoint<coord_type, D> axis_[D];

    gtMatlabImageHeader();
    gtMatlabImageHeader(const ImageType& im);
    virtual ~gtMatlabImageHeader();

    /// for the axis, it will be a D*D rotation matrix
    /// every column is a oritentation vector for a dimension
    virtual bool toMatlab(mxArray*& header);
    virtual bool fromMatlab(const mxArray* header);

protected:

    // the header field names
    std::vector<char*> header_fields_;

    void set_header_fields()
    {
        size_t num = 3; // origin, pixelSize, axis
        header_fields_.resize(3);
        header_fields_[0] = "origin";
        header_fields_[1] = "pixelSize";
        header_fields_[2] = "axis";
    }
};

template <typename T, unsigned int D>
gtMatlabImageHeader<T, D>::gtMatlabImageHeader()
{
    unsigned int ii;
    for (ii=0;ii<D; ii++)
    {
        pixelSize_[ii] = 1;
        origin_[ii] = 0;
        axis_[ii].fill(0);
        axis_[ii][ii] = coord_type(1.0);
    }

    this->set_header_fields();
}

template <typename T, unsigned int D>
gtMatlabImageHeader<T, D>::gtMatlabImageHeader(const ImageType& im)
{
    std::vector<coord_type> pixelSize;
    im.get_pixel_size(pixelSize);

    std::vector<coord_type> origin;
    im.get_origin(origin);

    axis_type axis;
    im.get_axis(axis);

    unsigned int ii;
    for (ii=0;ii<D; ii++)
    {
        pixelSize_[ii] = pixelSize[ii];
        origin_[ii] = origin[ii];
        axis_[ii] = axis[ii];
    }

    this->set_header_fields();
}

template <typename T, unsigned int D>
gtMatlabImageHeader<T, D>::~gtMatlabImageHeader()
{

}

template <typename T, unsigned int D>
bool gtMatlabImageHeader<T, D>::toMatlab(mxArray*& header)
{
    try
    {
        unsigned int ii, jj;

        mwSize num[2] = {1, 1};
        header = mxCreateStructArray(2, num, (int)header_fields_.size(), const_cast<const char**>(&header_fields_[0]));

        mwSize dims[1];
        dims[0] = D;

        mxArray* aMx = mxCreateNumericArray(1, dims, mxSINGLE_CLASS, mxREAL);
        float* pr = static_cast<float*>(mxGetData(aMx));
        for ( ii=0; ii<D; ii++ )
        {
            pr[ii] = origin_[ii];
        }

        mxSetField(header, 0, header_fields_[0], aMx);

        aMx = mxCreateNumericArray(1, dims, mxSINGLE_CLASS, mxREAL);
        pr = static_cast<float*>(mxGetData(aMx));
        for ( ii=0; ii<D; ii++ )
        {
            pr[ii] = pixelSize_[ii];
        }

        mxSetField(header, 0, header_fields_[1], aMx);

        mwSize dimsAxis[2];
        dimsAxis[0] = D;
        dimsAxis[1] = D;

        aMx = mxCreateNumericMatrix(D, D, mxSINGLE_CLASS, mxREAL);
        pr = static_cast<float*>(mxGetData(aMx));
        for ( jj=0; jj<D; jj++ )
        {
            for ( ii=0; ii<D; ii++ )
            {
                pr[jj + ii*D] = axis_[jj][ii]; // stored in column-wise
            }
        }

        mxSetField(header, 0, header_fields_[2], aMx);
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Error happened in gtMatlabImageHeader<T, D>::toMatlab(mxArray*& header) ... ");
        return false;
    }

    return true;
}

template <typename T, unsigned int D>
bool gtMatlabImageHeader<T, D>::fromMatlab(const mxArray* header)
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(mxIsStruct(header));

        unsigned int ii, jj;

        mxArray* aMx = mxGetField(header, 0, header_fields_[0]);
        size_t N = mxGetNumberOfElements(aMx);

        if ( mxIsSingle(aMx) )
        {
            float* pr = static_cast<float*>(mxGetData(aMx));

            for ( ii=0; ii<GT_MIN(D, N); ii++ )
            {
                origin_[ii] = (coord_type)pr[ii];
            }
        }
        else
        {
            double* pr = static_cast<double*>(mxGetData(aMx));

            for ( ii=0; ii<GT_MIN(D, N); ii++ )
            {
                origin_[ii] = (coord_type)pr[ii];
            }
        }

        aMx = mxGetField(header, 0, header_fields_[1]);
        N = mxGetNumberOfElements(aMx);

        if ( mxIsSingle(aMx) )
        {
            float* pr = static_cast<float*>(mxGetData(aMx));

            for ( ii=0; ii<GT_MIN(D, N); ii++ )
            {
                pixelSize_[ii] = (coord_type)pr[ii];
            }
        }
        else
        {
            double* pr = static_cast<double*>(mxGetData(aMx));

            for ( ii=0; ii<GT_MIN(D, N); ii++ )
            {
                pixelSize_[ii] = (coord_type)pr[ii];
            }
        }

        aMx = mxGetField(header, 0, header_fields_[2]);

        if ( mxIsSingle(aMx) )
        {
            float* pr = static_cast<float*>(mxGetData(aMx));

            for ( jj=0; jj<GT_MIN(D, N); jj++ )
            {
                for ( ii=0; ii<GT_MIN(D, N); ii++ )
                {
                    axis_[jj][ii] = (coord_type)pr[jj + ii*D];
                }
            }
        }
        else
        {
            double* pr = static_cast<double*>(mxGetData(aMx));

            for ( jj=0; jj<GT_MIN(D, N); jj++ )
            {
                for ( ii=0; ii<GT_MIN(D, N); ii++ )
                {
                    axis_[jj][ii] = (coord_type)pr[jj + ii*D];
                }
            }
        }
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Error happened in gtMatlabImageHeader<T, D>::fromMatlab(const mxArray* header) ... ");
        return false;
    }

    return true;
}

}
