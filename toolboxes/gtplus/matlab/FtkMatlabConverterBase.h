/**
*  @file    FtkMatlabConverterBase.h
*  @brief   Ftk and Matlab converter base
*  @author  Hui Xue
*  @date    July 18, 2011
*  @Site    SCR, Princeton
*
*  Copyright (C) Siemens Corporate Research, Inc. 2011 All Rights Reserved
**/

#ifdef FTK_MATLAB_SUPPORT

#ifndef FTK_FTKMATLABCONVERTERBASE_H
#define FTK_FTKMATLABCONVERTERBASE_H

#include <interface/matlab/FtkMexExport.h>

#include <vector>
#include <typeinfo>
#include <core/basic/Clock.h>
#include <core/basic/Common.h>
#include <core/basic/Exception.h>
#include <core/math/MathMacros.h>
#include <core/basic/Allocate.h>
#include <core/basic/Array1d.h> 
#include <core/basic/RealMatrix.h> 
#include <core/image/Image.h>
#include "core/image/container/ImageContainerArray.h"
#include "core/image/container/ImageContainerMatrix.h"
#include "core/image/container/ImageContainerAllocation.h"

BEGIN_NAMESPACE_1(ftk)

struct ImageInfo
{
    int sizeX;
    int sizeY;
    int sizeZ;
    int sizeT;
    int sizeN;
    int sizeM;

    double spacingX;
    double spacingY;
    double spacingZ;
    double spacingT;
    double spacingN;
    double spacingM;

    double positionPatient[3];
    double orientationPatient[3][3];

    /// number of fields
    FTK_STATIC_CONST( numOfFields, IndexType, 20 );

    ImageInfo() 
    {
        initialize();
    }

    ImageInfo(const ImageBase<2>& aImage) 
    {
        initialize();

        sizeX = aImage.getSize(0);
        sizeY = aImage.getSize(1);

        spacingX = aImage.getPixelSpacing(0);
        spacingY = aImage.getPixelSpacing(1);

        positionPatient[0] = aImage.getPosition(0);
        positionPatient[1] = aImage.getPosition(1);
        positionPatient[2] = aImage.getPosition(2);

        orientationPatient[0][0] = aImage.getOrient3D(0, 0);
        orientationPatient[0][1] = aImage.getOrient3D(0, 1);
        orientationPatient[0][2] = aImage.getOrient3D(0, 2);

        orientationPatient[1][0] = aImage.getOrient3D(1, 0);
        orientationPatient[1][1] = aImage.getOrient3D(1, 1);
        orientationPatient[1][2] = aImage.getOrient3D(1, 2);

        orientationPatient[2][0] = aImage.getOrient3D(2, 0);
        orientationPatient[2][1] = aImage.getOrient3D(2, 1);
        orientationPatient[2][2] = aImage.getOrient3D(2, 2);
    }

    ImageInfo(const ImageBase<3>& aImage) 
    {
        initialize();

        sizeX = aImage.getSize(0);
        sizeY = aImage.getSize(1);
        sizeZ = aImage.getSize(2);

        spacingX = aImage.getPixelSpacing(0);
        spacingY = aImage.getPixelSpacing(1);
        spacingZ = aImage.getPixelSpacing(2);

        positionPatient[0] = aImage.getPosition(0);
        positionPatient[1] = aImage.getPosition(1);
        positionPatient[2] = aImage.getPosition(2);

        orientationPatient[0][0] = aImage.getOrient3D(0, 0);
        orientationPatient[0][1] = aImage.getOrient3D(0, 1);
        orientationPatient[0][2] = aImage.getOrient3D(0, 2);

        orientationPatient[1][0] = aImage.getOrient3D(1, 0);
        orientationPatient[1][1] = aImage.getOrient3D(1, 1);
        orientationPatient[1][2] = aImage.getOrient3D(1, 2);

        orientationPatient[2][0] = aImage.getOrient3D(2, 0);
        orientationPatient[2][1] = aImage.getOrient3D(2, 1);
        orientationPatient[2][2] = aImage.getOrient3D(2, 2);
    }

    ImageInfo(const ImageBase<4>& aImage) 
    {
        initialize();

        sizeX = aImage.getSize(0);
        sizeY = aImage.getSize(1);
        sizeZ = aImage.getSize(2);
        sizeT = aImage.getSize(3);

        spacingX = aImage.getPixelSpacing(0);
        spacingY = aImage.getPixelSpacing(1);
        spacingZ = aImage.getPixelSpacing(2);
        spacingT = aImage.getPixelSpacing(3);

        positionPatient[0] = aImage.getPosition(0);
        positionPatient[1] = aImage.getPosition(1);
        positionPatient[2] = aImage.getPosition(2);

        orientationPatient[0][0] = aImage.getOrient3D(0, 0);
        orientationPatient[0][1] = aImage.getOrient3D(0, 1);
        orientationPatient[0][2] = aImage.getOrient3D(0, 2);

        orientationPatient[1][0] = aImage.getOrient3D(1, 0);
        orientationPatient[1][1] = aImage.getOrient3D(1, 1);
        orientationPatient[1][2] = aImage.getOrient3D(1, 2);

        orientationPatient[2][0] = aImage.getOrient3D(2, 0);
        orientationPatient[2][1] = aImage.getOrient3D(2, 1);
        orientationPatient[2][2] = aImage.getOrient3D(2, 2);
    }

    ImageInfo(const ImageBase<5>& aImage) 
    {
        initialize();

        sizeX = aImage.getSize(0);
        sizeY = aImage.getSize(1);
        sizeZ = aImage.getSize(2);
        sizeT = aImage.getSize(3);
        sizeN = aImage.getSize(4);

        spacingX = aImage.getPixelSpacing(0);
        spacingY = aImage.getPixelSpacing(1);
        spacingZ = aImage.getPixelSpacing(2);
        spacingT = aImage.getPixelSpacing(3);
        spacingN = aImage.getPixelSpacing(4);

        positionPatient[0] = aImage.getPosition(0);
        positionPatient[1] = aImage.getPosition(1);
        positionPatient[2] = aImage.getPosition(2);

        orientationPatient[0][0] = aImage.getOrient3D(0, 0);
        orientationPatient[0][1] = aImage.getOrient3D(0, 1);
        orientationPatient[0][2] = aImage.getOrient3D(0, 2);

        orientationPatient[1][0] = aImage.getOrient3D(1, 0);
        orientationPatient[1][1] = aImage.getOrient3D(1, 1);
        orientationPatient[1][2] = aImage.getOrient3D(1, 2);

        orientationPatient[2][0] = aImage.getOrient3D(2, 0);
        orientationPatient[2][1] = aImage.getOrient3D(2, 1);
        orientationPatient[2][2] = aImage.getOrient3D(2, 2);
    }

    ImageInfo(const ImageBase<6>& aImage) 
    {
        initialize();

        sizeX = aImage.getSize(0);
        sizeY = aImage.getSize(1);
        sizeZ = aImage.getSize(2);
        sizeT = aImage.getSize(3);
        sizeN = aImage.getSize(4);
        sizeM = aImage.getSize(5);

        spacingX = aImage.getPixelSpacing(0);
        spacingY = aImage.getPixelSpacing(1);
        spacingZ = aImage.getPixelSpacing(2);
        spacingT = aImage.getPixelSpacing(3);
        spacingN = aImage.getPixelSpacing(4);
        spacingM = aImage.getPixelSpacing(5);

        positionPatient[0] = aImage.getPosition(0);
        positionPatient[1] = aImage.getPosition(1);
        positionPatient[2] = aImage.getPosition(2);

        orientationPatient[0][0] = aImage.getOrient3D(0, 0);
        orientationPatient[0][1] = aImage.getOrient3D(0, 1);
        orientationPatient[0][2] = aImage.getOrient3D(0, 2);

        orientationPatient[1][0] = aImage.getOrient3D(1, 0);
        orientationPatient[1][1] = aImage.getOrient3D(1, 1);
        orientationPatient[1][2] = aImage.getOrient3D(1, 2);

        orientationPatient[2][0] = aImage.getOrient3D(2, 0);
        orientationPatient[2][1] = aImage.getOrient3D(2, 1);
        orientationPatient[2][2] = aImage.getOrient3D(2, 2);
    }

    ~ImageInfo() {}

    void initialize()
    {
        sizeX = 1;
        sizeY = 1;
        sizeZ = 1;
        sizeT = 1;
        sizeN = 1;
        sizeM = 1;

        spacingX = 1.0;
        spacingY = 1.0;
        spacingZ = 1.0;
        spacingT = 1.0;
        spacingN = 1.0;
        spacingM = 1.0;

        positionPatient[0] = 0.0;
        positionPatient[1] = 0.0;
        positionPatient[2] = 0.0;

        orientationPatient[0][0] = 1.0;
        orientationPatient[0][1] = 0.0;
        orientationPatient[0][2] = 0.0;

        orientationPatient[1][0] = 0.0;
        orientationPatient[1][1] = 1.0;
        orientationPatient[1][2] = 0.0;

        orientationPatient[2][0] = 0.0;
        orientationPatient[2][1] = 0.0;
        orientationPatient[2][2] = 1.0;

        int ind = 0;
        fieldnames.resize(numOfFields);
        fieldnames[ind++] = "sizeX";
        fieldnames[ind++] = "sizeY";
        fieldnames[ind++] = "sizeZ";
        fieldnames[ind++] = "sizeT";
        fieldnames[ind++] = "sizeN";
        fieldnames[ind++] = "sizeM";
        fieldnames[ind++] = "spacingX";
        fieldnames[ind++] = "spacingY";
        fieldnames[ind++] = "spacingZ";
        fieldnames[ind++] = "spacingT";
        fieldnames[ind++] = "spacingN";
        fieldnames[ind++] = "spacingM";
        fieldnames[ind++] = "positionPatient";
        fieldnames[ind++] = "orientationPatient";
        fieldnames[ind++] = "xsize";
        fieldnames[ind++] = "ysize";
        fieldnames[ind++] = "zsize";
        fieldnames[ind++] = "xvoxelsize";
        fieldnames[ind++] = "yvoxelsize";
        fieldnames[ind++] = "zvoxelsize";
    }

    mxArray* convertToMatlab() const 
    {
        try
        {
            mwSize num[2] = {1, 1};
            mxArray* info = mxCreateStructArray(2, num, numOfFields, const_cast<const char**>(&fieldnames[0]));

            int ind = 0;

            mxSetField(info, 0, fieldnames[ind++], mxCreateDoubleScalar(sizeX));
            mxSetField(info, 0, fieldnames[ind++], mxCreateDoubleScalar(sizeY));
            mxSetField(info, 0, fieldnames[ind++], mxCreateDoubleScalar(sizeZ));
            mxSetField(info, 0, fieldnames[ind++], mxCreateDoubleScalar(sizeT));
            mxSetField(info, 0, fieldnames[ind++], mxCreateDoubleScalar(sizeN));
            mxSetField(info, 0, fieldnames[ind++], mxCreateDoubleScalar(sizeM));

            mxSetField(info, 0, fieldnames[ind++], mxCreateDoubleScalar(spacingX));
            mxSetField(info, 0, fieldnames[ind++], mxCreateDoubleScalar(spacingY));
            mxSetField(info, 0, fieldnames[ind++], mxCreateDoubleScalar(spacingZ));
            mxSetField(info, 0, fieldnames[ind++], mxCreateDoubleScalar(spacingT));
            mxSetField(info, 0, fieldnames[ind++], mxCreateDoubleScalar(spacingN));
            mxSetField(info, 0, fieldnames[ind++], mxCreateDoubleScalar(spacingM));

            mxArray* mxPositionPatient = mxCreateDoubleMatrix(1, 3, mxREAL);
            double* pPositionData = mxGetPr(mxPositionPatient);
            pPositionData[0] = positionPatient[0];
            pPositionData[1] = positionPatient[1];
            pPositionData[2] = positionPatient[2];
            mxSetField(info, 0, fieldnames[ind++], mxPositionPatient);

            mxArray* mxOrientationPatient = mxCreateDoubleMatrix(3, 3, mxREAL);
            double* pOrientationData = mxGetPr(mxOrientationPatient);
            pOrientationData[0] = orientationPatient[0][0];
            pOrientationData[1] = orientationPatient[1][0];
            pOrientationData[2] = orientationPatient[2][0];
            pOrientationData[3] = orientationPatient[0][1];
            pOrientationData[4] = orientationPatient[1][1];
            pOrientationData[5] = orientationPatient[2][1];
            pOrientationData[6] = orientationPatient[0][2];
            pOrientationData[7] = orientationPatient[1][2];
            pOrientationData[8] = orientationPatient[2][2];
            mxSetField(info, 0, fieldnames[ind++], mxOrientationPatient);

            mxSetField(info, 0, fieldnames[ind++], mxCreateDoubleScalar(sizeX));
            mxSetField(info, 0, fieldnames[ind++], mxCreateDoubleScalar(sizeY));
            mxSetField(info, 0, fieldnames[ind++], mxCreateDoubleScalar(sizeZ));

            mxSetField(info, 0, fieldnames[ind++], mxCreateDoubleScalar(spacingX));
            mxSetField(info, 0, fieldnames[ind++], mxCreateDoubleScalar(spacingY));
            mxSetField(info, 0, fieldnames[ind++], mxCreateDoubleScalar(spacingZ));

            return info;
        }
        catch(...)
        {
            mexErrMsgTxt("Exceptions happened in ImageInfo::convertToMatlab() ... ");
            throw;
        }

        return NULL;
    }

    bool convertFromMatlab(const mxArray* info)
    {
        try
        {
            int ind = 0;
            sizeX = static_cast<int>(mxGetScalar(mxGetField(info, 0, fieldnames[ind++])));
            sizeY = static_cast<int>(mxGetScalar(mxGetField(info, 0, fieldnames[ind++])));
            sizeZ = static_cast<int>(mxGetScalar(mxGetField(info, 0, fieldnames[ind++])));
            sizeT = static_cast<int>(mxGetScalar(mxGetField(info, 0, fieldnames[ind++])));
            sizeN = static_cast<int>(mxGetScalar(mxGetField(info, 0, fieldnames[ind++])));
            sizeM = static_cast<int>(mxGetScalar(mxGetField(info, 0, fieldnames[ind++])));

            spacingX = mxGetScalar(mxGetField(info, 0, fieldnames[ind++]));
            spacingY = mxGetScalar(mxGetField(info, 0, fieldnames[ind++]));
            spacingZ = mxGetScalar(mxGetField(info, 0, fieldnames[ind++]));
            spacingT = mxGetScalar(mxGetField(info, 0, fieldnames[ind++]));
            spacingN = mxGetScalar(mxGetField(info, 0, fieldnames[ind++]));
            spacingM = mxGetScalar(mxGetField(info, 0, fieldnames[ind++]));

            mxArray* mxPositionPatient = mxGetField(info, 0, fieldnames[ind++]);
            double* pPositionData = mxGetPr(mxPositionPatient);
            positionPatient[0] = pPositionData[0];
            positionPatient[1] = pPositionData[1];
            positionPatient[2] = pPositionData[2];
            
            mxArray* mxOrientationPatient = mxGetField(info, 0, fieldnames[ind++]);
            double* pOrientationData = mxGetPr(mxOrientationPatient);
            orientationPatient[0][0] = pOrientationData[0];
            orientationPatient[1][0] = pOrientationData[1];
            orientationPatient[2][0] = pOrientationData[2];
            orientationPatient[0][1] = pOrientationData[3];
            orientationPatient[1][1] = pOrientationData[4];
            orientationPatient[2][1] = pOrientationData[5];
            orientationPatient[0][2] = pOrientationData[6];
            orientationPatient[1][2] = pOrientationData[7];
            orientationPatient[2][2] = pOrientationData[8];
        }
        catch(...)
        {
            mexErrMsgTxt("Exceptions happened in ImageInfo::convertFromMatlab() ... ");
            return false;
        }

        return true;
    }

    Size<IndexType, 2> getSize2D() const { return Size<IndexType, 2>(sizeX, sizeY); }
    Spacing<2> getSpacing2D() const { return Spacing<2>(spacingX, spacingY); }

    Size<IndexType, 3> getSize3D() const { return Size<IndexType, 3>(sizeX, sizeY, sizeZ); }
    Spacing<3> getSpacing3D() const { return Spacing<3>(spacingX, spacingY, spacingZ); }

    Size<IndexType, 4> getSize4D() const 
    { 
        Size<IndexType, 4> aSize; 
        aSize[0] = sizeX;
        aSize[1] = sizeY;
        aSize[2] = sizeZ;
        aSize[3] = sizeT;
        return aSize;
    }

    Size<IndexType, 5> getSize5D() const 
    { 
        Size<IndexType, 5> aSize; 
        aSize[0] = sizeX;
        aSize[1] = sizeY;
        aSize[2] = sizeZ;
        aSize[3] = sizeT;
        aSize[4] = sizeN;
        return aSize;
    }

    Size<IndexType, 6> getSize6D() const 
    { 
        Size<IndexType, 6> aSize; 
        aSize[0] = sizeX;
        aSize[1] = sizeY;
        aSize[2] = sizeZ;
        aSize[3] = sizeT;
        aSize[4] = sizeN;
        aSize[5] = sizeM;
        return aSize;
    }

    Spacing<4> getSpacing4D() const 
    { 
        Spacing<4> spacing;
        spacing[0] = spacingX;
        spacing[1] = spacingY;
        spacing[2] = spacingZ;
        spacing[3] = spacingT;

        return spacing; 
    }

    Spacing<5> getSpacing5D() const 
    { 
        Spacing<5> spacing;
        spacing[0] = spacingX;
        spacing[1] = spacingY;
        spacing[2] = spacingZ;
        spacing[3] = spacingT;
        spacing[4] = spacingN;

        return spacing; 
    }

    Spacing<6> getSpacing6D() const 
    { 
        Spacing<6> spacing;
        spacing[0] = spacingX;
        spacing[1] = spacingY;
        spacing[2] = spacingZ;
        spacing[3] = spacingT;
        spacing[4] = spacingN;
        spacing[5] = spacingM;

        return spacing; 
    }

    Point3d<double> getPosition() const { return Point3d<double>(positionPatient[0], positionPatient[1], positionPatient[2]); }
    Point3d<double> getOrient3D(int i) const { return Point3d<double>(orientationPatient[i][0], orientationPatient[i][1], orientationPatient[i][2]); }

protected: 

    std::vector<char*> fieldnames;
};

template <typename ValueType> 
class FtkMatlabConverterBase : public Object
{
public:
    
    /** @name Typedefs */
    //@{
    /// 2D Image 
    typedef Image<ValueType, 2> Image2DType;
    /// 3D Image 
    typedef Image<ValueType, 3> Image3DType;
    /// 4D Image 
    typedef Image<ValueType, 4> Image4DType;
    /// 5D Image 
    typedef Image<ValueType, 5> Image5DType;
    /// 6D Image 
    typedef Image<ValueType, 6> Image6DType;
    /// Image array 
    typedef ImageContainerArray<Image2DType> ImageContainerArrayType;
    /// Image matrix
    typedef ImageContainerMatrix<ImageContainerArrayType> ImageContainerMatrixType;
    /// std vector type
    typedef std::vector<ValueType> StdVectorType;
    /// ftk vector type
    typedef Array1d<ValueType> FtkVectorType;
    /// ftk matrix type
    typedef RealMatrix<ValueType> FtkMatrixType;
    //@}

    /** @name Constructors and destructor */
    //@{
    FtkMatlabConverterBase() {}
    virtual ~FtkMatlabConverterBase() {}
    //@}

    /** @name functions to convert ftk to/from Matlab */
    //@{
    // 2D image
    virtual bool convertToMatlab(const Image2DType& aImage, mxArray*& aMxImage, mxArray*& aHeader) = 0;
    virtual bool convertFromMatlab(const mxArray* aMxImage, const mxArray* aHeader, Image2DType& aImage) = 0;
    // 3D image
    virtual bool convertToMatlab(const Image3DType& aImage, mxArray*& aMxImage, mxArray*& aHeader) = 0;
    virtual bool convertFromMatlab(const mxArray* aMxImage, const mxArray* aHeader, Image3DType& aImage) = 0;
    // image array
    virtual bool convertToMatlab(const ImageContainerArrayType& aImageArray, mxArray*& aMxImage, mxArray*& aHeader) = 0;
    virtual bool convertFromMatlab(const mxArray* aMxImage, const mxArray* aHeader, ImageContainerArrayType& aImage) = 0;
    // image matrix
    virtual bool convertToMatlab(const ImageContainerMatrixType& aImageMatrix, mxArray*& aMxImage, mxArray*& aHeader) = 0;
    virtual bool convertFromMatlab(const mxArray* aMxImage, const mxArray* aHeader, ImageContainerMatrixType& aImageMatrix) = 0;
    // ftk vector
    virtual bool convertToMatlab(const FtkVectorType& vec, mxArray*& aMxArray) = 0;
    virtual bool convertFromMatlab(const mxArray* aMxArray, FtkVectorType& vec) = 0;
    // ftk matrix
    virtual bool convertToMatlab(const FtkMatrixType& vec, mxArray*& aMxArray) = 0;
    virtual bool convertFromMatlab(const mxArray* aMxArray, FtkMatrixType& vec) = 0;

    // std vector
    virtual bool convertToMatlab(const StdVectorType& vec, mxArray*& aMxArray) = 0;
    virtual bool convertFromMatlab(const mxArray* aMxArray, StdVectorType& vec) = 0;
    // std string
    virtual bool convertToMatlab(const std::string& str, mxArray*& aMxStr);
    virtual bool convertFromMatlab(const mxArray* aMxStr, std::string& str);
    //@}

    virtual void print(std::ostream& os) const = 0;

protected:
    
};

// -------------------------------------------------------
// std string
// -------------------------------------------------------

template <typename ValueType> 
bool FtkMatlabConverterBase<ValueType>::
convertToMatlab(const std::string& str, mxArray*& aMxStr)
{
    aMxStr = mxCreateString(str.c_str());
    return (aMxStr != NULL);
}

template <typename ValueType> 
bool FtkMatlabConverterBase<ValueType>::
convertFromMatlab(const mxArray* aMxStr, std::string& str)
{
    FTK_CHECK_RETURN_FALSE(aMxStr!=NULL);

    int buflen = mxGetNumberOfElements(aMxStr) + 1;

    std::vector<char> buf(buflen, '\0');

    if (mxGetString(aMxStr, &buf[0], buflen) != 0)
        return false;

    str = std::string(&buf[0]);

    return true;
}

END_NAMESPACE_1(ftk)

#endif // FTK_FTKMATLABCONVERTERBASE_H 

#endif // FTK_MATLAB_SUPPORT 
