/** \file   GtPlusImageReconGadget.h
    \brief  The GtPlus image reconstruction gadget, used after GtPlus kspace reconstruction
    \author Hui Xue
*/

#pragma once

#include <complex>
#include "GtPlusGadgetExport.h"
#include "Gadget.h"
#include "hoNDArray.h"
#include "hoNDObjectArray.h"
#include "ismrmrd.h"
#include "GadgetIsmrmrdReadWrite.h"

#include "hoNDImageContainer2D.h"
#include "hoNDArray_utils.h"
#include "hoNDImage.h"

#include "GadgetronCommon.h"
#include "GtPlusGadgetImageArray.h"

#include "gtPlusIOAnalyze.h"
#include "gtPlusISMRMRDReconUtil.h"

#include "GadgetStreamController.h"
#include "GtPlusReconGadgetUtil.h"

#ifdef USE_OMP
    #include "omp.h"
#endif // USE_OMP

namespace Gadgetron
{

// the dimensionsal order of buffered images
// [Cha Slice E2 Con Phase Rep Set Ave]
//   0    1    2   3   4    5   6   7

class EXPORTGTPLUSGADGET GtPlusImageReconGadget : public Gadget1< hoNDObjectArray< hoNDImage<std::complex<float>, 2> > >
{
public:
    GADGET_DECLARE(GtPlusImageReconGadget);

    typedef float T;
    typedef std::complex<T> ValueType;

    typedef hoNDImage<ValueType, 2> ImageType;
    typedef hoNDImage<ValueType, 3> Image3DType;

    typedef hoNDImage<T, 2> ImageMagType;
    typedef hoNDImage<T, 3> Image3DMagType;

    // typedef hoNDArray<ImageType*> ImageBufferType;
    typedef hoNDObjectArray<ImageType> ImageBufferType;
    typedef hoNDArray<ValueType> ImgArrayType;

    typedef hoNDImageContainer2D<ImageType> ImageContainer2DType;
    typedef hoNDImageContainer2D<Image3DType> ImageContainer3DType;

    typedef hoNDImageContainer2D<ImageMagType> ImageContainer2DMagType;
    typedef hoNDImageContainer2D<Image3DMagType> ImageContainer3DMagType;

    typedef Gadget1< ImageBufferType > BaseClass;

    GtPlusImageReconGadget();
    ~GtPlusImageReconGadget();

    virtual int close(unsigned long flags);

    /// image series number
    int image_series_num_;

    // debug folder
    std::string debugFolder_;
    std::string debugFolder_fullPath_;

    // whether to perform timing
    bool performTiming_;

protected:

    // encoding space size
    ISMRMRD::EncodingCounters meas_max_idx_;

    // read in parameters
    bool readParameters();

    virtual int process_config(ACE_Message_Block* mb);
    virtual int process(GadgetContainerMessage<ImageBufferType>* m1);

    virtual int processImageBuffer(ImageBufferType& ori);

    /// fill image buffer with null point
    bool fillWithNULL(ImageBufferType& buf);

    /// release the image buffer
    bool releaseImageBuffer(ImageBufferType& buf);

    /// get the 3D images from 2D buffer
    Image3DType* getImage3DFromImage2D(ImageBufferType& ori, size_t cha, size_t slc, size_t con, size_t phs, size_t rep, size_t set, size_t ave);

    /// get the 2D image in buffer from a 3D image
    bool getImage2DFromImage3D(Image3DType& image3D, ImageBufferType& image2DBuf);

    /// compute the image number
    size_t computeSeriesImageNumber (ISMRMRD::ImageHeader& imheader, size_t nCHA, size_t cha, size_t nE2, size_t e2);

    /// send out the images as a Gadget3 message
    /// ISMRMRD::ImageHeader, hoNDArray< std::complex<float> >, GtImageAttribType
    /// windowCenter and windowWidth is for every SLC
    virtual bool sendOutImages(ImageBufferType& images, int seriesNum, const std::vector<std::string>& processStr, const std::vector<std::string>& dataRole, const std::vector<float>& windowCenter=std::vector<float>(), const std::vector<float>& windowWidth=std::vector<float>());

    /// utility function to export image container
    bool exportImageContainer2D(ImageContainer2DType& input, const std::string& prefix);
    bool exportImageContainer2D(ImageContainer2DMagType& input, const std::string& prefix);

    bool exportImageContainer3D(ImageContainer3DType& input, const std::string& prefix);
    bool exportImageContainer3D(ImageContainer3DMagType& input, const std::string& prefix);

    // util for gtplus
    Gadgetron::gtPlus::gtPlusISMRMRDReconUtil<GT_Complex8> gtPlus_util_;

    // clock for timing
    Gadgetron::GadgetronTimer gt_timer1_;
    Gadgetron::GadgetronTimer gt_timer2_;
    Gadgetron::GadgetronTimer gt_timer3_;

    // exporter
    Gadgetron::gtPlus::gtPlusIOAnalyze gt_exporter_;

    // in verbose mode, more info is printed out
    bool verboseMode_;
};

}
