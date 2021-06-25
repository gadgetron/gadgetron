/** \file   GenericImageReconGadget.h
    \brief  The GtPlus image reconstruction gadget, used after GtPlus kspace reconstruction
    \author Hui Xue
*/

#pragma once

#include <complex>
#include "gadgetron_mricore_export.h"

#include "ismrmrd/ismrmrd.h"
#include "ismrmrd/xml.h"
#include "ismrmrd/meta.h"

#include "Gadget.h"
#include "hoNDArray.h"
#include "ismrmrd/ismrmrd.h"

#include "hoNDImageContainer2D.h"
#include "hoNDArray_utils.h"
#include "hoNDArray_reductions.h"
#include "hoMRImage.h"
#include "hoNDObjectArray.h"
#include "mri_core_def.h"
#include "mri_core_kspace_filter.h"
#include "mri_core_data.h"
#include "mri_core_utility.h"
#include "ImageIOAnalyze.h"
#include "GadgetronTimer.h"

#ifdef USE_OMP
#include "omp.h"
#endif // USE_OMP

namespace Gadgetron {

    // the dimensionsal order of buffered images
    // [Cha Slice Con Phase Rep Set Ave]
    //   0    1    2   3     4   5   6

    class EXPORTGADGETSMRICORE GenericImageReconGadgetBase
    {
    public:
        typedef float T;
        typedef std::complex<T> ValueType;

        typedef hoMRImage<ValueType, 3> Image3DType;
        typedef hoMRImage<ValueType, 2> Image2DType;

        typedef hoMRImage<T, 3> Image3DMagType;
        typedef hoMRImage<T, 2> Image2DMagType;

        typedef hoNDObjectArray<Image3DType> Image3DBufferType;
        typedef hoNDObjectArray<Image2DType> Image2DBufferType;

        typedef hoNDObjectArray<Image3DMagType> Image3DMagBufferType;
        typedef hoNDObjectArray<Image2DMagType> Image2DMagBufferType;

        typedef hoNDArray<ValueType> ImgArrayType;
        typedef hoNDArray<T> ImgArrayMagType;

        typedef hoNDArray<uint16_t> ImgRGBType;

        typedef hoNDImageContainer2D<Image2DType> ImageContainer2DType;
        typedef hoNDImageContainer2D<Image3DType> ImageContainer3DType;

        typedef hoNDImageContainer2D<Image2DMagType> ImageContainer2DMagType;
        typedef hoNDImageContainer2D<Image3DMagType> ImageContainer3DMagType;

        typedef hoMRImage< std::complex<double>, 2> ImageDoubleType;
        typedef hoNDImageContainer2D<ImageDoubleType> ImageContainer2DDoubleType;
        typedef hoNDImageContainer2D< hoMRImage< double, 2> > ImageContainer2DMagDoubleType;

        typedef hoNDImageContainer2D< hoNDImage<ValueType, 2> > NDImageContainer2DType;
        typedef hoNDImageContainer2D< hoNDImage<ValueType, 3> > NDImageContainer3DType;

        typedef hoNDImageContainer2D< hoNDImage<T, 2> > NDImageContainer2DMagType;
        typedef hoNDImageContainer2D< hoNDImage<T, 2> > NDImageContainer3DMagType;

        typedef hoNDImageContainer2D< hoNDImage< std::complex<double>, 2> > NDImageContainer2DDoubleType;
        typedef hoNDImageContainer2D< hoNDImage<double, 2> > NDImageContainer2DMagDoubleType;

        typedef Gadget1< Image2DBufferType > BaseClass;

        GenericImageReconGadgetBase();
        virtual ~GenericImageReconGadgetBase();

        /// get the 3D images from 2D buffer
        Image3DType* getImage3DFromImage2D(Image2DBufferType& ori, size_t cha, size_t slc, size_t con, size_t phs, size_t rep, size_t set, size_t ave);

        /// get the 2D image in buffer from a 3D image
        bool getImage2DFromImage3D(Image3DType& image3D, Image2DBufferType& image2DBuf);

        /// compute the image number
        size_t computeSeriesImageNumber(ISMRMRD::ISMRMRD_ImageHeader& imheader, size_t nCHA, size_t cha, size_t nE2, size_t e2);

        virtual void decorateImageHeader(ISMRMRD::ISMRMRD_ImageHeader& header, ISMRMRD::MetaContainer& attrib, int seriesNum, const std::vector<std::string>& processStr, const std::vector<std::string>& dataRole, const std::vector<float>& windowCenter, const std::vector<float>& windowWidth, bool resetImageCommentsParametricMaps, size_t slc, size_t SLC);

        /// convert 3D container to 2D and vice versa
        /// a [RO E1 E2] 3D image will be converted into E2 2D images
        /// All images in a row of the container must have the same size
        /// e2[r] gives the original E2 dimension of 3D images for the row r
        bool convertContainer3Dto2D(const ImageContainer3DType& input, std::vector<size_t>& e2, ImageContainer2DType& output);
        /// every e2[r] 2D images are converted into a 3D images
        bool convertContainer2Dto3D(const ImageContainer2DType& input, const std::vector<size_t>& e2, ImageContainer3DType& output);

        /// if E2 in 3D image equals 1, cast it to 2D image
        //void cast3DImageTo2D(const Image3DType& im3D, Image2DType& im2D);
        //void cast2DImageTo3D(const Image2DType& im2D, Image3DType& im3D);

        /// utility function to export image container
        template<typename ImageContainerType>
        bool exportImageContainer(ImageContainerType& input, const std::string& prefix)
        {
            if (!debug_folder_full_path_.empty())
            {
                size_t R = input.rows();

                size_t r;

                hoNDArray<ValueType> outArray;

                for (r = 0; r < R; r++)
                {
                    input.to_NDArray(r, outArray);

                    std::ostringstream ostr;
                    ostr << prefix << "_" << r;

                    if (!debug_folder_full_path_.empty()) gt_exporter_.export_array_complex(outArray, debug_folder_full_path_ + ostr.str());
                }
            }

            return true;
        }

        bool exportImageContainer2D(ImageContainer2DType& input, const std::string& prefix);
        bool exportImageContainer2D(ImageContainer2DMagType& input, const std::string& prefix);

        bool exportImageContainer2D(ImageContainer2DDoubleType& input, const std::string& prefix);
        bool exportImageContainer2D(ImageContainer2DMagDoubleType& input, const std::string& prefix);

        bool exportImageContainer3D(ImageContainer3DType& input, const std::string& prefix);
        bool exportImageContainer3D(ImageContainer3DMagType& input, const std::string& prefix);

        bool exportImageContainer2D(NDImageContainer2DType& input, const std::string& prefix);
        bool exportImageContainer2D(NDImageContainer2DMagType& input, const std::string& prefix);

        bool exportImageContainer2D(NDImageContainer2DDoubleType& input, const std::string& prefix);
        bool exportImageContainer2D(NDImageContainer2DMagDoubleType& input, const std::string& prefix);

        bool exportImageContainer3D(NDImageContainer3DType& input, const std::string& prefix);
        bool exportImageContainer3D(NDImageContainer3DMagType& input, const std::string& prefix);

        /// set pixel value range of images
        /// if ignore_empty_image == true and im is empty, do not do the filling
        bool set_pixel_value_range(hoNDArray<ValueType>& im, T min_V, T max_V, bool ignore_empty_image = true);
        bool set_pixel_value_range(hoNDArray<T>& im, T min_V, T max_V, bool ignore_empty_image = true);

        /// store gfactor map
        /// the gmap will be stored if it is for a new slice or e2
        bool storeGFactorMap(const Image2DType& gmap);
        bool storeGFactorMap(const Image3DType& gmap);

        /// find gfactor map for certain slice and E2
        /// true means the gmap is found; false means the gmap is not found
        bool findGFactorMap(size_t slc, Image2DType& gmap);
        bool findGFactorMap(size_t slc, Image3DType& gmap);

        /// image series number
        int image_series_num_;

        // encoding space size
        ISMRMRD::EncodingCounters meas_max_idx_;

        // clock for timing
        Gadgetron::GadgetronTimer gt_timer_local_;
        Gadgetron::GadgetronTimer gt_timer_;

        // debug folder
        std::string debug_folder_full_path_;

        // exporter
        Gadgetron::ImageIOAnalyze gt_exporter_;

        // protocol encoding and matrix size
        size_t matrix_size_encoding_[3];
        size_t matrix_size_recon_[3];

        // forward gfactor map
        bool send_out_gfactor_map_;
        bool send_out_snr_map_;
        bool send_out_std_map_;
        bool send_out_wrap_around_map_;

        // whether interpolation is on
        bool inteprolation_on_;

        // whether to send out PSIR images as real type
        bool send_out_PSIR_as_real_;

        // the buffer to store incoming gfactor map for every slice
        std::vector< Image2DType > gfactor_buf_;
        std::vector< Image3DType > gfactor_buf_3D_;
    };

    class EXPORTGADGETSMRICORE GenericImageReconGadget : public BasicPropertyGadget, public GenericImageReconGadgetBase
    {
    public:
        GADGET_DECLARE(GenericImageReconGadget);

        typedef BasicPropertyGadget BaseClass;

        GenericImageReconGadget();
        ~GenericImageReconGadget();

        virtual int close(unsigned long flags);

        GADGET_PROPERTY(image_series_num, int, "Image series number", 100);
        GADGET_PROPERTY(add_original_series_num, bool, "If true, the original image series index is added to the set seriesNum", false);

        GADGET_PROPERTY(verbose, bool, "Whether to print more information", false);
        GADGET_PROPERTY(debug_folder, std::string, "If set, the debug output will be written out", "");
        GADGET_PROPERTY(perform_timing, bool, "Whether to perform timing on some computational steps", false);

        GADGET_PROPERTY(send_out_image_array, bool, "Whether to send out image array; if false, send out individual image", false);
        GADGET_PROPERTY(send_out_gfactor_map, bool, "Whether to send out gfactor maps", false);
        GADGET_PROPERTY(send_out_snr_map, bool, "Whether to send out snr maps", true);
        GADGET_PROPERTY(send_out_std_map, bool, "Whether to send out std maps", false);
        GADGET_PROPERTY(send_out_wrap_around_map, bool, "Whether to send out wrap around maps", false);
        GADGET_PROPERTY(send_out_PSIR_as_real, bool, "Whether to send out PSIR image as real type", false);

        GADGET_PROPERTY(timeStampResolution, float, "Time tick resolution in second", 0.0025f);

    protected:

        // read in parameters
        bool readParameters();

        virtual int process_config(ACE_Message_Block* mb);

        int process(ACE_Message_Block* mb) override;

        virtual int process_image(GadgetContainerMessage< ISMRMRD::ImageHeader >* m1);
        int process2D(GadgetContainerMessage<Image2DBufferType>* m1);
        int process3D(GadgetContainerMessage<Image3DBufferType>* m1);

        virtual int processImageBuffer(Image2DBufferType& ori);
        virtual int processImageBuffer(Image3DBufferType& ori);

        ///// send out the images as a Gadget3 message
        ///// windowCenter and windowWidth is for every SLC
        bool sendOutImages(Image2DBufferType& images, int seriesNum, const std::vector<std::string>& processStr, const std::vector<std::string>& dataRole, const std::vector<float>& windowCenter = std::vector<float>(), const std::vector<float>& windowWidth = std::vector<float>(), bool resetImageCommentsParametricMaps = true, Gadget* anchor = NULL);

        bool sendOutImageBuffer(Image2DBufferType& images, int seriesNum, const std::vector<std::string>& processStr, const std::vector<std::string>& dataRole, const std::vector<float>& windowCenter = std::vector<float>(), const std::vector<float>& windowWidth = std::vector<float>(), bool resetImageCommentsParametricMaps = true, Gadget* anchor = NULL);

        bool sendOutImages(Image3DBufferType& images, int seriesNum, const std::vector<std::string>& processStr, const std::vector<std::string>& dataRole, const std::vector<float>& windowCenter = std::vector<float>(), const std::vector<float>& windowWidth = std::vector<float>(), bool resetImageCommentsParametricMaps = true, Gadget* anchor = NULL);

        bool sendOutImageBuffer(Image3DBufferType& images, int seriesNum, const std::vector<std::string>& processStr, const std::vector<std::string>& dataRole, const std::vector<float>& windowCenter = std::vector<float>(), const std::vector<float>& windowWidth = std::vector<float>(), bool resetImageCommentsParametricMaps = true, Gadget* anchor = NULL);

        bool sendOutRGBImages(Image3DMagBufferType& images, int seriesNum, const std::vector<std::string>& processStr, const std::vector<std::string>& dataRole, bool resetImageCommentsParametricMaps = true, Gadget* anchor = NULL);
        bool sendOutRGBImageBuffer(Image3DMagBufferType& images, int seriesNum, const std::vector<std::string>& processStr, const std::vector<std::string>& dataRole, bool resetImageCommentsParametricMaps = true, Gadget* anchor = NULL);

        bool releaseImageBuffer(hoNDObjectArray< hoMRImage<ValueType, 2> >& buf);
        bool releaseImageBuffer(hoNDObjectArray< hoMRImage<ValueType, 3> >& buf);
    };
}
