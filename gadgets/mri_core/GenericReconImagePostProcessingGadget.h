/** \file   GenericReconImagePostProcessingGadget.h
    \brief  This gadget can serve as the base class for the image stage reconstruction compuation
            Examples include e.g. fat-water seperation, motion correction etc.
    \author Hui Xue
*/

#pragma once

#include <complex>
#include "gadgetron_mricore_export.h"
#include "Gadget.h"
#include "GadgetronTimer.h"

#include "ismrmrd/ismrmrd.h"
#include "ismrmrd/xml.h"
#include "ismrmrd/meta.h"

#include "mri_core_def.h"
#include "mri_core_data.h"
#include "mri_core_utility.h"

// #include "gtPlusIOAnalyze.h"

#include "hoNDObjectArray.h"
#include "hoNDImageContainer2D.h"

namespace Gadgetron { 

// the dimensionsal order of buffered images
// [Cha Slice Con Phase Rep Set Ave]
//   0    1    2   3     4   5   6

template <typename T, int D> 
class EXPORTGADGETSMRICORE GenericReconImagePostProcessingGadget : public Gadget1< hoNDObjectArray< hoNDImage<std::complex<T>, D> > >
{
public:
    GADGET_DECLARE(GenericReconImagePostProcessingGadget);

    typedef std::complex<T> ValueType;
    typedef hoNDImage<ValueType, D> ImageType;
    typedef hoNDArray<ValueType> ImgArrayType;
    typedef hoNDObjectArray<ImageType> ImageBufferType;
    typedef hoNDImageContainer2D<ImageType> ImageContainerType;

    typedef hoNDImage<T, D> ImageMagType;
    typedef hoNDImageContainer2D<ImageMagType> ImageContainerMagType;

    typedef Gadget1< ImageBufferType > BaseClass;

    GenericReconImagePostProcessingGadget();
    virtual ~GenericReconImagePostProcessingGadget();

    virtual int close(unsigned long flags);

    /// debug and timing
    GADGET_PROPERTY(verbose, bool, "Whether to print more information", false);
    GADGET_PROPERTY(debug_folder, std::string, "If set, the debug output will be written out", "");
    GADGET_PROPERTY(perform_timing, bool, "Whether to perform timing on some computational steps", false);

    /// whether to send out functional maps
    GADGET_PROPERTY(send_out_gfactor_map, bool, "Whether to send out gfactor maps", false);
    GADGET_PROPERTY(send_out_snr_map, bool, "Whether to send out snr maps", true);
    GADGET_PROPERTY(send_out_std_map, bool, "Whether to send out std maps", false);

    GADGET_PROPERTY(processing_result_image_series_num, int, "Starting image series number for processing results", 100);

protected:

    virtual int process_config(ACE_Message_Block* mb);
    virtual int process(GadgetContainerMessage<ImageBufferType>* m1);
    virtual int process_image_buffer(ImageBufferType& ori);

    /// compute the image number
    size_t compute_image_number (ISMRMRD::ImageHeader& imheader, size_t nCHA, size_t cha);

    /// release the image buffer
    int release_image_buffer(ImageBufferType& buf)
    {
        try
        {
            size_t N = buf.get_number_of_elements();
            size_t ii;
            for (ii = 0; ii<N; ii++)
            {
                ImageType* pImage = buf(ii);
                if (buf.delete_data_on_destruct() && pImage != NULL)
                {
                    delete pImage;
                    buf(ii) = NULL;
                }
            }
        }
        catch (...)
        {
            GERROR_STREAM("Errors happened in GenericReconImagePostProcessingGadget::release_image_buffer(...) ... ");
            return GADGET_FAIL;
        }

        return GADGET_OK;
    }

    /// send out the images
    virtual int send_out_images(ImageBufferType& images, int seriesNum, const std::vector<std::string>& processStr, const std::vector<std::string>& dataRole, const std::vector<float>& windowCenter=std::vector<float>(), const std::vector<float>& windowWidth=std::vector<float>());

    /// store gfactor map, the gmap will be stored if it is for a new slice
    int store_gfactor_map(const ImageType& gmap);

    /// find gfactor map for certain slice
    /// if cannot find the gmap, return GADGET_FAIL
    int find_gfactor_map(size_t slc, ImageType& gmap);

    // number of times the process function is called
    size_t process_called_times_;

    // encoding space size
    ISMRMRD::EncodingCounters meas_max_idx_;

    // the buffer to store incoming gfactor map for every slice
    std::vector< ImageType > gfactor_buf_;

    // --------------------------------------------------
    // function and variables for debug and timing
    // --------------------------------------------------

    /*void export_image_container(const hoNDImageContainer2D< hoNDImage<std::complex<T>, D> >& input, const std::string& name)
    {
        if (!this->debug_folder_full_path_.empty())
        {
            size_t R = input.rows();

            size_t r;

            hoNDArray<ValueType> outArray;

            for (r = 0; r<R; r++)
            {
                input.to_NDArray(r, outArray);

                std::ostringstream ostr;
                ostr << name << "_" << r;

                if (!this->debug_folder_full_path_.empty()) gt_exporter_.exportArrayComplex(outArray, this->debug_folder_full_path_ + ostr.str());
            }
        }
    }*/

    /*void export_image_container(const hoNDImageContainer2D< hoNDImage<T, D> >& input, const std::string& name)
    {
        if (!this->debug_folder_full_path_.empty())
        {
            size_t R = input.rows();

            size_t r;

            hoNDArray<T> outArray;

            for (r = 0; r<R; r++)
            {
                input.to_NDArray(r, outArray);

                std::ostringstream ostr;
                ostr << name << "_" << r;

                if (!this->debug_folder_full_path_.empty()) gt_exporter_.exportArray(outArray, this->debug_folder_full_path_ + ostr.str());
            }
        }
    }*/

    // clock for timing
    Gadgetron::GadgetronTimer gt_timer_local_;
    Gadgetron::GadgetronTimer gt_timer_;

    // debug folder
    // std::string debug_folder_full_path_;

    // exporter
    // Gadgetron::gtPlus::gtPlusIOAnalyze gt_exporter_;
};

class EXPORTGADGETSMRICORE GenericReconImagePostProcessing2DGadget : public GenericReconImagePostProcessingGadget<float, 2>
{
public:
    GADGET_DECLARE(GenericReconImagePostProcessing2DGadget);

    typedef GenericReconImagePostProcessingGadget<float, 2> BaseClass;

    GenericReconImagePostProcessing2DGadget();
    virtual ~GenericReconImagePostProcessing2DGadget();
};

class EXPORTGADGETSMRICORE GenericReconImagePostProcessing3DGadget : public GenericReconImagePostProcessingGadget<float, 3>
{
public:
    GADGET_DECLARE(GenericReconImagePostProcessing3DGadget);

    typedef GenericReconImagePostProcessingGadget<float, 3> BaseClass;

    GenericReconImagePostProcessing3DGadget();
    virtual~GenericReconImagePostProcessing3DGadget();
};

}
