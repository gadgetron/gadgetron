
#include "GenericReconFieldOfViewAdjustmentGadget.h"
#include <iomanip>

#include "hoNDArray_utils.h"
#include "hoNDArray_elemwise.h"
#include "hoNDFFT.h"

#include "mri_core_utility.h"
#include "mri_core_kspace_filter.h"
#include "hoNDArray_reductions.h"
#include "mri_core_def.h"

namespace Gadgetron {

    GenericReconFieldOfViewAdjustmentGadget::GenericReconFieldOfViewAdjustmentGadget() : BaseClass()
    {
    }

    GenericReconFieldOfViewAdjustmentGadget::~GenericReconFieldOfViewAdjustmentGadget()
    {
    }

    int GenericReconFieldOfViewAdjustmentGadget::process_config(const mrd::Header& header)
    {
        GADGET_CHECK_RETURN(BaseClass::process_config(header) == GADGET_OK, GADGET_FAIL);

        auto& h = header;

        if (!h.acquisition_system_information)
        {
            GDEBUG("acquisitionSystemInformation not found in header. Bailing out");
            return GADGET_FAIL;
        }

        // -------------------------------------------------

        size_t NE = h.encoding.size();

        num_encoding_spaces_ = NE;

        GDEBUG_CONDITION_STREAM(verbose.value(), "Number of encoding spaces: " << NE);

        // get the encoding FOV and recon FOV

        encoding_FOV_.resize(NE);
        recon_FOV_.resize(NE);
        recon_size_.resize(NE);

        size_t e;
        for (e = 0; e < NE; e++)
        {
            encoding_FOV_[e].resize(3, 0);
            encoding_FOV_[e][0] = h.encoding[e].encoded_space.field_of_view_mm.x;
            encoding_FOV_[e][1] = h.encoding[e].encoded_space.field_of_view_mm.y;
            encoding_FOV_[e][2] = h.encoding[e].encoded_space.field_of_view_mm.z;

            recon_FOV_[e].resize(3, 0);
            recon_FOV_[e][0] = h.encoding[e].recon_space.field_of_view_mm.x;
            recon_FOV_[e][1] = h.encoding[e].recon_space.field_of_view_mm.y;
            recon_FOV_[e][2] = h.encoding[e].recon_space.field_of_view_mm.z;

            recon_size_[e].resize(3, 0);
            recon_size_[e][0] = h.encoding[e].recon_space.matrix_size.x;
            recon_size_[e][1] = h.encoding[e].recon_space.matrix_size.y;
            recon_size_[e][2] = h.encoding[e].recon_space.matrix_size.z;

            GDEBUG_CONDITION_STREAM(verbose.value(), "Encoding space : " << e << " - encoding FOV : [" << encoding_FOV_[e][0] << " " << encoding_FOV_[e][1] << " " << encoding_FOV_[e][2] << " ]");
            GDEBUG_CONDITION_STREAM(verbose.value(), "Encoding space : " << e << " - recon    FOV : [" << recon_FOV_[e][0]    << " " << recon_FOV_[e][1]    << " " << recon_FOV_[e][2] << " ]");
            GDEBUG_CONDITION_STREAM(verbose.value(), "Encoding space : " << e << " - recon    size : [" << recon_size_[e][0] << " " << recon_size_[e][1] << " " << recon_size_[e][2] << " ]");
        }

        return GADGET_OK;
    }

    int GenericReconFieldOfViewAdjustmentGadget::process(Gadgetron::GadgetContainerMessage< mrd::ImageArray >* m1)
    {
        if (perform_timing.value()) { gt_timer_.start("GenericReconFieldOfViewAdjustmentGadget::process"); }

        GDEBUG_CONDITION_STREAM(verbose.value(), "GenericReconFieldOfViewAdjustmentGadget::process(...) starts ... ");

        process_called_times_++;

        mrd::ImageArray* recon_res_ = m1->getObjectPtr();

        // print out recon info
        if (verbose.value())
        {
            GDEBUG_STREAM("----> GenericReconFieldOfViewAdjustmentGadget::process(...) has been called " << process_called_times_ << " times ...");
            std::stringstream os;
            recon_res_->data.print(os);
            GDEBUG_STREAM(os.str());
        }

        if (!debug_folder_full_path_.empty()) { gt_exporter_.export_array_complex(recon_res_->data, debug_folder_full_path_ + "data_before_FOV_adjustment"); }

        // ----------------------------------------------------------
        // FOV adjustment
        // ----------------------------------------------------------

        GADGET_CHECK_RETURN(this->adjust_FOV(*recon_res_) == GADGET_OK, GADGET_FAIL);

        if (!debug_folder_full_path_.empty()) { gt_exporter_.export_array_complex(recon_res_->data, debug_folder_full_path_ + "data_after_FOV_adjustment"); }

        this->gt_streamer_.stream_to_mrd_image_buffer(GENERIC_RECON_STREAM_RECONED_COMPLEX_IMAGE_AFTER_POSTPROCESSING, recon_res_->data, recon_res_->headers, recon_res_->meta);

        GDEBUG_CONDITION_STREAM(verbose.value(), "GenericReconFieldOfViewAdjustmentGadget::process(...) ends ... ");

        // ----------------------------------------------------------
        // send out results
        // ----------------------------------------------------------
        if (this->next()->putq(m1) == -1)
        {
            GERROR("GenericReconFieldOfViewAdjustmentGadget::process, passing data on to next gadget");
            return GADGET_FAIL;
        }

        if (perform_timing.value()) { gt_timer_.stop(); }

        return GADGET_OK;
    }

    void GenericReconFieldOfViewAdjustmentGadget::perform_fft(size_t E2, const hoNDArray< std::complex<float> >& input, hoNDArray< std::complex<float> >& output)
    {
        if (E2>1)
        {
            if (&input == &output)
                Gadgetron::hoNDFFT<float>::instance()->fft3c(output);
            else
                Gadgetron::hoNDFFT<float>::instance()->fft3c(input, output);
        }
        else
        {
            if (&input == &output)
                Gadgetron::hoNDFFT<float>::instance()->fft2c(output);
            else
                Gadgetron::hoNDFFT<float>::instance()->fft2c(input, output);
        }
    }

    void GenericReconFieldOfViewAdjustmentGadget::perform_ifft(size_t E2, const hoNDArray< std::complex<float> >& input, hoNDArray< std::complex<float> >& output)
    {
        if (E2>1)
        {
            if (&input == &output)
                Gadgetron::hoNDFFT<float>::instance()->ifft3c(output);
            else
                Gadgetron::hoNDFFT<float>::instance()->ifft3c(input, output);
        }
        else
        {
            if (&input == &output)
                Gadgetron::hoNDFFT<float>::instance()->ifft2c(output);
            else
                Gadgetron::hoNDFFT<float>::instance()->ifft2c(input, output);
        }
    }

    int GenericReconFieldOfViewAdjustmentGadget::adjust_FOV(mrd::ImageArray& recon_res)
    {
//        try
        {

            size_t RO = recon_res.data.get_size(0);
            size_t E1 = recon_res.data.get_size(1);
            size_t E2 = recon_res.data.get_size(2);

            auto& encoding_fov = recon_res.meta[0]["encoding_FOV"];
            double encodingFOV_RO = std::get<double>(encoding_fov[0]);
            double encodingFOV_E1 = std::get<double>(encoding_fov[1]);
            double encodingFOV_E2 = std::get<double>(encoding_fov[2]);

            auto& recon_fov = recon_res.meta[0]["recon_FOV"];
            double reconFOV_RO = std::get<double>(recon_fov[0]);
            double reconFOV_E1 = std::get<double>(recon_fov[1]);
            double reconFOV_E2 = std::get<double>(recon_fov[2]);

            long encoding = std::get<long>(recon_res.meta[0]["encoding"].front());

            size_t reconSizeRO = recon_size_[encoding][0];
            size_t reconSizeE1 = recon_size_[encoding][1];
            size_t reconSizeE2 = recon_size_[encoding][2];

            // if 2D reconstruction, no need to process along E2
            if (E2 <= 1)
            {
                reconSizeE2 = E2;
                reconFOV_E2 = encodingFOV_E2;
            }

            // if encoded FOV are the same as recon FOV
            if ((std::abs(encodingFOV_RO / 2 - reconFOV_RO)<0.1) && (std::abs(encodingFOV_E1 - reconFOV_E1)<0.1) && (std::abs(encodingFOV_E2 - reconFOV_E2)<0.1))
            {
                if (RO <= reconSizeRO && E1 <= reconSizeE1 && E2 <= reconSizeE2)
                {
                    /** BUG ALERT: This is INCORRECT. The result is stored in `this->res_`, but never used again,
                     * so in the end, `recon_res.data` is NOT adjusted
                     */
                    Gadgetron::zero_pad_resize(recon_res.data, reconSizeRO, reconSizeE1, reconSizeE2, res_);
                }
                else if (RO >= reconSizeRO && E1 >= reconSizeE1 && E2 >= reconSizeE2)
                {
                    this->perform_fft(E2, recon_res.data, kspace_buf_);
                    Gadgetron::crop(reconSizeRO, reconSizeE1, reconSizeE2, kspace_buf_, res_);
                    this->perform_ifft(E2, res_, recon_res.data);
                }
                else
                {
                    GDEBUG_STREAM("Inconsistent image size [" << RO << " " << E1 << " " << E2 << "]; recon image size [" << reconSizeRO << " " << reconSizeE1 << " " << reconSizeE2 << "] ... ");
                    return GADGET_FAIL;
                }
            }
            else if ((encodingFOV_E1 >= reconFOV_E1) && (encodingFOV_E2 >= reconFOV_E2))
            {
                size_t encodingE1 = reconSizeE1;
                if (encodingFOV_E1 > reconFOV_E1)
                {
                    double spacingE1 = reconFOV_E1 / reconSizeE1;
                    encodingE1 = (size_t)2*std::lround(encodingFOV_E1 / (2*spacingE1));
                }

                size_t encodingE2 = reconSizeE2;
                if (encodingFOV_E2 > reconFOV_E2)
                {
                    double spacingE2 = reconFOV_E2 / reconSizeE2;
                    encodingE2 = (size_t)2*std::lround(encodingFOV_E2 / (2*spacingE2));
                }

                hoNDArray< std::complex<float> >* pSrc = &recon_res.data;
                hoNDArray< std::complex<float> >* pDst = &res_;
                hoNDArray< std::complex<float> >* pTmp;

                // adjust E1
                if (encodingE1 >= E1 + 1)
                {
                    Gadgetron::zero_pad_resize(*pSrc, RO, encodingE1, E2, *pDst);
                    pTmp = pSrc; pSrc = pDst; pDst = pTmp;
                }
                else if (encodingE1 <= E1 - 1)
                {
                    this->perform_fft(E2, *pSrc, kspace_buf_);
                    Gadgetron::crop(RO, encodingE1, E2, kspace_buf_, *pDst);
                    this->perform_ifft(E2, *pDst, *pDst);

                    pTmp = pSrc; pSrc = pDst; pDst = pTmp;
                }

                // adjust E2
                if (encodingE2 >= E2 + 1)
                {
                    Gadgetron::zero_pad_resize(*pSrc, RO, pSrc->get_size(1), encodingE2, *pDst);
                    pTmp = pSrc; pSrc = pDst; pDst = pTmp;
                }
                else if (encodingE2 <= E2 - 1)
                {
                    this->perform_fft(E2, *pSrc, kspace_buf_);
                    Gadgetron::crop(RO, pSrc->get_size(1), encodingE2, kspace_buf_, *pDst);
                    this->perform_ifft(E2, *pDst, *pDst);

                    pTmp = pSrc; pSrc = pDst; pDst = pTmp;
                }

                //adjust RO
                if (RO < reconSizeRO)
                {
                    Gadgetron::zero_pad_resize(*pSrc, reconSizeRO, pSrc->get_size(1), pSrc->get_size(2), *pDst);
                    pTmp = pSrc; pSrc = pDst; pDst = pTmp;
                }
                else if (RO > reconSizeRO)
                {
                    this->perform_fft(E2, *pSrc, kspace_buf_);
                    Gadgetron::crop(reconSizeRO, pSrc->get_size(1), pSrc->get_size(2), kspace_buf_, *pDst);
                    this->perform_ifft(E2, *pDst, *pDst);

                    pTmp = pSrc; pSrc = pDst; pDst = pTmp;
                }

                // final cut on image
                Gadgetron::crop(reconSizeRO, reconSizeE1, reconSizeE2, *pSrc, *pDst);

                if (pDst != &recon_res.data)
                {
                    recon_res.data = *pDst;
                }
            }
        }
//        catch (...)
//        {
//            GERROR_STREAM("Errors in GenericReconFieldOfViewAdjustmentGadget::adjust_FOV(ImageArray& data) ... ");
//            return GADGET_FAIL;
//        }

        return GADGET_OK;
    }


    int GenericReconFieldOfViewAdjustmentGadget::close(unsigned long flags)
    {
        GDEBUG_CONDITION_STREAM(this->verbose.value(), "GenericReconFieldOfViewAdjustmentGadget - close(flags) : " << flags);
        if (BaseClass::close(flags) != GADGET_OK) return GADGET_FAIL;
        this->gt_streamer_.close_stream_buffer();
        return GADGET_OK;
    }

    // ----------------------------------------------------------------------------------------

    GADGET_FACTORY_DECLARE(GenericReconFieldOfViewAdjustmentGadget)

}
