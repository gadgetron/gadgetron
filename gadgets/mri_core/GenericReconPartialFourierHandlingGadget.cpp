
#include "GenericReconPartialFourierHandlingGadget.h"
#include <iomanip>

#include "mri_core_kspace_filter.h"
#include "hoNDArray_reductions.h"
#include "mri_core_def.h"

namespace Gadgetron { 

    GenericReconPartialFourierHandlingGadget::GenericReconPartialFourierHandlingGadget() : BaseClass()
    {
        startRO_ = 0;
        endRO_ = 0;

        startE1_ = 0;
        endE1_ = 0;

        startE2_ = 0;
        endE2_ = 0;
    }

    GenericReconPartialFourierHandlingGadget::~GenericReconPartialFourierHandlingGadget()
    {
    }

    int GenericReconPartialFourierHandlingGadget::process_config(ACE_Message_Block* mb)
    {
        GADGET_CHECK_RETURN(BaseClass::process_config(mb) == GADGET_OK, GADGET_FAIL);

        ISMRMRD::IsmrmrdHeader h;
        try
        {
            deserialize(mb->rd_ptr(),h);
        }
        catch (...)
        {
            GDEBUG("Error parsing ISMRMRD Header");
        }

        if (!h.acquisitionSystemInformation)
        {
            GDEBUG("acquisitionSystemInformation not found in header. Bailing out");
            return GADGET_FAIL;
        }

        // -------------------------------------------------

        size_t NE = h.encoding.size();

        num_encoding_spaces_ = NE;

        GDEBUG_CONDITION_STREAM(verbose.value(), "Number of encoding spaces: " << NE);

        acceFactorE1_.resize(NE, 1);
        acceFactorE2_.resize(NE, 1);

        size_t e;
        for (e = 0; e < h.encoding.size(); e++)
        {
            if (!h.encoding[e].parallelImaging)
            {
                GDEBUG_STREAM("Parallel Imaging section not found in header for encoding " << e);
                acceFactorE1_[e] = 1;
                acceFactorE2_[e] = 1;
            }
            else
            {
                ISMRMRD::ParallelImaging p_imaging = *h.encoding[0].parallelImaging;

                acceFactorE1_[e] = p_imaging.accelerationFactor.kspace_encoding_step_1;
                acceFactorE2_[e] = p_imaging.accelerationFactor.kspace_encoding_step_2;
                GDEBUG_CONDITION_STREAM(verbose.value(), "acceFactorE1 is " << acceFactorE1_[e]);
                GDEBUG_CONDITION_STREAM(verbose.value(), "acceFactorE2 is " << acceFactorE2_[e]);
            }
        }

        return GADGET_OK;
    }

    int GenericReconPartialFourierHandlingGadget::process(Gadgetron::GadgetContainerMessage< IsmrmrdImageArray >* m1)
    {
        if (perform_timing.value()) { gt_timer_local_.start("GenericReconPartialFourierHandlingGadget::process"); }

        GDEBUG_CONDITION_STREAM(verbose.value(), "GenericReconPartialFourierHandlingGadget::process(...) starts ... ");

        process_called_times_++;

        IsmrmrdImageArray* recon_res_ = m1->getObjectPtr();

        // print out recon info
        if (verbose.value())
        {
            GDEBUG_STREAM("----> GenericReconPartialFourierHandlingGadget::process(...) has been called " << process_called_times_ << " times ...");
            std::stringstream os;
            recon_res_->data_.print(os);
            GDEBUG_STREAM(os.str());
        }

        // some images do not need partial fourier handling processing
        if (recon_res_->meta_[0].length(skip_processing_meta_field.value().c_str())>0)
        {
            if (this->next()->putq(m1) == -1)
            {
                GERROR("GenericReconPartialFourierHandlingGadget::process, passing incoming image array on to next gadget");
                return GADGET_FAIL;
            }

            if (perform_timing.value()) { gt_timer_local_.stop(); }

            return GADGET_OK;
        }

        // call the partial foureir

        size_t encoding = (size_t)recon_res_->meta_[0].as_long("encoding", 0);
        GADGET_CHECK_RETURN(encoding<num_encoding_spaces_, GADGET_FAIL);

        std::string dataRole = std::string(recon_res_->meta_[0].as_str(GADGETRON_DATA_ROLE));

        std::stringstream os;
        os << "encoding_" << encoding << "_" << dataRole;
        std::string str = os.str();

        size_t RO = recon_res_->data_.get_size(0);
        size_t E1 = recon_res_->data_.get_size(1);
        size_t E2 = recon_res_->data_.get_size(2);
        size_t CHA = recon_res_->data_.get_size(3);
        size_t N = recon_res_->data_.get_size(4);
        size_t S = recon_res_->data_.get_size(5);
        size_t SLC = recon_res_->data_.get_size(6);

        // perform SNR unit scaling
        SamplingLimit sampling_limits[3];

        if (recon_res_->meta_[0].length("sampling_limits_RO")>0)
        {
            sampling_limits[0].min_     = (uint16_t)recon_res_->meta_[0].as_long("sampling_limits_RO", 0);
            sampling_limits[0].center_  = (uint16_t)recon_res_->meta_[0].as_long("sampling_limits_RO", 1);
            sampling_limits[0].max_     = (uint16_t)recon_res_->meta_[0].as_long("sampling_limits_RO", 2);
        }

        if (!((sampling_limits[0].min_ >= 0) && (sampling_limits[0].max_ < RO) && (sampling_limits[0].min_ <= sampling_limits[0].max_)))
        {
            sampling_limits[0].min_     = 0;
            sampling_limits[0].center_  = RO / 2;
            sampling_limits[0].max_     = RO - 1;
        }

        if (recon_res_->meta_[0].length("sampling_limits_E1") > 0)
        {
            sampling_limits[1].min_     = (uint16_t)recon_res_->meta_[0].as_long("sampling_limits_E1", 0);
            sampling_limits[1].center_  = (uint16_t)recon_res_->meta_[0].as_long("sampling_limits_E1", 1);
            sampling_limits[1].max_     = (uint16_t)recon_res_->meta_[0].as_long("sampling_limits_E1", 2);
        }

        if (!((sampling_limits[1].min_ >= 0) && (sampling_limits[1].max_ < E1) && (sampling_limits[1].min_ <= sampling_limits[1].max_)))
        {
            sampling_limits[1].min_     = 0;
            sampling_limits[1].center_  = E1 / 2;
            sampling_limits[1].max_     = E1 - 1;
        }

        if (recon_res_->meta_[0].length("sampling_limits_E2") > 0)
        {
            sampling_limits[2].min_     = (uint16_t)recon_res_->meta_[0].as_long("sampling_limits_E2", 0);
            sampling_limits[2].center_  = (uint16_t)recon_res_->meta_[0].as_long("sampling_limits_E2", 1);
            sampling_limits[2].max_     = (uint16_t)recon_res_->meta_[0].as_long("sampling_limits_E2", 2);
        }

        if (!((sampling_limits[2].min_ >= 0) && (sampling_limits[2].max_ < E2) && (sampling_limits[2].min_ <= sampling_limits[2].max_)))
        {
            sampling_limits[2].min_     = 0;
            sampling_limits[2].center_  = E2 / 2;
            sampling_limits[2].max_     = E2 - 1;
        }

        // ----------------------------------------------------------
        // pf kspace sampling range
        // ----------------------------------------------------------
        // if image padding is performed, those dimension may not need partial fourier handling

        startRO_ = sampling_limits[0].min_;
        endRO_ = sampling_limits[0].max_;

        startE1_ = 0;
        endE1_ = E1 - 1;

        startE2_ = 0;
        endE2_ = E2 - 1;

        if (std::abs((double)(sampling_limits[1].max_ - E1 / 2) - (double)(E1 / 2 - sampling_limits[1].min_)) > acceFactorE1_[encoding])
        {
            startE1_ = sampling_limits[1].min_;
            endE1_ = sampling_limits[1].max_;
        }

        if ((E2>1) && (std::abs((double)(sampling_limits[2].max_ - E2 / 2) - (double)(E2 / 2 - sampling_limits[2].min_)) > acceFactorE2_[encoding]))
        {
            startE2_ = sampling_limits[2].min_;
            endE2_ = sampling_limits[2].max_;
        }

        long lenRO = endRO_ - startRO_ + 1;
        long lenE1 = endE1_ - startE1_ + 1;
        long lenE2 = endE2_ - startE2_ + 1;

        if (lenRO == RO && lenE1 == E1 && lenE2 == E2)
        {
            GDEBUG_CONDITION_STREAM(verbose.value(), "lenRO == RO && lenE1 == E1 && lenE2 == E2");

            if (this->next()->putq(m1) == -1)
            {
                GERROR("GenericReconPartialFourierHandlingGadget::process, passing data on to next gadget");
                return GADGET_FAIL;
            }

            if (perform_timing.value()) { gt_timer_local_.stop(); }

            return GADGET_OK;
        }

        // ----------------------------------------------------------
        // go to kspace
        // ----------------------------------------------------------
        if (E2 > 1)
        {
            Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft3c(recon_res_->data_, kspace_buf_);
        }
        else
        {
            Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft2c(recon_res_->data_, kspace_buf_);
        }

        if (!debug_folder_full_path_.empty())
        {
            gt_exporter_.export_array_complex(kspace_buf_, debug_folder_full_path_ + "kspace_before_pf");
        }

        // ----------------------------------------------------------
        // pf handling
        // ----------------------------------------------------------
        GADGET_CHECK_RETURN(this->perform_partial_fourier_handling() == GADGET_OK, GADGET_FAIL);

        if (!debug_folder_full_path_.empty())
        {
            gt_exporter_.export_array_complex(pf_res_, debug_folder_full_path_ + "kspace_after_pf");
        }

        // ----------------------------------------------------------
        // go back to image domain
        // ----------------------------------------------------------
        if (E2 > 1)
        {
            Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft3c(pf_res_, recon_res_->data_);
        }
        else
        {
            Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft2c(pf_res_, recon_res_->data_);
        }

        if (!debug_folder_full_path_.empty())
        {
            gt_exporter_.export_array_complex(recon_res_->data_, debug_folder_full_path_ + "data_after_pf");
        }

        GDEBUG_CONDITION_STREAM(verbose.value(), "GenericReconPartialFourierHandlingGadget::process(...) ends ... ");

        // ----------------------------------------------------------
        // send out results
        // ----------------------------------------------------------
        if (this->next()->putq(m1) == -1)
        {
            GERROR("GenericReconPartialFourierHandlingGadget::process, passing data on to next gadget");
            return GADGET_FAIL;
        }

        if (perform_timing.value()) { gt_timer_local_.stop(); }

        return GADGET_OK;
    }

    int GenericReconPartialFourierHandlingGadget::close(unsigned long flags)
    {
        GDEBUG_CONDITION_STREAM(true, "GenericReconPartialFourierHandlingGadget - close(flags) : " << flags);

        if ( BaseClass::close(flags) != GADGET_OK ) return GADGET_FAIL;

        if ( flags != 0 )
        {
        }

        return GADGET_OK;
    }

    // ----------------------------------------------------------------------------------------
}
