
#include "GenericReconPartialFourierHandlingGadget.h"
#include <iomanip>

#include "mri_core_kspace_filter.h"
#include "hoNDArray_reductions.h"
#include "mri_core_def.h"

namespace Gadgetron { 

    GenericReconPartialFourierHandlingGadget::GenericReconPartialFourierHandlingGadget()
    {
        num_encoding_spaces_ = 1;

        process_called_times_ = 0;
    }

    GenericReconPartialFourierHandlingGadget::~GenericReconPartialFourierHandlingGadget()
    {
    }

    int GenericReconPartialFourierHandlingGadget::process_config(ACE_Message_Block* mb)
    {
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
                GDEBUG("Parallel Imaging section not found in header");
                return GADGET_FAIL;
            }

            ISMRMRD::ParallelImaging p_imaging = *h.encoding[0].parallelImaging;

            acceFactorE1_[e] = p_imaging.accelerationFactor.kspace_encoding_step_1;
            acceFactorE2_[e] = p_imaging.accelerationFactor.kspace_encoding_step_2;
            GDEBUG_CONDITION_STREAM(verbose.value(), "acceFactorE1 is " << acceFactorE1_[e]);
            GDEBUG_CONDITION_STREAM(verbose.value(), "acceFactorE2 is " << acceFactorE2_[e]);
        }

        // ---------------------------------------------------------------------------------------------------------
        // generate the destination folder
        /*if (!debug_folder.value().empty())
        {
            Gadgetron::get_debug_folder_path(debug_folder.value(), debug_folder_full_path_);
            GDEBUG_CONDITION_STREAM(verbose.value(), "Debug folder is " << debug_folder_full_path_);
        }
        else
        {
            GDEBUG_CONDITION_STREAM(verbose.value(), "Debug folder is not set ... ");
        }*/

        return GADGET_OK;
    }

    int GenericReconPartialFourierHandlingGadget::process(Gadgetron::GadgetContainerMessage< IsmrmrdImageArray >* m1)
    {
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

        std::string dataRole = std::string(recon_res_->meta_[0].as_str(GADGETRON_DATA_ROLE));

        // some images do not need partial fourier handling processing
        if (dataRole == GADGETRON_IMAGE_GFACTOR)
        {
            if (this->next()->putq(m1) == -1)
            {
                GERROR("GenericReconPartialFourierHandlingGadget::process, passing gfactor on to next gadget");
                return GADGET_FAIL;
            }

            return GADGET_OK;
        }

        // call the partial foureir

        size_t encoding = (size_t)recon_res_->meta_[0].as_long("encoding", 0);
        GADGET_CHECK_RETURN(encoding<num_encoding_spaces_, GADGET_FAIL);

        // perform SNR unit scaling
        SamplingLimit sampling_limits[3];

        sampling_limits[0].min_    = (uint16_t)recon_res_->meta_[0].as_long("sampling_limits_RO", 0);
        sampling_limits[0].center_ = (uint16_t)recon_res_->meta_[0].as_long("sampling_limits_RO", 1);
        sampling_limits[0].max_    = (uint16_t)recon_res_->meta_[0].as_long("sampling_limits_RO", 2);

        sampling_limits[1].min_    = (uint16_t)recon_res_->meta_[0].as_long("sampling_limits_E1", 0);
        sampling_limits[1].center_ = (uint16_t)recon_res_->meta_[0].as_long("sampling_limits_E1", 1);
        sampling_limits[1].max_    = (uint16_t)recon_res_->meta_[0].as_long("sampling_limits_E1", 2);

        sampling_limits[2].min_    = (uint16_t)recon_res_->meta_[0].as_long("sampling_limits_E2", 0);
        sampling_limits[2].center_ = (uint16_t)recon_res_->meta_[0].as_long("sampling_limits_E2", 1);
        sampling_limits[2].max_    = (uint16_t)recon_res_->meta_[0].as_long("sampling_limits_E2", 2);

        size_t RO  = recon_res_->data_.get_size(0);
        size_t E1  = recon_res_->data_.get_size(1);
        size_t E2  = recon_res_->data_.get_size(2);
        size_t CHA = recon_res_->data_.get_size(3);
        size_t N   = recon_res_->data_.get_size(4);
        size_t S   = recon_res_->data_.get_size(5);
        size_t SLC = recon_res_->data_.get_size(6);

        // ----------------------------------------------------------
        // pf SNR scaling
        // ----------------------------------------------------------
        real_value_type partialFourierCompensationFactor = 1;

        long lenRO = sampling_limits[0].max_ - sampling_limits[0].min_ + 1;
        long lenE1 = sampling_limits[1].max_ - sampling_limits[1].min_ + 1;
        long lenE2 = sampling_limits[2].max_ - sampling_limits[2].min_ + 1;

        if (lenRO == RO && lenE1 == E1 && lenE2 == E2)
        {
            GDEBUG_CONDITION_STREAM(verbose.value(), "lenRO == RO && lenE1 == E1 && lenE2 == E2");

            if (this->next()->putq(m1) == -1)
            {
                GERROR("GenericReconPartialFourierHandlingGadget::process, passing data on to next gadget");
                return GADGET_FAIL;
            }

            return GADGET_OK;
        }

        if (lenRO < RO)
        {
            partialFourierCompensationFactor *= (real_value_type)(RO) / (real_value_type)(lenRO);
        }

        if (lenE1 < E1)
        {
            partialFourierCompensationFactor *= (real_value_type)(E1) / (real_value_type)(lenE1);
        }

        if (E2>1 && lenE2 < E2)
        {
            partialFourierCompensationFactor *= (real_value_type)(E2) / (real_value_type)(lenE2);
        }

        partialFourierCompensationFactor = std::sqrt(partialFourierCompensationFactor);
        if (verbose.value())
        {
            GDEBUG_STREAM("Partial fourier scaling factor : " << partialFourierCompensationFactor);
        }

        /*if (!debug_folder_full_path_.empty())
        {
            gt_exporter_.exportArrayComplex(recon_res_->data_, debug_folder_full_path_ + "data_before_pf");
        }*/

        if (partialFourierCompensationFactor>1)
        {
            Gadgetron::scal(partialFourierCompensationFactor, recon_res_->data_);
        }

        if (partial_fourier_algo.value() == "None")
        {
            GDEBUG_CONDITION_STREAM(verbose.value(), "partial_fourier_algo.value() == None");

            if (this->next()->putq(m1) == -1)
            {
                GERROR("GenericReconPartialFourierHandlingGadget::process, passing data on to next gadget");
                return GADGET_FAIL;
            }

            return GADGET_OK;
        }

        // not all PF handling methods support SNR images
        // only zero-filling filter PF handling supports SNR images
        if (dataRole == GADGETRON_IMAGE_SNR_MAP && partial_fourier_algo.value() != "ZeroFillingFilter")
        {
            if (this->next()->putq(m1) == -1)
            {
                GERROR("GenericReconPartialFourierHandlingGadget::process, passing SNR images on to next gadget");
                return GADGET_FAIL;
            }

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

        /*if (!debug_folder_full_path_.empty())
        {
            gt_exporter_.exportArrayComplex(kspace_buf_, debug_folder_full_path_ + "kspace_before_pf");
        }*/

        // ----------------------------------------------------------
        // pf handling
        // ----------------------------------------------------------

        // if image padding is performed, those dimension may not need partial fourier handling
        size_t startE1(0), endE1(E1 - 1);
        size_t startE2(0), endE2(E2 - 1);

        if (std::abs((double)(sampling_limits[1].max_ - E1 / 2) - (double)(E1 / 2 - sampling_limits[1].min_)) > acceFactorE1_[encoding])
        {
            startE1 = sampling_limits[1].min_;
            endE1 = sampling_limits[1].max_;
        }

        if ( (E2>1) && (std::abs((double)(sampling_limits[2].max_ - E2 / 2) - (double)(E2 / 2 - sampling_limits[2].min_)) > acceFactorE2_[encoding]) ) 
        {
            startE2 = sampling_limits[2].min_;
            endE2 = sampling_limits[2].max_;
        }

        if (partial_fourier_algo.value() == "ZeroFillingFilter")
        {
            if (perform_timing.value()) { gt_timer_.start("GenericReconPartialFourierHandlingGadget, partial_fourier_filter"); }
            GADGET_CHECK_EXCEPTION_RETURN(Gadgetron::partial_fourier_filter(kspace_buf_,
                sampling_limits[0].min_, sampling_limits[0].max_, startE1, endE1, startE2, endE2,
                partial_fourier_filter_RO_width.value(), partial_fourier_filter_E1_width.value(),
                partial_fourier_filter_E2_width.value(), partial_fourier_filter_densityComp.value(),
                filter_pf_RO_, filter_pf_E1_, filter_pf_E2_, pf_res_), GADGET_FAIL);
            if (perform_timing.value()) { gt_timer_.stop(); }
        }
        else if (partial_fourier_algo.value() == "POCS")
        {
            if (perform_timing.value()) { gt_timer_.start("GenericReconPartialFourierHandlingGadget, partial_fourier_POCS"); }
            GADGET_CHECK_EXCEPTION_RETURN(Gadgetron::partial_fourier_POCS(kspace_buf_,
                sampling_limits[0].min_, sampling_limits[0].max_, startE1, endE1, startE2, endE2,
                partial_fourier_POCS_transitBand.value(), partial_fourier_POCS_transitBand.value(),
                partial_fourier_POCS_transitBand_E2.value(), partial_fourier_POCS_iters.value(),
                partial_fourier_POCS_thres.value(), pf_res_), GADGET_FAIL);
            if (perform_timing.value()) { gt_timer_.stop(); }
        }
        else
        {
            GERROR_STREAM("Incorrect partial fourier handling algorithm type ... ");

            if (this->next()->putq(m1) == -1)
            {
                GERROR("GenericReconPartialFourierHandlingGadget::process, passing data on to next gadget");
                return GADGET_FAIL;
            }

            return GADGET_OK;
        }

        /*if (!debug_folder_full_path_.empty())
        {
            gt_exporter_.exportArrayComplex(pf_res_, debug_folder_full_path_ + "kspace_after_pf");
        }*/

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

        /*if (!debug_folder_full_path_.empty())
        {
            gt_exporter_.exportArrayComplex(recon_res_->data_, debug_folder_full_path_ + "data_after_pf");
        }*/

        GDEBUG_CONDITION_STREAM(verbose.value(), "GenericReconPartialFourierHandlingGadget::process(...) ends ... ");

        // ----------------------------------------------------------
        // send out results
        // ----------------------------------------------------------
        if (this->next()->putq(m1) == -1)
        {
            GERROR("GenericReconPartialFourierHandlingGadget::process, passing data on to next gadget");
            return GADGET_FAIL;
        }

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

    GADGET_FACTORY_DECLARE(GenericReconPartialFourierHandlingGadget)

}
