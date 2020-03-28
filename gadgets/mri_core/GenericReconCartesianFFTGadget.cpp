
#include "GenericReconCartesianFFTGadget.h"
#include "hoNDArray_reductions.h"

/*
    The input is IsmrmrdReconData and output is single 2D or 3D ISMRMRD images


    The image number computation logic is implemented in compute_image_number function, which can be overloaded
*/

namespace Gadgetron {

    GenericReconCartesianFFTGadget::GenericReconCartesianFFTGadget() : BaseClass()
    {
    }

    GenericReconCartesianFFTGadget::~GenericReconCartesianFFTGadget()
    {
    }

    int GenericReconCartesianFFTGadget::process_config(ACE_Message_Block* mb)
    {
        GADGET_CHECK_RETURN(BaseClass::process_config(mb) == GADGET_OK, GADGET_FAIL);

        // -------------------------------------------------

        ISMRMRD::IsmrmrdHeader h;
        try
        {
            deserialize(mb->rd_ptr(), h);
        }
        catch (...)
        {
            GDEBUG("Error parsing ISMRMRD Header");
        }

        size_t NE = h.encoding.size();
        num_encoding_spaces_ = NE;
        GDEBUG_CONDITION_STREAM(verbose.value(), "Number of encoding spaces: " << NE);

        recon_obj_.resize(NE);

        return GADGET_OK;
    }

    int GenericReconCartesianFFTGadget::process(Gadgetron::GadgetContainerMessage< IsmrmrdReconData >* m1)
    {
        if (perform_timing.value()) { gt_timer_local_.start("GenericReconCartesianFFTGadget::process"); }

        process_called_times_++;

        IsmrmrdReconData* recon_bit_ = m1->getObjectPtr();
        if (recon_bit_->rbit_.size() > num_encoding_spaces_)
        {
            GWARN_STREAM("Incoming recon_bit has more encoding spaces than the protocol : " << recon_bit_->rbit_.size() << " instead of " << num_encoding_spaces_);
        }

        GadgetContainerMessage< std::vector<ISMRMRD::Waveform> > * wav = AsContainerMessage< std::vector<ISMRMRD::Waveform>  >(m1->cont());
        if (wav)
        {
            if(verbose.value())
            {
                GDEBUG_STREAM("Incoming recon_bit with " << wav->getObjectPtr()->size() << " wave form samples ");
            }
        }

        // for every encoding space
        for (size_t e = 0; e < recon_bit_->rbit_.size(); e++)
        {
            std::stringstream os;
            os << "_encoding_" << e;

            GDEBUG_CONDITION_STREAM(verbose.value(), "Calling " << process_called_times_ << " , encoding space : " << e);
            GDEBUG_CONDITION_STREAM(verbose.value(), "======================================================================");

            // ---------------------------------------------------------------
            // export incoming data

            if (!debug_folder_full_path_.empty())
            {
                gt_exporter_.export_array_complex(recon_bit_->rbit_[e].data_.data_, debug_folder_full_path_ + "data" + os.str());
            }

            if (!debug_folder_full_path_.empty() && recon_bit_->rbit_[e].data_.trajectory_)
            {
                if (recon_bit_->rbit_[e].ref_->trajectory_->get_number_of_elements() > 0)
                {
                    gt_exporter_.export_array(*(recon_bit_->rbit_[e].data_.trajectory_), debug_folder_full_path_ + "data_traj" + os.str());
                }
            }

            // ---------------------------------------------------------------

            if (recon_bit_->rbit_[e].ref_)
            {
                if (!debug_folder_full_path_.empty())
                {
                    gt_exporter_.export_array_complex(recon_bit_->rbit_[e].ref_->data_, debug_folder_full_path_ + "ref" + os.str());
                }

                if (!debug_folder_full_path_.empty() && recon_bit_->rbit_[e].ref_->trajectory_)
                {
                    if (recon_bit_->rbit_[e].ref_->trajectory_->get_number_of_elements() > 0)
                    {
                        gt_exporter_.export_array(*(recon_bit_->rbit_[e].ref_->trajectory_), debug_folder_full_path_ + "ref_traj" + os.str());
                    }
                }

                // ---------------------------------------------------------------

                // after this step, the recon_obj_[e].ref_calib_ and recon_obj_[e].ref_coil_map_ are set

                if (perform_timing.value()) { gt_timer_.start("GenericReconCartesianFFTGadget::make_ref_coil_map"); }
                this->make_ref_coil_map(*recon_bit_->rbit_[e].ref_,*recon_bit_->rbit_[e].data_.data_.get_dimensions(), recon_obj_[e].ref_calib_, recon_obj_[e].ref_coil_map_, e);
                if (perform_timing.value()) { gt_timer_.stop(); }

                // ----------------------------------------------------------
                // export prepared ref for calibration and coil map
                if (!debug_folder_full_path_.empty())
                {
                    this->gt_exporter_.export_array_complex(recon_obj_[e].ref_calib_, debug_folder_full_path_ + "ref_calib" + os.str());
                }

                if (!debug_folder_full_path_.empty())
                {
                    this->gt_exporter_.export_array_complex(recon_obj_[e].ref_coil_map_, debug_folder_full_path_ + "ref_coil_map" + os.str());
                }


                // ---------------------------------------------------------------

                recon_bit_->rbit_[e].ref_ = boost::none;
            }

            if (recon_bit_->rbit_[e].data_.data_.get_number_of_elements() > 0)
            {
                if (!debug_folder_full_path_.empty())
                {
                    gt_exporter_.export_array_complex(recon_bit_->rbit_[e].data_.data_, debug_folder_full_path_ + "data_before_unwrapping" + os.str());
                }

                if (!debug_folder_full_path_.empty() && recon_bit_->rbit_[e].data_.trajectory_)
                {
                    if (recon_bit_->rbit_[e].data_.trajectory_->get_number_of_elements() > 0)
                    {
                        gt_exporter_.export_array(*(recon_bit_->rbit_[e].data_.trajectory_), debug_folder_full_path_ + "data_before_unwrapping_traj" + os.str());
                    }
                }

                // ---------------------------------------------------------------

                if (perform_timing.value()) { gt_timer_.start("GenericReconCartesianFFTGadget::perform_fft_combine"); }
                this->perform_fft_combine(recon_bit_->rbit_[e], recon_obj_[e], e);
                if (perform_timing.value()) { gt_timer_.stop(); }

                // ---------------------------------------------------------------

                if (perform_timing.value()) { gt_timer_.start("GenericReconCartesianFFTGadget::compute_image_header"); }
                this->compute_image_header(recon_bit_->rbit_[e], recon_obj_[e].recon_res_, e);
                if (perform_timing.value()) { gt_timer_.stop(); }

                // ---------------------------------------------------------------

                if (!debug_folder_full_path_.empty())
                {
                    this->gt_exporter_.export_array_complex(recon_obj_[e].recon_res_.data_, debug_folder_full_path_ + "recon_res" + os.str());
                }

                if (wav)
                {
                    recon_obj_[e].recon_res_.waveform_ = *(wav->getObjectPtr());
                }

                recon_obj_[e].recon_res_.acq_headers_ = recon_bit_->rbit_[e].data_.headers_;

                if (perform_timing.value()) { gt_timer_.start("GenericReconCartesianFFTGadget::send_out_image_array"); }
                this->send_out_image_array(recon_obj_[e].recon_res_, e, image_series.value() + ((int)e + 1), GADGETRON_IMAGE_REGULAR);
                if (perform_timing.value()) { gt_timer_.stop(); }
            }

            recon_obj_[e].recon_res_.data_.clear();
            recon_obj_[e].recon_res_.headers_.clear();
            recon_obj_[e].recon_res_.meta_.clear();
        }

        m1->release();

        if (perform_timing.value()) { gt_timer_local_.stop(); }

        return GADGET_OK;
    }

   


    void GenericReconCartesianFFTGadget::perform_fft_combine(IsmrmrdReconBit& recon_bit, ReconObjType& recon_obj, size_t e)
    {
        try
        {
            typedef std::complex<float> T;

            size_t RO = recon_bit.data_.data_.get_size(0);
            size_t E1 = recon_bit.data_.data_.get_size(1);
            size_t E2 = recon_bit.data_.data_.get_size(2);
            size_t dstCHA = recon_bit.data_.data_.get_size(3);
            size_t N = recon_bit.data_.data_.get_size(4);
            size_t S = recon_bit.data_.data_.get_size(5);
            size_t SLC = recon_bit.data_.data_.get_size(6);

            hoNDArray< std::complex<float> >& src = recon_obj.ref_calib_;

            size_t ref_RO = src.get_size(0);
            size_t ref_E1 = src.get_size(1);
            size_t ref_E2 = src.get_size(2);
            size_t srcCHA = src.get_size(3);
            size_t ref_N = src.get_size(4);
            size_t ref_S = src.get_size(5);
            size_t ref_SLC = src.get_size(6);

            recon_obj.recon_res_.data_.create(RO, E1, E2, 1, N, S, SLC);

            if (!debug_folder_full_path_.empty())
            {
                std::stringstream os;
                os << "encoding_" << e;
                std::string suffix = os.str();
                gt_exporter_.export_array_complex(recon_bit.data_.data_, debug_folder_full_path_ + "data_src_" + suffix);
            }

            // compute aliased images
            data_recon_buf_.create(RO, E1, E2, dstCHA, N, S, SLC);
	    

            if (E2>1)
            {
                Gadgetron::hoNDFFT<float>::instance()->ifft3c(recon_bit.data_.data_, complex_im_recon_buf_, data_recon_buf_);
            }
            else
            {
                Gadgetron::hoNDFFT<float>::instance()->ifft2c(recon_bit.data_.data_, complex_im_recon_buf_, data_recon_buf_);
            }

          

            if (!debug_folder_full_path_.empty())
            {
                std::stringstream os;
                os << "encoding_" << e;
                std::string suffix = os.str();
                gt_exporter_.export_array_complex(complex_im_recon_buf_, debug_folder_full_path_ + "aliasedIm_" + suffix);
            }

	    // combine channels
	    recon_obj.recon_res_.data_.create(RO, E1, E2, 1, N, S, SLC);
		
	    Gadgetron::sum_over_dimension(complex_im_recon_buf_, recon_obj.recon_res_.data_, 3);

        }
        catch (...)
        {
            GADGET_THROW("Errors happened in GenericReconCartesianFFTGadget::perform_fft_combine(...) ... ");
        }
    }

   

    GADGET_FACTORY_DECLARE(GenericReconCartesianFFTGadget)
}
