
#include "GenericReconCartesianFFTGadget.h"
#include "hoNDArray_reductions.h"

/*
    The input is IsmrmrdReconData and output is single 2D or 3D ISMRMRD images


    The image number computation logic is implemented in compute_image_number function, which can be overloaded
*/

namespace Gadgetron {


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


        // for every encoding space
        for (size_t e = 0; e < recon_bit_->rbit_.size(); e++)
        {
            std::stringstream os;
            os << "_encoding_" << e;

            GDEBUG_CONDITION_STREAM(verbose.value(), "Calling " << process_called_times_ << " , encoding space : " << e);
            GDEBUG_CONDITION_STREAM(verbose.value(), "======================================================================");

            ReconObjType  recon_obj;
            if (recon_bit_->rbit_[e].ref_)
            {
                this->make_ref_coil_map(*recon_bit_->rbit_[e].ref_, recon_bit_->rbit_[e].data_.data_.get_dimensions(), recon_obj.ref_calib_, recon_obj.ref_coil_map_, e);
                this->perform_coil_map_estimation(recon_obj.ref_coil_map_, recon_obj.coil_map_,e);
            } else {
                this->perform_coil_map_estimation(recon_bit_->rbit_[e].data_.data_,recon_obj.coil_map_,e);
            }

            if (recon_bit_->rbit_[e].data_.data_.get_number_of_elements() > 0)
            {
                if (perform_timing.value()) { gt_timer_.start("GenericReconCartesianFFTGadget::perform_fft_combine"); }
                this->perform_fft_combine(recon_bit_->rbit_[e], recon_obj, e);
                if (perform_timing.value()) { gt_timer_.stop(); }

                // ---------------------------------------------------------------

                if (perform_timing.value()) { gt_timer_.start("GenericReconCartesianFFTGadget::compute_image_header"); }
                this->compute_image_header(recon_bit_->rbit_[e], recon_obj.recon_res_, e);
                if (perform_timing.value()) { gt_timer_.stop(); }

                // ---------------------------------------------------------------
                recon_obj.recon_res_.acq_headers_ = recon_bit_->rbit_[e].data_.headers_;

                if (perform_timing.value()) { gt_timer_.start("GenericReconCartesianFFTGadget::send_out_image_array"); }
                this->send_out_image_array(recon_obj.recon_res_, e, image_series.value() + ((int)e + 1), GADGETRON_IMAGE_REGULAR);
                if (perform_timing.value()) { gt_timer_.stop(); }
            }
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


            recon_obj.recon_res_.data_.create(RO, E1, E2, 1, N, S, SLC);


            // compute aliased images
            data_recon_buf_.create(RO, E1, E2, dstCHA, N, S, SLC);
	    

            if (E2>1)
            {
                Gadgetron::hoNDFFT<float>::instance()->ifft3c(recon_bit.data_.data_, complex_im_recon_buf_, data_recon_buf_ );
            }
            else
            {
                Gadgetron::hoNDFFT<float>::instance()->ifft2c(recon_bit.data_.data_, complex_im_recon_buf_, data_recon_buf_);
            }

            coil_combine(complex_im_recon_buf_,recon_obj.coil_map_,3,recon_obj.recon_res_.data_);

        }
        catch (...)
        {
            GADGET_THROW("Errors happened in GenericReconCartesianFFTGadget::perform_fft_combine(...) ... ");
        }
    }

   

    GADGET_FACTORY_DECLARE(GenericReconCartesianFFTGadget)
}
