
#include "GenericReconKSpaceFilteringGadget.h"
#include <iomanip>

#include "mri_core_kspace_filter.h"
#include "hoNDArray_reductions.h"
#include "mri_core_def.h"
#include "hoNDArray_utils.h"
#include "hoNDArray_elemwise.h"
#include "hoNDFFT.h"
#include "mri_core_utility.h"

namespace Gadgetron {

    GenericReconKSpaceFilteringGadget::GenericReconKSpaceFilteringGadget() : BaseClass()
    {
    }

    GenericReconKSpaceFilteringGadget::~GenericReconKSpaceFilteringGadget()
    {
    }

    int GenericReconKSpaceFilteringGadget::process_config(const mrd::Header& header)
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

            GDEBUG_CONDITION_STREAM(verbose.value(), "Encoding space : " << e << " - encoding FOV : [" << encoding_FOV_[e][0] << " " << encoding_FOV_[e][1] << " " << encoding_FOV_[e][2] << " ]");
            GDEBUG_CONDITION_STREAM(verbose.value(), "Encoding space : " << e << " - recon    FOV : [" << recon_FOV_[e][0]    << " " << recon_FOV_[e][1]    << " " << recon_FOV_[e][2] << " ]");
        }

        filter_RO_.resize(NE);
        filter_E1_.resize(NE);
        filter_E2_.resize(NE);

        return GADGET_OK;
    }

    int GenericReconKSpaceFilteringGadget::process(Gadgetron::GadgetContainerMessage< mrd::ImageArray >* m1)
    {
        if (perform_timing.value()) { gt_timer_.start("GenericReconKSpaceFilteringGadget::process"); }

        GDEBUG_CONDITION_STREAM(verbose.value(), "GenericReconKSpaceFilteringGadget::process(...) starts ... ");

        process_called_times_++;

        mrd::ImageArray* recon_res_ = m1->getObjectPtr();

        // print out recon info
        if (verbose.value())
        {
            GDEBUG_STREAM("----> GenericReconKSpaceFilteringGadget::process(...) has been called " << process_called_times_ << " times ...");
            std::stringstream os;
            recon_res_->data.print(os);
            GDEBUG_STREAM(os.str());
        }

        // some images do not need kspace filter
        if (recon_res_->meta[0].count(skip_processing_meta_field) && recon_res_->meta[0][skip_processing_meta_field].size() > 0)
        {
            GDEBUG_CONDITION_STREAM(verbose.value(), "Skip kspace filtering for this image array ... ");

            if (this->next()->putq(m1) == -1)
            {
                GERROR("GenericReconKSpaceFilteringGadget::process, passing incoming image array on to next gadget");
                return GADGET_FAIL;
            }

            return GADGET_OK;
        }

        size_t encoding = (size_t)std::get<long>(recon_res_->meta[0]["encoding"].front());
        GADGET_CHECK_RETURN(encoding<num_encoding_spaces_, GADGET_FAIL);

        std::string dataRole = std::get<std::string>(recon_res_->meta[0][GADGETRON_DATA_ROLE].front());

        std::stringstream os;
        os << "encoding_" << encoding << "_" << dataRole;
        std::string str = os.str();

        size_t RO = recon_res_->data.get_size(0);
        size_t E1 = recon_res_->data.get_size(1);
        size_t E2 = recon_res_->data.get_size(2);
        size_t CHA = recon_res_->data.get_size(3);
        size_t N = recon_res_->data.get_size(4);
        size_t S = recon_res_->data.get_size(5);
        size_t SLC = recon_res_->data.get_size(6);

        // perform SNR unit scaling
        mrd::SamplingLimits sampling_limits;

        if (recon_res_->meta[0].count("sampling_limits_RO")) {
            auto& sl = recon_res_->meta[0]["sampling_limits_RO"];
            sampling_limits.kspace_encoding_step_0.minimum = (uint32_t)std::get<long>(sl[0]);
            sampling_limits.kspace_encoding_step_0.center = (uint32_t)std::get<long>(sl[1]);
            sampling_limits.kspace_encoding_step_0.maximum = (uint32_t)std::get<long>(sl[2]);
        }

        if ( !( (sampling_limits.kspace_encoding_step_0.minimum >= 0) && (sampling_limits.kspace_encoding_step_0.maximum < RO) && (sampling_limits.kspace_encoding_step_0.minimum <= sampling_limits.kspace_encoding_step_0.maximum)) )
        {
            sampling_limits.kspace_encoding_step_0.minimum    = 0;
            sampling_limits.kspace_encoding_step_0.center = RO / 2;
            sampling_limits.kspace_encoding_step_0.maximum    = RO - 1;
        }

        if (recon_res_->meta[0].count("sampling_limits_E1")) {
            auto& sl = recon_res_->meta[0]["sampling_limits_E1"];
            sampling_limits.kspace_encoding_step_1.minimum = (uint32_t)std::get<long>(sl[0]);
            sampling_limits.kspace_encoding_step_1.center = (uint32_t)std::get<long>(sl[1]);
            sampling_limits.kspace_encoding_step_1.maximum = (uint32_t)std::get<long>(sl[2]);
        }

        if (!((sampling_limits.kspace_encoding_step_1.minimum >= 0) && (sampling_limits.kspace_encoding_step_1.maximum < E1) && (sampling_limits.kspace_encoding_step_1.minimum <= sampling_limits.kspace_encoding_step_1.maximum)))
        {
            sampling_limits.kspace_encoding_step_1.minimum    = 0;
            sampling_limits.kspace_encoding_step_1.center = E1 / 2;
            sampling_limits.kspace_encoding_step_1.maximum    = E1 - 1;
        }

        if (recon_res_->meta[0].count("sampling_limits_E2")) {
            auto& sl = recon_res_->meta[0]["sampling_limits_E2"];
            sampling_limits.kspace_encoding_step_2.minimum = (uint32_t)std::get<long>(sl[0]);
            sampling_limits.kspace_encoding_step_2.center = (uint32_t)std::get<long>(sl[1]);
            sampling_limits.kspace_encoding_step_2.maximum = (uint32_t)std::get<long>(sl[2]);
        }

        if (!((sampling_limits.kspace_encoding_step_2.minimum >= 0) && (sampling_limits.kspace_encoding_step_2.maximum < E2) && (sampling_limits.kspace_encoding_step_2.minimum <= sampling_limits.kspace_encoding_step_2.maximum)))
        {
            sampling_limits.kspace_encoding_step_2.minimum    = 0;
            sampling_limits.kspace_encoding_step_2.center = E2 / 2;
            sampling_limits.kspace_encoding_step_2.maximum    = E2 - 1;
        }

        // ----------------------------------------------------------
        // filter create if needed
        // ----------------------------------------------------------
        size_t ii;

        if (filter_RO_[encoding].get_number_of_elements() != RO)
        {
            if (sampling_limits.kspace_encoding_step_0.minimum == 0 || sampling_limits.kspace_encoding_step_0.maximum == RO - 1)
            {
                if (filterRO.value() != "None")
                {
                    filter_RO_[encoding].create(RO);
                    Gadgetron::generate_symmetric_filter(RO, filter_RO_[encoding], Gadgetron::get_kspace_filter_type(filterRO.value()), filterRO_sigma.value(), (size_t)std::ceil(filterRO_width.value()*RO));
                }
            }
            else
            {
                if (filterRO.value() != "None")
                {
                    size_t len;
                    this->find_kspace_sampled_range(sampling_limits.kspace_encoding_step_0.minimum, sampling_limits.kspace_encoding_step_0.maximum, RO, len);

                    hoNDArray< std::complex<float> > f;
                    Gadgetron::generate_symmetric_filter(len, f, Gadgetron::get_kspace_filter_type(filterRO.value()), filterRO_sigma.value(), (size_t)std::ceil(filterRO_width.value()*len));
                    Gadgetron::pad(RO, f, filter_RO_[encoding]);
                }
            }

            if (filter_RO_[encoding].get_number_of_elements() == RO)
            {
                if (sampling_limits.kspace_encoding_step_0.minimum != 0 || sampling_limits.kspace_encoding_step_0.maximum != RO - 1)
                {
                    // compensate the sacling from min_ to max_
                    std::complex<float> sos = 0.0f;
                    for (ii = sampling_limits.kspace_encoding_step_0.minimum; ii <= sampling_limits.kspace_encoding_step_0.maximum; ii++)
                    {
                        sos += filter_RO_[encoding](ii)* std::conj(filter_RO_[encoding](ii));
                    }

                    Gadgetron::scal((float)(1.0 / std::sqrt(std::abs(sos) / (sampling_limits.kspace_encoding_step_0.maximum - sampling_limits.kspace_encoding_step_0.minimum + 1))), filter_RO_[encoding]);
                }
            }

            if (!debug_folder_full_path_.empty())
            {
                if(filter_RO_[encoding].get_number_of_elements()>0)
                    gt_exporter_.export_array_complex(filter_RO_[encoding], debug_folder_full_path_ + "filterRO_" + str);
            }
        }

        if (filter_E1_[encoding].get_number_of_elements() != E1)
        {
            if (sampling_limits.kspace_encoding_step_1.minimum == 0 || sampling_limits.kspace_encoding_step_1.maximum == E1 - 1)
            {
                if (filterE1.value() != "None")
                {
                    filter_E1_[encoding].create(E1);
                    Gadgetron::generate_symmetric_filter(E1, filter_E1_[encoding], Gadgetron::get_kspace_filter_type(filterE1.value()), filterE1_sigma.value(), (size_t)std::ceil(filterE1_width.value()*E1));
                }
            }
            else
            {
                if (filterE1.value() != "None")
                {
                    size_t len;
                    this->find_kspace_sampled_range(sampling_limits.kspace_encoding_step_1.minimum, sampling_limits.kspace_encoding_step_1.maximum, E1, len);

                    hoNDArray< std::complex<float> > f;
                    Gadgetron::generate_symmetric_filter(len, f, Gadgetron::get_kspace_filter_type(filterE1.value()), filterE1_sigma.value(), (size_t)std::ceil(filterE1_width.value()*len));
                    Gadgetron::pad(E1, f, filter_E1_[encoding]);
                }
            }

            if (filter_E1_[encoding].get_number_of_elements() == E1)
            {
                if (sampling_limits.kspace_encoding_step_1.minimum != 0 || sampling_limits.kspace_encoding_step_1.maximum != E1 - 1)
                {
                    // compensate the sacling from min_ to max_
                    std::complex<float> sos = 0.0f;
                    for (ii = sampling_limits.kspace_encoding_step_1.minimum; ii <= sampling_limits.kspace_encoding_step_1.maximum; ii++)
                    {
                        sos += filter_E1_[encoding](ii)* std::conj(filter_E1_[encoding](ii));
                    }

                    Gadgetron::scal((float)(1.0 / std::sqrt(std::abs(sos) / (sampling_limits.kspace_encoding_step_1.maximum - sampling_limits.kspace_encoding_step_1.minimum + 1))), filter_E1_[encoding]);
                }
            }

            if (!debug_folder_full_path_.empty())
            {
                if (filter_E1_[encoding].get_number_of_elements()>0)
                    gt_exporter_.export_array_complex(filter_E1_[encoding], debug_folder_full_path_ + "filterE1_" + str);
            }
        }

        if (E2>1 && filter_E2_[encoding].get_number_of_elements() != E2)
        {
            if (sampling_limits.kspace_encoding_step_2.minimum == 0 || sampling_limits.kspace_encoding_step_2.maximum == E1 - 1)
            {
                if (filterE2.value() != "None")
                {
                    filter_E2_[encoding].create(E2);
                    Gadgetron::generate_symmetric_filter(E2, filter_E2_[encoding], Gadgetron::get_kspace_filter_type(filterE2.value()), filterE2_sigma.value(), (size_t)std::ceil(filterE2_width.value()*E2));
                }
            }
            else
            {
                if (filterE2.value() != "None")
                {
                    size_t len;
                    this->find_kspace_sampled_range(sampling_limits.kspace_encoding_step_2.minimum, sampling_limits.kspace_encoding_step_2.maximum, E2, len);

                    hoNDArray< std::complex<float> > f;
                    Gadgetron::generate_symmetric_filter(len, f, Gadgetron::get_kspace_filter_type(filterE2.value()), filterE2_sigma.value(), (size_t)std::ceil(filterE2_width.value()*len));
                    Gadgetron::pad(E2, f, filter_E2_[encoding]);
                }
            }

            if (filter_E2_[encoding].get_number_of_elements() == E2)
            {
                if (sampling_limits.kspace_encoding_step_2.minimum != 0 || sampling_limits.kspace_encoding_step_2.maximum != E1 - 1)
                {
                    // compensate the sacling from min_ to max_
                    std::complex<float> sos = 0.0f;
                    for (ii = sampling_limits.kspace_encoding_step_2.minimum; ii <= sampling_limits.kspace_encoding_step_2.maximum; ii++)
                    {
                        sos += filter_E2_[encoding](ii)* std::conj(filter_E2_[encoding](ii));
                    }

                    Gadgetron::scal((float)(1.0 / std::sqrt(std::abs(sos) / (sampling_limits.kspace_encoding_step_2.maximum - sampling_limits.kspace_encoding_step_2.minimum + 1))), filter_E2_[encoding]);
                }
            }

            if (!debug_folder_full_path_.empty())
            {
                if (filter_E2_[encoding].get_number_of_elements()>0)
                    gt_exporter_.export_array_complex(filter_E2_[encoding], debug_folder_full_path_ + "filterE2_" + str);
            }
        }

        // ----------------------------------------------------------
        // perform kspace filter
        // ----------------------------------------------------------
        if ( (filter_RO_[encoding].get_number_of_elements() == RO)
                || (filter_E1_[encoding].get_number_of_elements() == E1)
                || ((E2>1) && (filter_E2_[encoding].get_number_of_elements() == E2)) )
        {
            // ----------------------------------------------------------
            // go to kspace
            // ----------------------------------------------------------
            if (!debug_folder_full_path_.empty()) { gt_exporter_.export_array_complex(recon_res_->data, debug_folder_full_path_ + "image_before_filtering_" + str); }

            if (perform_timing.value()) { gt_timer_.start("GenericReconKSpaceFilteringGadget: fftc"); }
            if (E2 > 1)
            {
                Gadgetron::hoNDFFT<float>::instance()->fft3c(recon_res_->data, kspace_buf_);
            }
            else
            {
                Gadgetron::hoNDFFT<float>::instance()->fft2c(recon_res_->data, kspace_buf_);
            }
            if (perform_timing.value()) { gt_timer_.stop(); }

            if (!debug_folder_full_path_.empty()) { gt_exporter_.export_array_complex(kspace_buf_, debug_folder_full_path_ + "kspace_before_filtering_" + str); }

            // ----------------------------------------------------------
            // filtering
            // ----------------------------------------------------------

            bool inKSpace = true;

            if ( (filter_RO_[encoding].get_number_of_elements() == RO)
                && (filter_E1_[encoding].get_number_of_elements() == E1)
                && (E2>1) && (filter_E2_[encoding].get_number_of_elements() == E2) )
            {
                if (perform_timing.value()) { gt_timer_.start("GenericReconKSpaceFilteringGadget, apply_kspace_filter_ROE1E2 ... "); }
                Gadgetron::apply_kspace_filter_ROE1E2(kspace_buf_, filter_RO_[encoding], filter_E1_[encoding], filter_E2_[encoding], filter_res_);
                if (perform_timing.value()) { gt_timer_.stop(); }

                if (!debug_folder_full_path_.empty()) { gt_exporter_.export_array_complex(filter_res_, debug_folder_full_path_ + "kspace_after_filtered_" + str); }
                inKSpace = true;
            }
            else if ( (filter_RO_[encoding].get_number_of_elements() == RO) && (filter_E1_[encoding].get_number_of_elements() == E1) )
            {
                if (perform_timing.value()) { gt_timer_.start("GenericReconKSpaceFilteringGadget, apply_kspace_filter_ROE1 ... "); }
                Gadgetron::apply_kspace_filter_ROE1(kspace_buf_, filter_RO_[encoding], filter_E1_[encoding], filter_res_);
                if (perform_timing.value()) { gt_timer_.stop(); }

                if (!debug_folder_full_path_.empty()) { gt_exporter_.export_array_complex(filter_res_, debug_folder_full_path_ + "kspace_after_filtered_" + str); }

                inKSpace = true;
            }
            else
            {
                filter_res_ = kspace_buf_;

                hoNDArray< std::complex<float> >* pSrc = &kspace_buf_;
                hoNDArray< std::complex<float> >* pDst = &filter_res_;

                bool filterPerformed = false;

                if (filter_RO_[encoding].get_number_of_elements() == RO)
                {
                    Gadgetron::apply_kspace_filter_RO(*pSrc, filter_RO_[encoding], *pDst);
                    std::swap(pSrc, pDst);
                    filterPerformed = true;
                }

                if (filter_E1_[encoding].get_number_of_elements() == E1)
                {
                    Gadgetron::apply_kspace_filter_E1(*pSrc, filter_E1_[encoding], *pDst);
                    std::swap(pSrc, pDst);
                    filterPerformed = true;
                }

                if (filter_E2_[encoding].get_number_of_elements() == E2)
                {
                    Gadgetron::apply_kspace_filter_E2(*pSrc, filter_E2_[encoding], *pDst);
                    std::swap(pSrc, pDst);
                    filterPerformed = true;
                }

                if (filterPerformed)
                {
                    if (pDst != &filter_res_)
                    {
                        filter_res_ = *pDst;
                    }
                }

                inKSpace = true;
            }

            if (!debug_folder_full_path_.empty()) { gt_exporter_.export_array_complex(filter_res_, debug_folder_full_path_ + "kspace_after_filtering_" + str); }

            // ----------------------------------------------------------
            // go back to image domain
            // ----------------------------------------------------------
            if (E2 > 1)
            {
                Gadgetron::hoNDFFT<float>::instance()->ifft3c(filter_res_, recon_res_->data);
            }
            else
            {
                Gadgetron::hoNDFFT<float>::instance()->ifft2c(filter_res_, recon_res_->data);
            }

            if (!debug_folder_full_path_.empty()) { gt_exporter_.export_array_complex(recon_res_->data, debug_folder_full_path_ + "image_after_filtering_" + str); }

            GDEBUG_CONDITION_STREAM(verbose.value(), "GenericReconKSpaceFilteringGadget::process(...) ends ... ");
        }

        // ----------------------------------------------------------
        // send out results
        // ----------------------------------------------------------
        if (this->next()->putq(m1) == -1)
        {
            GERROR("GenericReconKSpaceFilteringGadget::process, passing data on to next gadget");
            return GADGET_FAIL;
        }

        if (perform_timing.value()) { gt_timer_.stop(); }

        return GADGET_OK;
    }

    void GenericReconKSpaceFilteringGadget::find_kspace_sampled_range(size_t min, size_t max, size_t len, size_t& r)
    {
        GADGET_CHECK_THROW(min < max);
        GADGET_CHECK_THROW(min < len);
        GADGET_CHECK_THROW(max < len);

        size_t hmin = len / 2 - min;
        size_t hmax = max - len / 2;
        r = (hmax>hmin) ? 2*hmax+1 : 2*hmin+1;

        if (r > len) r = len;
    }


    // ----------------------------------------------------------------------------------------

    GADGET_FACTORY_DECLARE(GenericReconKSpaceFilteringGadget)

}
