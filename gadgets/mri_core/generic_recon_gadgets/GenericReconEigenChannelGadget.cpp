
#include "GenericReconEigenChannelGadget.h"
#include <iomanip>

#include "hoNDArray_reductions.h"
#include "mri_core_def.h"

namespace Gadgetron {

    GenericReconEigenChannelGadget::GenericReconEigenChannelGadget() : BaseClass()
    {
    }

    GenericReconEigenChannelGadget::~GenericReconEigenChannelGadget()
    {
    }

    int GenericReconEigenChannelGadget::process_config(const mrd::Header& header)
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

        calib_mode_.resize(NE, mrd::CalibrationMode::kNoacceleration);

        KLT_.resize(NE);

        for (size_t e = 0; e < h.encoding.size(); e++)
        {
            mrd::EncodingSpaceType e_space = h.encoding[e].encoded_space;
            mrd::EncodingSpaceType r_space = h.encoding[e].recon_space;
            mrd::EncodingLimitsType e_limits = h.encoding[e].encoding_limits;

            calib_mode_[e] = mrd::CalibrationMode::kNoacceleration;

            if (!h.encoding[e].parallel_imaging)
            {
                GDEBUG_STREAM("Parallel Imaging section not found in header");
            }
            else
            {
                mrd::ParallelImagingType p_imaging = *h.encoding[0].parallel_imaging;

                if (p_imaging.acceleration_factor.kspace_encoding_step_1 > 1 || p_imaging.acceleration_factor.kspace_encoding_step_2 > 1)
                {
                    calib_mode_[e] = p_imaging.calibration_mode.value_or(mrd::CalibrationMode::kNoacceleration);
                }
            }
        }

        return GADGET_OK;
    }

    int GenericReconEigenChannelGadget::process(Gadgetron::GadgetContainerMessage< mrd::ReconData >* m1)
    {
        if (perform_timing.value()) { gt_timer_.start("GenericReconEigenChannelGadget::process"); }

        process_called_times_++;

        mrd::ReconData* recon_data = m1->getObjectPtr();
        if (recon_data->buffers.size() > num_encoding_spaces_)
        {
            GWARN_STREAM("Incoming recon_bit has more encoding spaces than the protocol : " << recon_data->buffers.size() << " instead of " << num_encoding_spaces_);
        }

        // for every encoding space, prepare the recon_data->buffers[e].ref
        size_t e, n, s, slc;
        for (e = 0; e < recon_data->buffers.size(); e++)
        {
            auto & rbit = recon_data->buffers[e];

            hoNDArray< std::complex<float> >& data = rbit.data.data;

            if (data.get_number_of_elements()==0) continue;

            size_t RO = data.get_size(0);
            size_t E1 = data.get_size(1);
            size_t E2 = data.get_size(2);
            size_t CHA = data.get_size(3);
            size_t N = data.get_size(4);
            size_t S = data.get_size(5);
            size_t SLC = data.get_size(6);

            GDEBUG_STREAM("GenericReconEigenChannelGadget - incoming data array : [RO E1 E2 CHA N S SLC] - [" << RO << " " << E1 << " " << E2 << " " << CHA << " " << N << " " << S << " " << SLC << "]");

            if(data.get_number_of_elements()==0)
            {
                m1->release();
                return GADGET_OK;
            }

            // whether it is needed to update coefficients
            bool recompute_coeff = false;
            if ( (KLT_[e].size()!=SLC) || update_eigen_channel_coefficients.value() )
            {
                recompute_coeff = true;
            }
            else
            {
                if(KLT_[e].size() == SLC)
                {
                    for (slc = 0; slc < SLC; slc++)
                    {
                        if (KLT_[e][slc].size() != S)
                        {
                            recompute_coeff = true;
                            break;
                        }
                        else
                        {
                            for (s = 0; s < S; s++)
                            {
                                if (KLT_[e][slc][s].size() != N)
                                {
                                    recompute_coeff = true;
                                    break;
                                }
                            }
                        }
                    }
                }
            }

                bool average_N = average_all_ref_N.value();
                bool average_S = average_all_ref_S.value();

            if(recompute_coeff)
            {
                if(rbit.ref)
                {
                    // use ref to compute coefficients
                    Gadgetron::compute_eigen_channel_coefficients(rbit.ref->data, average_N, average_S,
                        (calib_mode_[e] == mrd::CalibrationMode::kInterleaved), N, S, upstream_coil_compression_thres.value(), upstream_coil_compression_num_modesKept.value(), KLT_[e]);
                }
                else
                {
                    // use data to compute coefficients
                    Gadgetron::compute_eigen_channel_coefficients(data, average_N, average_S,
                        (calib_mode_[e] == mrd::CalibrationMode::kInterleaved), N, S, upstream_coil_compression_thres.value(), upstream_coil_compression_num_modesKept.value(), KLT_[e]);
                }

                if (verbose.value())
                {
                    hoNDArray< std::complex<float> > E;

                    for (slc = 0; slc < SLC; slc++)
                    {
                        for (s = 0; s < S; s++)
                        {
                            for (n = 0; n < N; n++)
                            {
                                KLT_[e][slc][s][n].eigen_value(E);
                            }
                        }
                    }
                }
                else
                {
                    if(average_N && average_S)
                    {
                        for (slc = 0; slc < SLC; slc++)
                        {
                            GDEBUG_STREAM("GenericReconEigenChannelGadget - Number of modes kept, SLC : " << slc << " - " << KLT_[e][slc][0][0].output_length() << " out of " << CHA);
                        }
                    }
                    else if(average_N && !average_S)
                    {
                        for (slc = 0; slc < SLC; slc++)
                        {
                            for (s = 0; s < S; s++)
                            {
                                GDEBUG_STREAM("GenericReconEigenChannelGadget - Number of modes kept, [SLC S] : [" << slc << " " << s << "] - " << KLT_[e][slc][s][0].output_length() << " out of " << CHA);
                            }
                        }
                    }
                    else if(!average_N && average_S)
                    {
                        for (slc = 0; slc < SLC; slc++)
                        {
                            for (n = 0; n < N; n++)
                            {
                                GDEBUG_STREAM("GenericReconEigenChannelGadget - Number of modes kept, [SLC N] : [" << slc << " " << n << "] - " << KLT_[e][slc][0][n].output_length() << " out of " << CHA);
                            }
                        }
                    }
                    else if(!average_N && !average_S)
                    {
                        for (slc = 0; slc < SLC; slc++)
                        {
                            for (s = 0; s < S; s++)
                            {
                                for (n = 0; n < N; n++)
                                {
                                    GDEBUG_STREAM("GenericReconEigenChannelGadget - Number of modes kept, [SLC S N] : [" << slc << " " << s << " " << n << "] - " << KLT_[e][slc][s][n].output_length() << " out of " << CHA);
                                }
                            }
                        }
                    }
                }
            }

            std::stringstream os;
            os << "_encoding_" << e;

            if (!debug_folder_full_path_.empty())
            {
                gt_exporter_.export_array_complex(data, debug_folder_full_path_ + "data_before_KLT" + os.str());
            }

            // apply KL coefficients
            Gadgetron::apply_eigen_channel_coefficients(KLT_[e], data);

            if (!debug_folder_full_path_.empty())
            {
                gt_exporter_.export_array_complex(data, debug_folder_full_path_ + "data_after_KLT" + os.str());
            }

            if (rbit.ref)
            {
                if (!debug_folder_full_path_.empty())
                {
                    gt_exporter_.export_array_complex(rbit.ref->data, debug_folder_full_path_ + "ref_before_KLT" + os.str());
                }

                Gadgetron::apply_eigen_channel_coefficients(KLT_[e], rbit.ref->data);

                if (!debug_folder_full_path_.empty())
                {
                    gt_exporter_.export_array_complex(rbit.ref->data, debug_folder_full_path_ + "ref_after_KLT" + os.str());
                }
            }
        }

        if (perform_timing.value()) { gt_timer_.stop(); }

        if (this->next()->putq(m1) < 0)
        {
            GERROR_STREAM("Put ReconData to Q failed ... ");
            return GADGET_FAIL;
        }

        return GADGET_OK;
    }

    GADGET_FACTORY_DECLARE(GenericReconEigenChannelGadget)
}
