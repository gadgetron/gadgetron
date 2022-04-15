
#include "GenericReconNoiseStdMapComputingGadget.h"
#include <iomanip>

#include "mri_core_kspace_filter.h"
#include "hoNDArray_reductions.h"
#include "mri_core_def.h"
#include "hoNDArray_utils.h"
#include "mri_core_utility.h"
#include "hoNDArray_elemwise.h"

namespace Gadgetron {

    GenericReconNoiseStdMapComputingGadget::GenericReconNoiseStdMapComputingGadget() : BaseClass()
    {
    }

    GenericReconNoiseStdMapComputingGadget::~GenericReconNoiseStdMapComputingGadget()
    {
    }

    int GenericReconNoiseStdMapComputingGadget::process_config(ACE_Message_Block* mb)
    {
        GADGET_CHECK_RETURN(BaseClass::process_config(mb) == GADGET_OK, GADGET_FAIL);

        ISMRMRD::IsmrmrdHeader h;
        try
        {
            deserialize(mb->rd_ptr(), h);
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

        return GADGET_OK;
    }

    int GenericReconNoiseStdMapComputingGadget::process(Gadgetron::GadgetContainerMessage< IsmrmrdImageArray >* m1)
    {
        GDEBUG_CONDITION_STREAM(verbose.value(), "GenericReconNoiseStdMapComputingGadget::process(...) starts ... ");

        process_called_times_++;

        IsmrmrdImageArray* recon_res_ = m1->getObjectPtr();

        // print out recon info
        if (verbose.value())
        {
            GDEBUG_STREAM("----> GenericReconNoiseStdMapComputingGadget::process(...) has been called " << process_called_times_ << " times ...");
            std::stringstream os;
            recon_res_->data_.print(os);
            GDEBUG_STREAM(os.str());
        }

        std::string dataRole = std::string(recon_res_->meta_[0].as_str(GADGETRON_DATA_ROLE));
        if (dataRole != GADGETRON_IMAGE_SNR_MAP)
        {
            if (this->next()->putq(m1) == -1)
            {
                GERROR("GenericReconNoiseStdMapComputingGadget::process, passing images on to next gadget");
                return GADGET_FAIL;
            }

            return GADGET_OK;
        }

        size_t encoding = (size_t)recon_res_->meta_[0].as_long("encoding", 0);
        GADGET_CHECK_RETURN(encoding<num_encoding_spaces_, GADGET_FAIL);

        std::stringstream os;
        os << "encoding_" << encoding;
        std::string str = os.str();

        size_t RO = recon_res_->data_.get_size(0);
        size_t E1 = recon_res_->data_.get_size(1);
        size_t E2 = recon_res_->data_.get_size(2);
        size_t CHA = recon_res_->data_.get_size(3);
        size_t N = recon_res_->data_.get_size(4);
        size_t S = recon_res_->data_.get_size(5);
        size_t SLC = recon_res_->data_.get_size(6);

        // perform std map computation
        if (N < start_N_for_std_map.value())
        {
            GWARN_STREAM("GenericReconNoiseStdMapComputingGadget, N < start_N_for_std_map.value() - " << N << " - " << start_N_for_std_map.value());

            if (this->next()->putq(m1) == -1)
            {
                GERROR("GenericReconNoiseStdMapComputingGadget::process, passing images on to next gadget");
                return GADGET_FAIL;
            }

            return GADGET_OK;
        }

        // make a copy for results
        Gadgetron::GadgetContainerMessage<IsmrmrdImageArray>* cm1 = new Gadgetron::GadgetContainerMessage<IsmrmrdImageArray>();
        *(cm1->getObjectPtr()) = *(m1->getObjectPtr());

        // pass on the incoming image array
        if (this->next()->putq(m1) == -1)
        {
            GERROR("GenericReconNoiseStdMapComputingGadget::process, passing images on to next gadget");
            return GADGET_FAIL;
        }

        if (!debug_folder_full_path_.empty())
        {
            gt_exporter_.export_array_complex(cm1->getObjectPtr()->data_, debug_folder_full_path_ + "incoming_SNR_images_" + str);
        }

        // compute std map
        size_t startN = start_N_for_std_map.value();

        hoNDArray<T> repBuf(RO, E1, E2, CHA, N - startN);
        hoNDArray<real_value_type> repBufMag(RO, E1, E2, CHA, N - startN);
        hoNDArray<real_value_type> repBufMagMean(RO, E1, E2, CHA, 1);
        hoNDArray<real_value_type> diffMag(RO, E1, E2, CHA);
        hoNDArray<real_value_type> stdMap(RO, E1, E2, CHA, 1, S, SLC);
        Gadgetron::clear(stdMap);

        hoNDArray<T>& snrMap = cm1->getObjectPtr()->data_;

        size_t n, s, slc;

        for (slc = 0; slc < SLC; slc++)
        {
            for (s = 0; s < S; s++)
            {
                for (n = startN; n < N; n++)
                {
                    T* pSNRMap = &(snrMap(0, 0, 0, 0, n, s, slc));
                    memcpy(repBuf.begin() + (n - startN)*RO*E1*E2*CHA, pSNRMap, sizeof(T)*RO*E1*E2*CHA);
                }

                Gadgetron::abs(repBuf, repBufMag);

                // compute mean
                Gadgetron::sum_over_dimension(repBufMag, repBufMagMean, 4);
                Gadgetron::scal((real_value_type)(1.0 / repBufMag.get_size(4)), repBufMagMean);

                hoNDArray<real_value_type> stdMapCurr;
                stdMapCurr.create(RO, E1, E2, CHA, &(stdMap(0, 0, 0, 0, 0, s, slc)));

                for (n = startN; n < N; n++)
                {
                    hoNDArray<real_value_type> snr;
                    snr.create(RO, E1, E2, CHA, &(repBufMag(0, 0, 0, 0, n - startN)));

                    Gadgetron::subtract(snr, repBufMagMean, diffMag);
                    Gadgetron::multiply(diffMag, diffMag, diffMag);

                    Gadgetron::add(stdMapCurr, diffMag, stdMapCurr);
                }

                Gadgetron::scal((real_value_type)(1.0 / (N - startN)), stdMapCurr);
                Gadgetron::sqrt(stdMapCurr, stdMapCurr);
            }
        }

        if (!debug_folder_full_path_.empty())
        {
            gt_exporter_.export_array(stdMap, debug_folder_full_path_ + "std_map_" + str);
        }

        GDEBUG_CONDITION_STREAM(verbose.value(), "GenericReconNoiseStdMapComputingGadget::process(...) ends ... ");

        // update image headers
        cm1->getObjectPtr()->data_.clear();
        Gadgetron::real_to_complex(stdMap, cm1->getObjectPtr()->data_);

        size_t num = cm1->getObjectPtr()->meta_.size();
        for (size_t n = 0; n < num; n++)
        {
            cm1->getObjectPtr()->meta_[n].set(GADGETRON_DATA_ROLE, GADGETRON_IMAGE_STD_MAP);

            cm1->getObjectPtr()->meta_[n].set(GADGETRON_IMAGECOMMENT, "GT");
            cm1->getObjectPtr()->meta_[n].append(GADGETRON_IMAGECOMMENT, "SNR");
            cm1->getObjectPtr()->meta_[n].append(GADGETRON_IMAGECOMMENT, GADGETRON_IMAGE_STD_MAP);

            cm1->getObjectPtr()->meta_[n].set(GADGETRON_SEQUENCEDESCRIPTION, "_STD_MAP");
            cm1->getObjectPtr()->meta_[n].set(GADGETRON_USE_DEDICATED_SCALING_FACTOR, (long)1);
        }

        num = cm1->getObjectPtr()->headers_.get_number_of_elements();
        for (size_t n = 0; n < num; n++)
        {
            cm1->getObjectPtr()->headers_(n).image_series_index *= 10;
        }

        // ----------------------------------------------------------
        // send out results
        // ----------------------------------------------------------
        if (this->next()->putq(cm1) == -1)
        {
            GERROR("GenericReconNoiseStdMapComputingGadget::process, passing data on to next gadget");
            return GADGET_FAIL;
        }

        return GADGET_OK;
    }



    // ----------------------------------------------------------------------------------------

    GADGET_FACTORY_DECLARE(GenericReconNoiseStdMapComputingGadget)

}
