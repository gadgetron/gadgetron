
#include "PartialFourierAdjustROGadget.h"
#include "GadgetIsmrmrdReadWrite.h"

namespace Gadgetron
{

PartialFourierAdjustROGadget::PartialFourierAdjustROGadget() : maxRO_(0)
{

}

int PartialFourierAdjustROGadget::process_config(ACE_Message_Block* mb)
{
    boost::shared_ptr<ISMRMRD::ismrmrdHeader> cfg = parseIsmrmrdXMLHeader(std::string(mb->rd_ptr()));

    ISMRMRD::ismrmrdHeader::encoding_sequence e_seq = cfg->encoding();
    if (e_seq.size() != 1)
    {
        GADGET_DEBUG2("Number of encoding spaces: %d\n", e_seq.size());
        GADGET_DEBUG1("This simple partial fourier gadget only supports one encoding space\n");
        return GADGET_FAIL;
    }

    ISMRMRD::encodingSpaceType e_space = (*e_seq.begin()).encodedSpace();
    ISMRMRD::encodingSpaceType r_space = (*e_seq.begin()).reconSpace();
    ISMRMRD::encodingLimitsType e_limits = (*e_seq.begin()).encodingLimits();

    maxRO_ = e_space.matrixSize().x();

    GADGET_MSG("max RO : " << maxRO_);

    return GADGET_OK;
}

int addPrePostZeros(unsigned long long centre_column, unsigned long long samples)
{
    // 1 : pre zeros
    // 2 : post zeros
    // 0 : no zeros
    if ( 2*centre_column == samples )
    {
        return 0;
    }

    if ( 2*centre_column < samples )
    {
        return 1;
    }

    if ( 2*centre_column > samples )
    {
        return 2;
    }

    return 0;
}

int PartialFourierAdjustROGadget
::process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1,
        GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2)
{

    bool is_noise = ISMRMRD::FlagBit(ISMRMRD::ACQ_IS_NOISE_MEASUREMENT).isSet(m1->getObjectPtr()->flags);
    unsigned long long channels = m1->getObjectPtr()->active_channels;
    unsigned long long samples = m1->getObjectPtr()->number_of_samples;
    unsigned long long centre_column = m1->getObjectPtr()->center_sample;

    if (!is_noise) 
    {
        // adjust the center echo
        int az = addPrePostZeros(centre_column, samples);

        if ( az!= 0 && samples < maxRO_ )
        {
            GadgetContainerMessage< hoNDArray< std::complex<float> > >* m3 = new GadgetContainerMessage< hoNDArray< std::complex<float> > >();
            if (!m3)
            {
                return GADGET_FAIL;
            }

            std::vector<unsigned long long> data_out_dims = *m2->getObjectPtr()->get_dimensions();
            data_out_dims[0] = maxRO_;
            try
            {
                m3->getObjectPtr()->create(&data_out_dims);
            }
            catch(...)
            {
                GADGET_DEBUG1("Unable to create new data array for downsampled data\n");
                return GADGET_FAIL;
            }
            m3->getObjectPtr()->fill(0);

            std::complex<float>* pM3 = m3->getObjectPtr()->get_data_ptr();
            std::complex<float>* pM2 = m2->getObjectPtr()->get_data_ptr();

            unsigned long long c;
            if ( az == 1 ) // pre zeros
            {
                for ( c=0; c<channels; c++ )
                {
                    memcpy(pM3+c*maxRO_+maxRO_-samples, pM2+c*samples, sizeof( std::complex<float> )*samples);
                }
            }

            if ( az == 2 ) // post zeros
            {
                for ( c=0; c<channels; c++ )
                {
                    memcpy(pM3+c*maxRO_, pM2+c*samples, sizeof( std::complex<float> )*samples);
                }
            }

            m2->release(); //We are done with this data

            m1->cont(m3);
            m1->getObjectPtr()->number_of_samples = data_out_dims[0];
        }

        if (this->next()->putq(m1) == -1) 
        {
            ACE_ERROR_RETURN( (LM_ERROR,
                    ACE_TEXT("%p\n"),
                    ACE_TEXT("NoiseAdjustGadget::process, passing data on to next gadget")),
                    -1);
        }
    }
    else
    {
        if (this->next()->putq(m1) == -1) 
        {
            ACE_ERROR_RETURN( (LM_ERROR,
                    ACE_TEXT("%p\n"),
                    ACE_TEXT("NoiseAdjustGadget::process, passing data on to next gadget")),
                    -1);
        }
    }

    return GADGET_OK;
}

GADGET_FACTORY_DECLARE(PartialFourierAdjustROGadget)

}
