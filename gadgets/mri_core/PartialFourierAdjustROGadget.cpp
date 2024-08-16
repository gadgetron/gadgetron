#include "PartialFourierAdjustROGadget.h"
#include "ismrmrd/xml.h"

namespace Gadgetron
{

PartialFourierAdjustROGadget::PartialFourierAdjustROGadget() : maxRO_(0)
{

}

int PartialFourierAdjustROGadget::process_config(ACE_Message_Block* mb)
{
  ISMRMRD::IsmrmrdHeader h;
  deserialize(mb->rd_ptr(),h);

  if (h.encoding.size() != 1) {
    GDEBUG("Number of encoding spaces: %d\n", h.encoding.size());
    GDEBUG("This partial fourier gadget only supports one encoding space\n");
    return GADGET_FAIL;
  }

  ISMRMRD::EncodingSpaceType e_space = h.encoding[0].encodedSpace;
  maxRO_ = e_space.matrixSize.x;
  GDEBUG_STREAM("max RO : " << maxRO_);
  return GADGET_OK;
}

int addPrePostZeros(size_t centre_column, size_t samples)
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

    bool is_noise = ISMRMRD::FlagBit(ISMRMRD::ISMRMRD_ACQ_IS_NOISE_MEASUREMENT).isSet(m1->getObjectPtr()->flags);
    size_t channels = m1->getObjectPtr()->active_channels;
    size_t samples = m1->getObjectPtr()->number_of_samples;
    size_t centre_column = m1->getObjectPtr()->center_sample;

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

            std::vector<size_t> data_out_dims = m2->getObjectPtr()->get_dimensions();
            data_out_dims[0] = maxRO_;
            try
            {
                m3->getObjectPtr()->create(&data_out_dims);
            }
            catch(...)
            {
                GDEBUG("Unable to create new data array for downsampled data\n");
                return GADGET_FAIL;
            }
            m3->getObjectPtr()->fill(0);

            std::complex<float>* pM3 = m3->getObjectPtr()->get_data_ptr();
            std::complex<float>* pM2 = m2->getObjectPtr()->get_data_ptr();

            size_t c;
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
	  GERROR("NoiseAdjustGadget::process, passing data on to next gadget");
	  return -1;
        }
    }
    else
    {
        if (this->next()->putq(m1) == -1) 
        {
	  GERROR("NoiseAdjustGadget::process, passing data on to next gadget");
	  return -1;
        }
    }

    return GADGET_OK;
}

GADGET_FACTORY_DECLARE(PartialFourierAdjustROGadget)

}
