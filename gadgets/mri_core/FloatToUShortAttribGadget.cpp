/*
*       FloatToUShortAttribGadget.cpp
*
*       Created on: March 10, 2014
*       Author: Hui Xue
*/

#include "GadgetIsmrmrdReadWrite.h"
#include "FloatToUShortAttribGadget.h"
#include "mri_core_def.h"

namespace Gadgetron
{
    FloatToUShortAttribGadget::FloatToUShortAttribGadget() : max_intensity_value_(4095), intensity_offset_value_(2048)
    {
    }

    FloatToUShortAttribGadget::~FloatToUShortAttribGadget()
    {
    }

    int FloatToUShortAttribGadget::process_config(ACE_Message_Block* mb)
    {
        // gadget parameters
         max_intensity_value_ = max_intensity.value();
        if ( max_intensity_value_ == 0 ) max_intensity_value_ = 4095;

        intensity_offset_value_ = intensity_offset.value();

        return GADGET_OK;
    }

    int FloatToUShortAttribGadget::process(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1, GadgetContainerMessage< hoNDArray< float > >* m2, GadgetContainerMessage<ISMRMRD::MetaContainer>* m3)
    {
        GadgetContainerMessage<hoNDArray< ACE_UINT16 > > *cm2 =
            new GadgetContainerMessage<hoNDArray< ACE_UINT16 > >();

        boost::shared_ptr< std::vector<size_t> > dims = m2->getObjectPtr()->get_dimensions();

        try {cm2->getObjectPtr()->create(dims);}
        catch (std::runtime_error &err){
            GEXCEPTION(err,"Unable to create unsigned short storage in Extract Magnitude Gadget");
            return GADGET_FAIL;
        }

        float* src = m2->getObjectPtr()->get_data_ptr();
        ACE_UINT16* dst = cm2->getObjectPtr()->get_data_ptr();

        long long i;
        long long numOfPixels = (long long)cm2->getObjectPtr()->get_number_of_elements();

        switch (m1->getObjectPtr()->image_type)
        {
            case ISMRMRD::ISMRMRD_IMTYPE_MAGNITUDE:
            {
                #pragma omp parallel for default(none) private(i) shared(numOfPixels, src, dst)
                for (i=0; i<numOfPixels; i++)
                {
                    float pix_val = src[i];
                    pix_val = std::abs(pix_val);
                    if (pix_val > max_intensity_value_) pix_val = max_intensity_value_;
                    dst[i] = static_cast<unsigned short>(pix_val+0.5);
                }
            }
            break;

            case ISMRMRD::ISMRMRD_IMTYPE_REAL:
            case ISMRMRD::ISMRMRD_IMTYPE_IMAG:
            {
                #pragma omp parallel for default(none) private(i) shared(numOfPixels, src, dst)
                for (i=0; i<numOfPixels; i++)
                {
                    float pix_val = src[i];
                    pix_val = pix_val + intensity_offset_value_;
                    if (pix_val < 0) pix_val = 0;
                    if (pix_val > max_intensity_value_) pix_val = max_intensity_value_;
                    dst[i] = static_cast<unsigned short>(pix_val+0.5);
                }

                if ( m3->getObjectPtr()->length(GADGETRON_IMAGE_WINDOWCENTER) > 0 )
                {
                    long windowCenter;
                    windowCenter = m3->getObjectPtr()->as_long(GADGETRON_IMAGE_WINDOWCENTER, 0);
                    m3->getObjectPtr()->set(GADGETRON_IMAGE_WINDOWCENTER, windowCenter+(long)intensity_offset_value_);
                }
            }
            break;

            case ISMRMRD::ISMRMRD_IMTYPE_PHASE:
            {
                #pragma omp parallel for default(none) private(i) shared(numOfPixels, src, dst)
                for (i=0; i<numOfPixels; i++)
                {
                    float pix_val = src[i];
                    pix_val *= (float)(intensity_offset_value_/3.14159265);
                    pix_val += intensity_offset_value_;
                    if (pix_val < 0) pix_val = 0;
                    if (pix_val > max_intensity_value_) pix_val = max_intensity_value_;
                    dst[i] = static_cast<unsigned short>(pix_val);
                }
            }
            break;

            default:
                GDEBUG("Unknown image type %d, bailing out\n",m1->getObjectPtr()->image_type);
                m1->release();
                cm2->release();
                return GADGET_FAIL;
        }

        m1->cont(cm2);
        cm2->cont(m3);

        m2->cont(NULL);
        m2->release();

        m1->getObjectPtr()->data_type = ISMRMRD::ISMRMRD_USHORT;

        if (this->next()->putq(m1) == -1)
        {
            m1->release();
            GDEBUG("Unable to put unsigned short magnitude image on next gadgets queue");
            return GADGET_FAIL;
        }

        return GADGET_OK;
    }

    GADGET_FACTORY_DECLARE(FloatToUShortAttribGadget)
}
