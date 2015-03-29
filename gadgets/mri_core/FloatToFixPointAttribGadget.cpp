/*
*       FloatToFixPointAttribGadget.cpp
*
*       Created on: March 10, 2014
*       Author: Hui Xue
*/

#include "GadgetIsmrmrdReadWrite.h"
#include "FloatToFixPointAttribGadget.h"
#include "mri_core_def.h"

namespace Gadgetron
{
    template <typename T> 
    FloatToFixPointAttribGadget<T>::FloatToFixPointAttribGadget() 
        : max_intensity_value_(std::numeric_limits<T>::max()), 
          min_intensity_value_(std::numeric_limits<T>::min()), 
          intensity_offset_value_(0)
    {
    }

    template <typename T> 
    FloatToFixPointAttribGadget<T>::~FloatToFixPointAttribGadget()
    {
    }

    template <typename T> 
    int FloatToFixPointAttribGadget<T>::process_config(ACE_Message_Block* mb)
    {
        // gadget parameters
        max_intensity_value_ = max_intensity.value();
        min_intensity_value_ = min_intensity.value();
        intensity_offset_value_ = intensity_offset.value();

        return GADGET_OK;
    }

    template <typename T> 
    int FloatToFixPointAttribGadget<T>::process(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1, GadgetContainerMessage< hoNDArray< float > >* m2, GadgetContainerMessage<ISMRMRD::MetaContainer>* m3)
    {
        GadgetContainerMessage<hoNDArray< T > > *cm2 =
            new GadgetContainerMessage<hoNDArray< T > >();

        boost::shared_ptr< std::vector<size_t> > dims = m2->getObjectPtr()->get_dimensions();

        try {cm2->getObjectPtr()->create(dims);}
        catch (std::runtime_error &err){
            GEXCEPTION(err,"Unable to create unsigned short storage in Extract Magnitude Gadget");
            return GADGET_FAIL;
        }

        float* src = m2->getObjectPtr()->get_data_ptr();
        T* dst = cm2->getObjectPtr()->get_data_ptr();

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
                    if (pix_val < (float)min_intensity_value_) pix_val = (float)min_intensity_value_;
                    if (pix_val > (float)max_intensity_value_) pix_val = (float)max_intensity_value_;
                    dst[i] = static_cast<T>(pix_val+0.5);
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
                    if (pix_val < (float)min_intensity_value_) pix_val = (float)min_intensity_value_;
                    if (pix_val > (float)max_intensity_value_) pix_val = (float)max_intensity_value_;
                    dst[i] = static_cast<T>(pix_val+0.5);
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
                    if (pix_val < (float)min_intensity_value_) pix_val = (float)min_intensity_value_;
                    if (pix_val > (float)max_intensity_value_) pix_val = (float)max_intensity_value_;
                    dst[i] = static_cast<T>(pix_val);
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

        if (typeid(T) == typeid(unsigned short))
        {
            m1->getObjectPtr()->data_type = ISMRMRD::ISMRMRD_USHORT;
        }
        else if (typeid(T) == typeid(short))
        {
            m1->getObjectPtr()->data_type = ISMRMRD::ISMRMRD_SHORT;
        }
        else if (typeid(T) == typeid(unsigned int))
        {
            m1->getObjectPtr()->data_type = ISMRMRD::ISMRMRD_UINT;
        }
        else if (typeid(T) == typeid(int))
        {
            m1->getObjectPtr()->data_type = ISMRMRD::ISMRMRD_INT;
        }
        else
        {
            GDEBUG("Unknown data type, bailing out\n");
            m1->release();
            cm2->release();
            return GADGET_FAIL;
        }

        if (this->next()->putq(m1) == -1)
        {
            m1->release();
            GDEBUG("Unable to put unsigned short magnitude image on next gadgets queue");
            return GADGET_FAIL;
        }

        return GADGET_OK;
    }

    FloatToUShortAttribGadget::FloatToUShortAttribGadget()
    {
        max_intensity.value(4095);
        min_intensity.value(0);
        intensity_offset.value(2048);

        max_intensity_value_ = 4095;
        min_intensity_value_ = 0;
        intensity_offset_value_ = 2048;
    }

    FloatToUShortAttribGadget::~FloatToUShortAttribGadget()
    {
    }

    FloatToShortAttribGadget::FloatToShortAttribGadget()
    {
    }

    FloatToShortAttribGadget::~FloatToShortAttribGadget()
    {
    }

    FloatToUIntAttribGadget::FloatToUIntAttribGadget()
    {
    }

    FloatToUIntAttribGadget::~FloatToUIntAttribGadget()
    {
    }

    FloatToIntAttribGadget::FloatToIntAttribGadget()
    {
    }

    FloatToIntAttribGadget::~FloatToIntAttribGadget()
    {
    }

    GADGET_FACTORY_DECLARE(FloatToUShortAttribGadget)
    GADGET_FACTORY_DECLARE(FloatToShortAttribGadget)
    GADGET_FACTORY_DECLARE(FloatToIntAttribGadget)
    GADGET_FACTORY_DECLARE(FloatToUIntAttribGadget)
}
