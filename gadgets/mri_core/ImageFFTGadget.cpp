/*
*       ImageFFTGadget.cpp
*       Author: Hui Xue
*/

#include "ImageFFTGadget.h"
#include "hoNDArray_elemwise.h"
#include "hoNDImage_util.h"
#include "hoNDFFT.h"
#include "mri_core_def.h"

namespace Gadgetron
{
    ImageFFTGadget::ImageFFTGadget()
    {
    }

    ImageFFTGadget::~ImageFFTGadget()
    {
    }

    int ImageFFTGadget::process(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1, GadgetContainerMessage< hoNDArray< ValueType > >* m2, GadgetContainerMessage <ISMRMRD::MetaContainer>* m3)
    {
        ArrayType* input_array = m2->getObjectPtr();

        ArrayType output_array(*input_array);

        std::vector<size_t> dims;
        input_array->get_dimensions(dims);

        size_t RO = dims[0];
        size_t E1 = dims[1];
        size_t E2 = 1;

        if (E2 > 1)
        {
            Gadgetron::hoNDFFT<float>::instance()->fft3c(*input_array, output_array);
        }
        else
        {
            Gadgetron::hoNDFFT<float>::instance()->fft2c(*input_array, output_array);
        }

        *m2->getObjectPtr() = output_array;
        m3->getObjectPtr()->append(GADGETRON_IMAGEPROCESSINGHISTORY, "FFT");

        if (this->next()->putq(m1) < 0)
        {
            GERROR_STREAM("ImageFFTGadget, failed to pass images to next gadget ... ");
            return GADGET_FAIL;
        }

        return GADGET_OK;
    }

    GADGET_FACTORY_DECLARE(ImageFFTGadget)
}
