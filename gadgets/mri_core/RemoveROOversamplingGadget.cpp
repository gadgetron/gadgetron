#include "GadgetIsmrmrdReadWrite.h"
#include "RemoveROOversamplingGadget.h"
#include "Gadgetron.h"
#include "hoNDFFT.h"

namespace Gadgetron{

    RemoveROOversamplingGadget::RemoveROOversamplingGadget() : constant_noise_variance_(false)
    {
    }

    RemoveROOversamplingGadget::~RemoveROOversamplingGadget()
    {
    }

    int RemoveROOversamplingGadget::process_config(ACE_Message_Block* mb)
    {
        constant_noise_variance_ = this->get_bool_value("constant_noise_variance");
        return GADGET_OK;
    }

    int RemoveROOversamplingGadget
        ::process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1,
        GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2)
    {
        GadgetContainerMessage< hoNDArray< std::complex<float> > >* m3 
            = new GadgetContainerMessage< hoNDArray< std::complex<float> > >();

        if (!m3)
        {
            return GADGET_FAIL;
        }

        std::vector<size_t> data_out_dims = *m2->getObjectPtr()->get_dimensions();
        data_out_dims[0] = data_out_dims[0]/2;

        try{ m3->getObjectPtr()->create(&data_out_dims);}
        catch (std::runtime_error &err)
        {
            GADGET_DEBUG_EXCEPTION(err,"Unable to create new data array for downsampled data\n");
            return GADGET_FAIL;
        }

        if ( constant_noise_variance_ )
        {
            hoNDFFT<float>::instance()->ifft1c(*m2->getObjectPtr());
        }
        else
        {
            hoNDFFT<float>::instance()->ifft(m2->getObjectPtr(), 0);
        }

        std::complex<float>* data_in  = m2->getObjectPtr()->get_data_ptr();
        std::complex<float>* data_out = m3->getObjectPtr()->get_data_ptr();

        int c;

        int CHA = (int)(data_out_dims[1]);

        size_t sRO = m2->getObjectPtr()->get_size(0);
        size_t start = (m2->getObjectPtr()->get_size(0)-data_out_dims[0])/2;

        size_t dRO = m3->getObjectPtr()->get_size(0);
        size_t numOfBytes = data_out_dims[0]*sizeof(std::complex<float>);

        #pragma omp parallel for default(none) private(c) shared(CHA, sRO, start, dRO, data_in, data_out, numOfBytes)
        for ( c=0; c<CHA; c++)
        {
            memcpy( data_out+c*dRO, data_in+c*sRO+start, numOfBytes );
        }

        if ( constant_noise_variance_ )
        {
            hoNDFFT<float>::instance()->fft1c(*m3->getObjectPtr());
        }
        else
        {
            hoNDFFT<float>::instance()->fft(m3->getObjectPtr(), 0);
        }

        m2->release(); //We are done with this data

        m1->cont(m3);
        m1->getObjectPtr()->number_of_samples = data_out_dims[0];
        m1->getObjectPtr()->center_sample /= 2;

        if (this->next()->putq(m1) == -1)
        {
            ACE_ERROR_RETURN( (LM_ERROR,
                ACE_TEXT("%p\n"),
                ACE_TEXT("RemoveROOversamplingGadget::process, passing data on to next gadget")),
                GADGET_FAIL);
        }

        return GADGET_OK;
    }


    GADGET_FACTORY_DECLARE(RemoveROOversamplingGadget)
}
