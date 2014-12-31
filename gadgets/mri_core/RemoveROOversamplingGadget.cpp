#include "RemoveROOversamplingGadget.h"
#include "Gadgetron.h"
#include "hoNDFFT.h"
#include "ismrmrd/xml.h"

#ifdef USE_OMP
    #include "omp.h"
#endif // USE_OMP

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

	ISMRMRD::IsmrmrdHeader h;
	ISMRMRD::deserialize(mb->rd_ptr(),h);

	if (h.encoding.size() == 0) {
	  GDEBUG("Number of encoding spaces: %d\n", h.encoding.size());
	  GDEBUG("This Gadget needs an encoding description\n");
	  return GADGET_FAIL;
	}

	
	ISMRMRD::EncodingSpace e_space = h.encoding[0].encodedSpace;
	ISMRMRD::EncodingSpace r_space = h.encoding[0].reconSpace;

        encodeNx_  = e_space.matrixSize.x;
        encodeFOV_ = e_space.fieldOfView_mm.x;
        reconNx_   = r_space.matrixSize.x;
        reconFOV_  = r_space.fieldOfView_mm.x;

        // limit the number of threads used to be 1
#ifdef USE_OMP
        omp_set_num_threads(1);
        GADGET_MSG("RemoveROOversamplingGadget:omp_set_num_threads(1) ... ");
#endif // USE_OMP

    // If the encoding and recon matrix size and FOV are the same
    // then the data is not oversampled and we can safely pass
    // the data onto the next gadget
    if ( (encodeNx_ == reconNx_) && (encodeFOV_ == reconFOV_) )
    {
      dowork_ = false;
    }
    else {
      dowork_ = true;
    }

        return GADGET_OK;
    }

    int RemoveROOversamplingGadget
        ::process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1,
        GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2)
    {

      // If we have work to do, do it, otherwise do nothing
      if (dowork_) {

        GadgetContainerMessage< hoNDArray< std::complex<float> > >* m3 
            = new GadgetContainerMessage< hoNDArray< std::complex<float> > >();

        if (!m3)
        {
            return GADGET_FAIL;
        }

        std::vector<size_t> data_out_dims = *m2->getObjectPtr()->get_dimensions();
        if ( !ifft_buf_.dimensions_equal(&data_out_dims) )
        {
            ifft_buf_.create(data_out_dims);
            ifft_res_.create(data_out_dims);
        }

        float ratioFOV = encodeFOV_/reconFOV_;

        data_out_dims[0] = (size_t)(data_out_dims[0]/ratioFOV);
        if ( !fft_buf_.dimensions_equal(&data_out_dims) )
        {
            fft_buf_.create(data_out_dims);
            fft_res_.create(data_out_dims);
        }

        try{ m3->getObjectPtr()->create(&data_out_dims);}
        catch (std::runtime_error &err)
        {
            GADGET_DEBUG_EXCEPTION(err,"Unable to create new data array for downsampled data\n");
            return GADGET_FAIL;
        }

        size_t sRO = m2->getObjectPtr()->get_size(0);
        size_t start = (size_t)( (m2->getObjectPtr()->get_size(0)-data_out_dims[0])/ratioFOV );

        size_t dRO = m3->getObjectPtr()->get_size(0);
        size_t numOfBytes = data_out_dims[0]*sizeof(std::complex<float>);

        int c;

        int CHA = (int)(data_out_dims[1]);

        std::complex<float>* data_in, *data_out;

        if ( constant_noise_variance_ )
        {
            hoNDFFT<float>::instance()->ifft1c(*m2->getObjectPtr(), ifft_res_, ifft_buf_);

            data_in  = ifft_res_.get_data_ptr();
            data_out = fft_res_.get_data_ptr();

            // #pragma omp parallel for default(none) private(c) shared(CHA, sRO, start, dRO, data_in, data_out, numOfBytes)
            for ( c=0; c<CHA; c++)
            {
                memcpy( data_out+c*dRO, data_in+c*sRO+start, numOfBytes );
            }

            hoNDFFT<float>::instance()->fft1c(fft_res_, *m3->getObjectPtr(), fft_buf_);
        }
        else
        {
            hoNDFFT<float>::instance()->ifft(m2->getObjectPtr(), 0);
            data_in  = m2->getObjectPtr()->get_data_ptr();
            data_out = m3->getObjectPtr()->get_data_ptr();

            // #pragma omp parallel for default(none) private(c) shared(CHA, sRO, start, dRO, data_in, data_out, numOfBytes)
            for ( c=0; c<CHA; c++)
            {
                memcpy( data_out+c*dRO, data_in+c*sRO+start, numOfBytes );
            }

            hoNDFFT<float>::instance()->fft(m3->getObjectPtr(), 0);
        }

        m2->release(); //We are done with this data

        m1->cont(m3);
        m1->getObjectPtr()->number_of_samples = data_out_dims[0];
        m1->getObjectPtr()->center_sample = (uint16_t)(m1->getObjectPtr()->center_sample/ratioFOV);

      } // end if (dowork_)

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
