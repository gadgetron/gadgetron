#include "CutXGadget.h"
#include "hoNDFFT.h"
#include "hoNDArray_utils.h"
#include "ismrmrd/xml.h"

namespace Gadgetron{

    CutXGadget::CutXGadget() {}
    CutXGadget::~CutXGadget() {}

    int CutXGadget::process_config(ACE_Message_Block* mb)
    {
      ISMRMRD::IsmrmrdHeader h;
      ISMRMRD::deserialize(mb->rd_ptr(),h);
      
      
      if (h.encoding.size() == 0) {
	GDEBUG("Number of encoding spaces: %d\n", h.encoding.size());
	GDEBUG("This Gadget needs an encoding description\n");
	return GADGET_FAIL;
      }
      
      // Get the encoding space and trajectory description
      ISMRMRD::EncodingSpace e_space = h.encoding[0].encodedSpace;
      ISMRMRD::EncodingSpace r_space = h.encoding[0].reconSpace;
      ISMRMRD::EncodingLimits e_limits = h.encoding[0].encodingLimits;
      ISMRMRD::TrajectoryDescription traj_desc;
      
      // Primary encoding space is for EPI
      encodeNx_  = e_space.matrixSize.x;
      encodeFOV_ = e_space.fieldOfView_mm.x;
      reconNx_   = r_space.matrixSize.x;
      reconFOV_  = r_space.fieldOfView_mm.x;
      
      cutNx_ = encodeNx_;

      return 0;
    }

    int CutXGadget::process( GadgetContainerMessage< ISMRMRD::AcquisitionHeader>* m1,
        GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2)
    {
        try
        {
            // cut the central half from the kspace line
            if ( m1->getObjectPtr()->number_of_samples > cutNx_ )
            {
                size_t RO = m1->getObjectPtr()->number_of_samples;

                uint16_t startX = m1->getObjectPtr()->center_sample - cutNx_/2;
                uint16_t endX = startX + cutNx_ - 1;

                float ratio = RO / (float)cutNx_;
                m1->getObjectPtr()->number_of_samples = cutNx_;
                m1->getObjectPtr()->center_sample = (uint16_t)(m1->getObjectPtr()->center_sample / ratio );

                GadgetContainerMessage< hoNDArray< std::complex<float> > >* m3 = new GadgetContainerMessage< hoNDArray< std::complex<float> > >();

                std::vector<size_t> dim(2);
                dim[0] = cutNx_;
                dim[1] = m2->getObjectPtr()->get_size(1);

                m3->getObjectPtr()->create(dim);

                size_t cha;
                for ( cha=0; cha<dim[1]; cha++ )
                {
                    memcpy(m3->getObjectPtr()->begin()+cha*cutNx_, 
                            m2->getObjectPtr()->begin()+cha*RO+startX, 
                            sizeof( std::complex<float> )*cutNx_);
                }

                m1->cont(m3);
                m2->release();
            }

            if (this->next()->putq(m1) < 0)
            {
                return GADGET_FAIL;
            }
        }
        catch(...)
        {
            return GADGET_FAIL;
        }

        return GADGET_OK;
    }

    GADGET_FACTORY_DECLARE(CutXGadget)
}
