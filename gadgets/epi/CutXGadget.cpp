#include "CutXGadget.h"
#include "hoNDFFT.h"
#include "hoNDArray_utils.h"

namespace Gadgetron{

    CutXGadget::CutXGadget() {}
    CutXGadget::~CutXGadget() {}

    int CutXGadget::process_config(const mrd::Header& header)
    {
        auto& h = header;

      if (h.encoding.size() == 0) {
	GDEBUG("Number of encoding spaces: %d\n", h.encoding.size());
	GDEBUG("This Gadget needs an encoding description\n");
	return GADGET_FAIL;
      }

      // Get the encoding space and trajectory description
      mrd::EncodingSpaceType e_space = h.encoding[0].encoded_space;
      mrd::EncodingSpaceType r_space = h.encoding[0].recon_space;
      mrd::EncodingLimitsType e_limits = h.encoding[0].encoding_limits;
      mrd::TrajectoryDescriptionType traj_desc;

      // Primary encoding space is for EPI
      encodeNx_  = e_space.matrix_size.x;
      encodeFOV_ = e_space.field_of_view_mm.x;
      reconNx_   = r_space.matrix_size.x;
      reconFOV_  = r_space.field_of_view_mm.x;

      cutNx_ = encodeNx_;

      return 0;
    }

    int CutXGadget::process( GadgetContainerMessage< mrd::Acquisition>* m1)
    {
        try
        {
            // cut the central half from the kspace line
            if ( m1->getObjectPtr()->Samples() > cutNx_ )
            {
                size_t RO = m1->getObjectPtr()->Samples();

                auto& head = m1->getObjectPtr()->head;

                uint32_t center_sample = head.center_sample.value_or(RO/2);
                uint16_t startX = center_sample - cutNx_/2;
                uint16_t endX = startX + cutNx_ - 1;

                float ratio = RO / (float)cutNx_;
                head.center_sample = (uint16_t)(center_sample / ratio );

                GadgetContainerMessage< hoNDArray< std::complex<float> > >* m3 = new GadgetContainerMessage< hoNDArray< std::complex<float> > >();

                auto& data = m1->getObjectPtr()->data;
                std::vector<size_t> dim(2);
                dim[0] = cutNx_;
                dim[1] = data.get_size(1);

                m3->getObjectPtr()->create(dim);

                size_t cha;
                for ( cha=0; cha<dim[1]; cha++ )
                {
                    memcpy(m3->getObjectPtr()->begin()+cha*cutNx_,
                            data.begin()+cha*RO+startX,
                            sizeof( std::complex<float> )*cutNx_);
                }

                m1->cont(m3);
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
