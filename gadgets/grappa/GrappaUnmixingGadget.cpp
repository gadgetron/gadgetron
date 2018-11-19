#include "GadgetIsmrmrdReadWrite.h"
#include "GrappaUnmixingGadget.h"
#include "hoNDFFT.h"
namespace Gadgetron{

  GrappaUnmixingGadget::GrappaUnmixingGadget() {
    // TODO Auto-generated constructor stub

  }

  GrappaUnmixingGadget::~GrappaUnmixingGadget() {
    // TODO Auto-generated destructor stub
  }

  int GrappaUnmixingGadget::process(GadgetContainerMessage<GrappaUnmixingJob>* m1,
                                    GadgetContainerMessage<ISMRMRD::ImageHeader>* m2, GadgetContainerMessage<hoNDArray<std::complex<float> > >* m3)
  {
    GadgetContainerMessage< hoNDArray<std::complex<float> > >* cm2 =
			new GadgetContainerMessage< hoNDArray<std::complex<float> > >();

    std::vector<size_t> combined_dims(3,0);
    combined_dims[0] = m2->getObjectPtr()->matrix_size[0];
    combined_dims[1] = m2->getObjectPtr()->matrix_size[1];
    combined_dims[2] = m2->getObjectPtr()->matrix_size[2];

    if (m2->getObjectPtr()->channels > 1) {
      combined_dims.push_back(m2->getObjectPtr()->channels);
    }

    try{cm2->getObjectPtr()->create(&combined_dims);}
    catch (std::runtime_error &err ){
      GEXCEPTION(err,"Unable to create combined image array\n");
      return GADGET_FAIL;
    }

    m1->cont(0);
    m2->cont(cm2);


    hoNDFFT<float>::instance()->ifft3c(*m3->getObjectPtr());
    /*
    hoNDFFT<float>::instance()->ifft(m3->getObjectPtr(),0);
    hoNDFFT<float>::instance()->ifft(m3->getObjectPtr(),1);
    hoNDFFT<float>::instance()->ifft(m3->getObjectPtr(),2);
    */

    if (!m1->getObjectPtr()->weights_) {
      GDEBUG("Weights are a NULL\n");
      return GADGET_FAIL;
    }

    float scale_factor = 1.0;
    int appl_result = m1->getObjectPtr()->weights_->apply(m3->getObjectPtr(), cm2->getObjectPtr(), scale_factor);
    if (appl_result < 0) {
      GDEBUG("Failed to apply GRAPPA weights: error code %d\n", appl_result);
      return GADGET_FAIL;
    }

    m1->release();
    m3->release();

    if (this->next()->putq(m2) < 0) {
      return GADGET_FAIL;
    }

    return GADGET_OK;
  }

  GADGET_FACTORY_DECLARE(GrappaUnmixingGadget)
}
