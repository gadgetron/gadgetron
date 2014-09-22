#ifndef CutXGADGET_H
#define CutXGADGET_H

#include "Gadget.h"
#include "hoNDArray.h"
#include "gadgetron_epi_export.h"

#include <ismrmrd/ismrmrd.h>
#include <complex>

namespace Gadgetron{

  class   EXPORTGADGETS_EPI CutXGadget : 
  public Gadget2<ISMRMRD::AcquisitionHeader, hoNDArray< std::complex<float> > >
  {
    public:
      CutXGadget();
      virtual ~CutXGadget();

    protected:
      virtual int process_config(ACE_Message_Block* mb);
      virtual int process( GadgetContainerMessage< ISMRMRD::AcquisitionHeader>* m1,
                       GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2);

      size_t encodeNx_;
      float encodeFOV_;
      size_t reconNx_;
      float reconFOV_;

      size_t cutNx_;
  };
}
#endif //CutXGADGET_H
