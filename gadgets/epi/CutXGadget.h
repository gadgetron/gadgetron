#ifndef CutXGADGET_H
#define CutXGADGET_H

#include "Gadget.h"
#include "hoNDArray.h"

#include <complex>

namespace Gadgetron{

  class   CutXGadget :
  public Gadget1<mrd::Acquisition>
  {
    public:
      CutXGadget();
      virtual ~CutXGadget();

    protected:
      virtual int process_config(const mrd::Header& header);
      virtual int process( GadgetContainerMessage< mrd::Acquisition>* m1);

      size_t encodeNx_;
      float encodeFOV_;
      size_t reconNx_;
      float reconFOV_;

      size_t cutNx_;
  };
}
#endif //CutXGADGET_H
