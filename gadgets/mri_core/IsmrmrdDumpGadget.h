#ifndef ISMRMRDDUMPGADGET_H
#define ISMRMRDDUMPGADGET_H

#include "Gadget.h"
#include "hoNDArray.h"
#include "gadgetron_mricore_export.h"

#include <ismrmrd/ismrmrd.h>
#include <ismrmrd/dataset.h>

#include <complex>

namespace Gadgetron{

  class EXPORTGADGETSMRICORE IsmrmrdDumpGadget : 
  public Gadgetron::Gadget2<ISMRMRD::AcquisitionHeader,hoNDArray< std::complex<float> > >
    {
    public:
      GADGET_DECLARE(IsmrmrdDumpGadget);

      IsmrmrdDumpGadget();

    protected:
      virtual int process_config(ACE_Message_Block* mb);

      virtual int process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1,
			  GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2);

    private:
      std::string file_prefix_;
      std::string ismrmrd_file_name_;
      boost::shared_ptr<ISMRMRD::Dataset>  ismrmrd_dataset_;
      bool append_timestamp_;
    };
}
#endif //ISMRMRDDUMPGADGET_H
