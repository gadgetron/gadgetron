#ifndef ISMRMRDDUMPGADGET_H
#define ISMRMRDDUMPGADGET_H

#include "Gadget.h"
#include "hoNDArray.h"
#include "gadgetron_mricore_export.h"

#include <ismrmrd/ismrmrd.h>
#include <ismrmrd/dataset.h>
#include <ismrmrd/xml.h>

#include <complex>

namespace Gadgetron{

  class EXPORTGADGETSMRICORE IsmrmrdDumpGadget : 
  public Gadgetron::Gadget2<ISMRMRD::AcquisitionHeader,hoNDArray< std::complex<float> > >
    {
    public:
      GADGET_DECLARE(IsmrmrdDumpGadget);

      IsmrmrdDumpGadget();
      virtual ~IsmrmrdDumpGadget();

    protected:
      GADGET_PROPERTY(folder,               std::string,    "Folder  for dump file", ".");
      GADGET_PROPERTY(file_prefix,          std::string,    "Prefix for dump file", "ISMRMRD_DUMP");
      GADGET_PROPERTY(append_id,            bool,           "ISMRMRD measurement ID to file name prefix (if available)", true);
      GADGET_PROPERTY(append_timestamp,     bool,           "Append timestamp to file name prefix", true);
      GADGET_PROPERTY(use_data_timestamp,   bool,           "Use the time stamp of incoming data set as file name", false);
      GADGET_PROPERTY(timestamp_tick,       double,         "Tick per unit of time stamp in ms", 2.5);

      virtual int process_config(ACE_Message_Block* mb);

      virtual int process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1,
                GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2);

    private:
      bool first_call_;
      ISMRMRD::IsmrmrdHeader ismrmrd_header_;
      std::string ismrmrd_xml_;
      boost::shared_ptr<ISMRMRD::Dataset>  ismrmrd_dataset_;

      int create_ismrmrd_dataset(ISMRMRD::AcquisitionHeader* acq = NULL);
    };
}
#endif //ISMRMRDDUMPGADGET_H
