#ifndef ISMRMRDDUMPGADGET_H
#define ISMRMRDDUMPGADGET_H

#include "Gadget.h"
#include "hoNDArray.h"
#include "gadgetron_mricore_export.h"

#include <ismrmrd/ismrmrd.h>
#include <ismrmrd/dataset.h>
#include <ismrmrd/xml.h>

#include <complex>

namespace Gadgetron {

    class EXPORTGADGETSMRICORE IsmrmrdDumpGadget :
        public Gadgetron::Gadget2<ISMRMRD::AcquisitionHeader, hoNDArray< std::complex<float> > >
    {
    public:
        typedef Gadgetron::Gadget2<ISMRMRD::AcquisitionHeader, hoNDArray< std::complex<float> > > BaseClass;

        GADGET_DECLARE(IsmrmrdDumpGadget);

        IsmrmrdDumpGadget();
        virtual ~IsmrmrdDumpGadget();

    protected:

#ifdef WIN32
        GADGET_PROPERTY(folder, std::string, "Folder for save dump file", "/tmp/gadgetron_data");
#else
        GADGET_PROPERTY(folder, std::string, "Folder for save dump file", "c:/temp/gadgetron_data");
#endif // WIN32

        GADGET_PROPERTY(file_prefix, std::string, "Prefix for dump file", "ISMRMRD_DUMP");

        // In some cases, data cannot be saved to a gadgetron server
        // if gadgetron ip equals to these preset ips, do not save data
        GadgetProperty<std::vector<std::string>, GadgetPropertyLimitsNoLimits<std::vector<std::string> > > ip_no_data_saving{ "ip_no_data_saving",
            "std::string",
            "If gadgetrion IP equals to this ip, do not save data",
            this,
            { "192.168.2.2", "192.168.56.2" },
            GadgetPropertyLimitsNoLimits<std::vector<std::string> >() };

        // if true, only save the xml header
        GADGET_PROPERTY(save_xml_header_only, bool, "If true, only save the xml header", false);

        virtual int process_config(ACE_Message_Block* mb);
        virtual int process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1, GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2);

    private:

        bool first_call_;
        bool save_ismrmrd_data_;
        ISMRMRD::IsmrmrdHeader ismrmrd_header_;
        std::string ismrmrd_xml_;
        boost::shared_ptr<ISMRMRD::Dataset>  ismrmrd_dataset_;

        int create_ismrmrd_dataset(ISMRMRD::AcquisitionHeader* acq = NULL);
    };
}
#endif //ISMRMRDDUMPGADGET_H
