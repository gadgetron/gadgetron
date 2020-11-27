#ifndef ISMRMRDDUMPGADGET_H
#define ISMRMRDDUMPGADGET_H

#include "Gadget.h"
#include "hoNDArray.h"
#include "gadgetron_mricore_export.h"

#include <ismrmrd/ismrmrd.h>
#include <ismrmrd/dataset.h>
#include <ismrmrd/xml.h>
#include <ismrmrd/waveform.h>

#include <complex>

namespace Gadgetron {

    class EXPORTGADGETSMRICORE IsmrmrdDumpGadget :
        public Gadgetron::Gadget1Of2<ISMRMRD::AcquisitionHeader, ISMRMRD::WaveformHeader >
    {
    public:
        typedef Gadgetron::Gadget1Of2<ISMRMRD::AcquisitionHeader, ISMRMRD::WaveformHeader > BaseClass;

        GADGET_DECLARE(IsmrmrdDumpGadget);

        IsmrmrdDumpGadget();
        virtual ~IsmrmrdDumpGadget();

    protected:

#ifndef WIN32
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

        // since the support to waveform is not fully implemented, this option is added for not passing waveform downstream
        // TODO: remove this option
        GADGET_PROPERTY(pass_waveform_downstream, bool, "If true, waveform data is passed downstream", false);

        int process_config(ACE_Message_Block* mb) override;
        int process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1) override ;
        int process(GadgetContainerMessage<ISMRMRD::WaveformHeader>* m1) override ;

    private:

        bool first_call_;
        bool save_ismrmrd_data_;
        ISMRMRD::IsmrmrdHeader ismrmrd_header_;
        std::string ismrmrd_xml_;
        boost::shared_ptr<ISMRMRD::Dataset>  ismrmrd_dataset_;

        int create_ismrmrd_dataset();
    };
}
#endif //ISMRMRDDUMPGADGET_H
