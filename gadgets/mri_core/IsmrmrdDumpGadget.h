#pragma once 
#include "Node.h"
#include "Types.h"
#include "hoNDArray.h"
#include "io/from_string.h"

#include <ismrmrd/ismrmrd.h>
#include <ismrmrd/dataset.h>
#include <ismrmrd/xml.h>
#include <ismrmrd/waveform.h>
#include <set>

#include <complex>



namespace Gadgetron {
       class IsmrmrdDumpGadget : public Core::ChannelGadget<Core::variant<Core::Acquisition, Core::Waveform>>
    {
    public:
      
        IsmrmrdDumpGadget(const Core::Context& context, const Core::GadgetProperties& props ); 
        virtual ~IsmrmrdDumpGadget() = default;

    protected:

#ifndef WIN32
        NODE_PROPERTY(folder, boost::filesystem::path, "Folder for save dump file", "/tmp/gadgetron_data");
#else
        NODE_PROPERTY(folder, boost::filesystem::path, "Folder for save dump file", "c:/temp/gadgetron_data");
#endif // WIN32

        NODE_PROPERTY(file_prefix, std::string, "Prefix for dump file", "ISMRMRD_DUMP");

        NODE_PROPERTY(env_var_to_control_dump, std::string,
                        "Environmental variable to control the dump, if empty, dump the data",
                        "GADGETRON_ISMRMRD_DUMP");

        // In some cases, data cannot be saved to a gadgetron server
        // if gadgetron ip equals to these preset ips, do not save data
        NODE_PROPERTY(ip_no_data_saving, std::set<std::string>, "If gadgetrion IP equals to this ip, do not save data",std::set<std::string>({ "192.168.2.2", "192.168.56.2" }));

        // if true, only save the xml header
        NODE_PROPERTY(save_xml_header_only, bool, "If true, only save the xml header", false);

        void process(Core::InputChannel<Core::variant<Core::Acquisition,Core::Waveform>>& input, Core::OutputChannel& output) override;

    private:

        bool save_ismrmrd_data_;
        ISMRMRD::Dataset create_ismrmrd_dataset() const;
        bool  is_ip_on_blacklist() const ; 
    };
}
