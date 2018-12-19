#pragma  once
#include "GadgetMRIHeaders.h"
#include "GadgetContainerMessage.h"
#include "hoNDArray.h"
#include "url_encode.h"
#include "gadgetron_mricore_export.h"
#include <ismrmrd/ismrmrd.h>
#include <ismrmrd/waveform.h>
#include <complex>
#include "Reader.h"
#include "Writer.h"
#include "mri_core_data.h"
#include "io/readers.h"
#include "io/writers.h"





namespace Gadgetron {

    class EXPORTGADGETSMRICORE GadgetIsmrmrdAcquisitionMessageWriter : public Core::Writer {

    public:
//        virtual void write(std::ostream &stream, std::unique_ptr<Core::Message> &&message) override final {
//            using namespace Core;
//            auto hm = dynamic_cast<Core::TypedMessage<Acquisition>*>(message.get());
//
//            if (!hm) {
//               throw std::runtime_error("GadgetAcquisitionMessageWriter, invalid acquisition message objects");
//            }
//
//            auto acquisition = hm->get_data();
//
//            ssize_t send_cnt = 0;
//
//            GadgetMessageIdentifier id;
//            id.id = GADGET_MESSAGE_ISMRMRD_ACQUISITION;
//
//            IO::write(stream,id);
//
//            ISMRMRD::AcquisitionHeader& acqHead = acquisition->header;
//
//            IO::write(stream,acqHead);
//
//
//            unsigned long trajectory_elements = acqHead.trajectory_dimensions * acqHead.number_of_samples;
//            unsigned long data_elements = acqHead.active_channels * acqHead.number_of_samples;
//
//
//            if (acquisition->trajectory) {
//                IO::write(stream,*acquisition->trajectory);
//            }
//
//            if (data_elements) {
//                IO::write(stream,acquisition->data);
//            }
//        }
//
//
//        virtual std::vector<std::type_index> supported_types() const override final {
//            return {std::type_index(typeid(Core::TypedMessage<Acquisition>))};
//        }
    };

    /**
    Default implementation of GadgetMessageReader for IsmrmrdAcquisition messages
    */
class EXPORTGADGETSMRICORE GadgetIsmrmrdAcquisitionMessageReader : public Core::Reader {

    public:

        virtual std::unique_ptr<Core::Message> read(std::istream& stream) override final;


        virtual uint16_t slot() override final;
        virtual ~GadgetIsmrmrdAcquisitionMessageReader(){};
    };

    // ------------------------------------------------------------------------------------------------------- //
    // ISMRMRD wave form reader/writer



class EXPORTGADGETSMRICORE GadgetIsmrmrdWaveformMessageReader : public Core::Reader {

    public:

        virtual std::unique_ptr<Core::Message> read(std::istream& stream) override final;
        virtual uint16_t slot() override final;
        virtual ~GadgetIsmrmrdWaveformMessageReader(){};
    };
}
