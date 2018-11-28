#pragma  once
#include "GadgetMRIHeaders.h"
#include "GadgetContainerMessage.h"
#include "GadgetMessageInterface.h"
#include "hoNDArray.h"
#include "url_encode.h"
#include "gadgetron_mricore_export.h"
#include <ismrmrd/ismrmrd.h>
#include <ismrmrd/waveform.h>
#include <complex>
#include "Stream.h"
#include "mri_core_data.h"
#include "io/readers.h"
#include "io/writers.h"

#include "NHLBICompression.h"

#if defined GADGETRON_COMPRESSION_ZFP
#include <zfp/zfp.h>
#endif //GADGETRON_COMPRESSION_ZFP



namespace Gadgetron {

    class EXPORTGADGETSMRICORE GadgetIsmrmrdAcquisitionMessageWriter : public Core::Writer {

    public:
        virtual void write(std::ostream &stream, std::unique_ptr<Core::Message> &&message) override final {
            using namespace Core;
            auto hm = dynamic_cast<Core::TypedMessage<Acquisition>*>(message.get());

            if (!hm) {
               throw std::runtime_error("GadgetAcquisitionMessageWriter, invalid acquisition message objects");
            }

            auto acquisition = hm->get_data();

            ssize_t send_cnt = 0;

            GadgetMessageIdentifier id;
            id.id = GADGET_MESSAGE_ISMRMRD_ACQUISITION;

            IO::write(stream,id);

            ISMRMRD::AcquisitionHeader& acqHead = acquisition->header;

            IO::write(stream,acqHead);


            unsigned long trajectory_elements = acqHead.trajectory_dimensions * acqHead.number_of_samples;
            unsigned long data_elements = acqHead.active_channels * acqHead.number_of_samples;


            if (acquisition->trajectory) {
                IO::write(stream,*acquisition->trajectory);
            }

            if (data_elements) {
                IO::write(stream,acquisition->data);
            }
        }

        virtual std::vector<std::type_index> supported_types() const override final {
            return {std::type_index(typeid(Core::TypedMessage<Acquisition>))};
        }
    };

    /**
    Default implementation of GadgetMessageReader for IsmrmrdAcquisition messages
    */
class EXPORTGADGETSMRICORE GadgetIsmrmrdAcquisitionMessageReader : public Core::Reader {

    public:

        virtual std::unique_ptr<Core::Message> read(std::istream& stream) {

            using namespace Core;
            using namespace std::literals;
            auto acquisition = std::make_unique<Acquisition>();
            auto& header = acquisition->header;

            IO::read(stream,header);

            if (header.trajectory_dimensions) {

                acquisition->trajectory = hoNDArray<float>( header.trajectory_dimensions, header.number_of_samples);
                IO::read(stream,*acquisition->trajectory);
            }

            acquisition->data = hoNDArray<std::complex<float>>(header.number_of_samples,header.active_channels);



            if (header.isFlagSet(ISMRMRD::ISMRMRD_ACQ_COMPRESSION1)) { //Is this ZFP compressed data

#if defined GADGETRON_COMPRESSION_ZFP

                uint32_t comp_size = 0;
                IO::read(stream, comp_size);

                std::vector<char> comp_buffer(comp_size);

                stream.read(comp_buffer.data(),comp_size);

                zfp_type type = zfp_type_float;
                auto  field = std::unique_ptr<zfp_field,void(zfp_field*)>(zfp_field_alloc(),zfp_field_free);

                auto zfp = std::unique_ptr<zfp_stream,void(zfp_stream*)>(zfp_stream_open(NULL),zfp_stream_close);
                size_t zfpsize = comp_size;
                

                auto cstream = std::unique_ptr<bitstream,void(bitstream*)>(
                        stream_open(comp_buffer.data(), comp_buffer.size()), stream_close);

                if (!cstream) {
                    throw std::runtime_error("Unable to open compressed stream");
                }

                zfp_stream_set_bit_stream(zfp.get(), cstream.get());
                
                zfp_stream_rewind(zfp.get());

                if (!zfp_read_header(zfp, field, ZFP_HEADER_FULL)) {
                    throw std::runtime_error("Unable to read compressed stream header");

                }
                
                size_t nx = std::max(field->nx, 1u);
                size_t ny = std::max(field->ny, 1u);
                size_t nz = std::max(field->nz, 1u);
                
                if (nx*ny*nz != (header.number_of_samples*2*header.active_channels)) {
                    std::stringstream errorstream;
                    errorstream << "Size of decompressed stream does not match the acquisition header";
                    errorstream << "nx=" << nx << ", ny=" << ny << ", nz=" << nz;
                    errorstream << ", number_of_samples=" << header.number_of_samples;
                    errorstream << "active_channels=" << header.active_channels;

                    throw std::runtime_error(errorstream.str());
                }

                zfp_field_set_pointer(field.get(), acquisition->data.get_data_ptr());
                
                if (!zfp_decompress(zfp, field)) {
                    throw std::runtime_error("Unable to decompress stream");
                }


                //At this point the data is no longer compressed and we should clear the flag
                header.clearFlag(ISMRMRD::ISMRMRD_ACQ_COMPRESSION1);

#else //GADGETRON COMPRESSION_ZFP

                //This is compressed data, but Gadgetron was not compiled with compression
                throw std::runtime_error("Receiving compressed (ZFP) data, but Gadgetron was not compiled with ZFP support");

#endif //GADGETRON_COMPRESSION_ZFP

            } else if (header.isFlagSet(ISMRMRD::ISMRMRD_ACQ_COMPRESSION2)) {
                //NHLBI Compression
                uint32_t comp_size;
                IO::read(stream,comp_size);

                std::vector<uint8_t> comp_buffer(comp_size);
                stream.read((char*)comp_buffer.data(),comp_size);

                CompressedBuffer<float> comp;
                comp.deserialize(comp_buffer);

                if (comp.size() != acquisition->data.get_number_of_elements() * 2) { //*2 for complex
                    std::stringstream error;
                    error << "Mismatch between uncompressed data samples " << comp.size();
                    error << " and expected number of samples" << acquisition->data.get_number_of_elements()*2;
                }

                float *d_ptr = (float *) acquisition->data.get_data_ptr();
                for (size_t i = 0; i < comp.size(); i++) {
                    d_ptr[i] = comp[i]; //This uncompresses sample by sample into the uncompressed array
                }

                //At this point the data is no longer compressed and we should clear the flag
                header.clearFlag(ISMRMRD::ISMRMRD_ACQ_COMPRESSION2);

            } else {
                //Uncompressed data
                IO::read(stream,acquisition->data);

            }

            return std::unique_ptr<Message>(new TypedMessage<Acquisition>(std::move(acquisition)));
        }

    };

    // ------------------------------------------------------------------------------------------------------- //
    // ISMRMRD wave form reader/writer

class EXPORTGADGETSMRICORE GadgetIsmrmrdWaveformMessageWriter : public Core::Writer {

    public:
        virtual void  write(std::ostream& stream, std::unique_ptr<Core::Message>&& message ) override final {

            using namespace Core;
            auto hm = dynamic_cast<Core::TypedMessage<Waveform>*>(message.get());

            if (!hm) {
               throw std::runtime_error("GadgetAcquisitionMessageWriter, invalid acquisition message objects");
            }
            auto h = hm->get_data();

            GadgetMessageIdentifier id;
            id.id = GADGET_MESSAGE_ISMRMRD_WAVEFORM;
            IO::write(stream,id);

            ISMRMRD::ISMRMRD_WaveformHeader& wavHead = h->header;
            IO::write(stream,wavHead);

            IO::write(stream,h->data);
        }
    };

class EXPORTGADGETSMRICORE GadgetIsmrmrdWaveformMessageReader : public Core::Reader {

    public:


        virtual std::unique_ptr<Core::Message> read(std::istream& stream) override final {
            using namespace Core;
            auto wave = std::make_unique<Waveform>();

            IO::read(stream,wave->header);

            wave->data = hoNDArray<uint32_t>(wave->header.number_of_samples,wave->header.channels);

            IO::read(stream,wave->data);

            return std::unique_ptr<Message>(new TypedMessage<Waveform>(std::move(wave)));

        }
    };
}
