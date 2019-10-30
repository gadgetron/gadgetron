#include "GadgetIsmrmrdReader.h"
#include "NHLBICompression.h"

#if defined GADGETRON_COMPRESSION_ZFP

#include <zfp/zfp.h>

#endif //GADGETRON_COMPRESSION_ZFP

namespace Gadgetron {

    Core::Message GadgetIsmrmrdAcquisitionMessageReader::read(std::istream &stream) {

        using namespace Core;
        using namespace std::literals;

        auto header = IO::read<ISMRMRD::AcquisitionHeader>(stream);

        optional<hoNDArray<float>> trajectory = boost::none;
        if (header.trajectory_dimensions) {
            trajectory = hoNDArray<float>(header.trajectory_dimensions,
                                               header.number_of_samples);
            IO::read(stream, trajectory->data(),trajectory->size());
        }
        auto data = hoNDArray<std::complex<float>>(header.number_of_samples,
                                                   header.active_channels);

        if (header.isFlagSet(ISMRMRD::ISMRMRD_ACQ_COMPRESSION1)) { //Is this ZFP compressed data

#if defined GADGETRON_COMPRESSION_ZFP

            uint32_t comp_size = IO::read<uint32_t>(stream);

            std::vector<char> comp_buffer(comp_size);

            stream.read(comp_buffer.data(), comp_size);

            zfp_type type = zfp_type_float;
            auto field = std::unique_ptr<zfp_field, decltype(&zfp_field_free)>(zfp_field_alloc(), &zfp_field_free);

            auto zfp = std::unique_ptr<zfp_stream, decltype(&zfp_stream_close)>(zfp_stream_open(NULL),
                                                                                &zfp_stream_close);
            size_t zfpsize = comp_size;


            auto cstream = std::unique_ptr<bitstream, decltype(&stream_close)>(
                    stream_open(comp_buffer.data(), comp_buffer.size()), &stream_close);

            if (!cstream) {
                throw std::runtime_error("Unable to open compressed stream");
            }

            zfp_stream_set_bit_stream(zfp.get(), cstream.get());

            zfp_stream_rewind(zfp.get());

            if (!zfp_read_header(zfp.get(), field.get(), ZFP_HEADER_FULL)) {
                throw std::runtime_error("Unable to read compressed stream header");

            }

            size_t nx = std::max(field->nx, 1u);
            size_t ny = std::max(field->ny, 1u);
            size_t nz = std::max(field->nz, 1u);

            if (nx * ny * nz != (header.number_of_samples * 2 * header.active_channels)) {
                std::stringstream errorstream;
                errorstream << "Size of decompressed stream does not match the acquisition header";
                errorstream << "nx=" << nx << ", ny=" << ny << ", nz=" << nz;
                errorstream << ", number_of_samples=" << header.number_of_samples;
                errorstream << "active_channels=" << header.active_channels;

                throw std::runtime_error(errorstream.str());
            }

            zfp_field_set_pointer(field.get(), data.get_data_ptr());

            if (!zfp_decompress(zfp.get(), field.get())) {
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
            uint32_t comp_size = IO::read<uint32_t>(stream);

            std::vector<uint8_t> comp_buffer(comp_size);
            stream.read((char *) comp_buffer.data(), comp_size);

            CompressedBuffer<float> comp;
            comp.deserialize(comp_buffer);

            if (comp.size() != data.get_number_of_elements() * 2) { //*2 for complex
                std::stringstream error;
                error << "Mismatch between uncompressed data samples " << comp.size();
                error << " and expected number of samples" << data.get_number_of_elements() * 2;
            }

            float *d_ptr = (float *) data.get_data_ptr();
            for (size_t i = 0; i < comp.size(); i++) {
                d_ptr[i] = comp[i]; //This uncompresses sample by sample into the uncompressed array
            }

            //At this point the data is no longer compressed and we should clear the flag
            header.clearFlag(ISMRMRD::ISMRMRD_ACQ_COMPRESSION2);

        } else {
            //Uncompressed data
            IO::read(stream, data.data(),data.size());

        }

        return Core::Message(std::move(header),std::move(data),std::move(trajectory));

    }

    uint16_t GadgetIsmrmrdAcquisitionMessageReader::slot() {
        return 1008;
    }

    Core::Message GadgetIsmrmrdWaveformMessageReader::read(std::istream &stream) {
        using namespace Core;
        using namespace std::literals;

        auto header = IO::read<ISMRMRD::WaveformHeader>(stream);
        auto data = hoNDArray<uint32_t>(header.number_of_samples, header.channels);

        IO::read(stream, data.data(), data.size());

        return Message(header, std::move(data));
    }


    uint16_t GadgetIsmrmrdWaveformMessageReader::slot() {
        return 1026;
    }

    GADGETRON_READER_EXPORT(GadgetIsmrmrdAcquisitionMessageReader)

    GADGETRON_READER_EXPORT(GadgetIsmrmrdWaveformMessageReader)

}
