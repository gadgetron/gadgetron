#include <string>
#include <memory>
#include <mutex>

#include <boost/asio.hpp>
#include <boost/filesystem.hpp>

#include <ismrmrd/ismrmrd.h>
#include <ismrmrd/dataset.h>
#include <ismrmrd/meta.h>
#include <ismrmrd/xml.h>
#include <ismrmrd/waveform.h>

#include "MessageID.h"


class IsmrmrdMsgAdapter
{
public:
    IsmrmrdMsgAdapter(){}
    ~IsmrmrdMsgAdapter(){}

    void convert(
        std::ostream& stream,
        const std::string& ismrmrd_filepath,
        const std::string& input_group = "/dataset")
    {
        auto dataset = std::make_shared<ISMRMRD::Dataset>(ismrmrd_filepath.c_str(), input_group.c_str(), false);

        write_ismrmrd_hdr(stream, dataset);
        write_data(stream, dataset);
        write_close(stream);
    }

private:

    void write_ismrmrd_hdr(std::ostream& stream, std::shared_ptr<ISMRMRD::Dataset> dataset)
    {
        std::string raw_ismrmrd_header;
        dataset->readHeader(raw_ismrmrd_header);
        uint32_t size = raw_ismrmrd_header.size();

        const auto id = ::Gadgetron::Core::MessageID::HEADER;
        stream.write(reinterpret_cast<const char*>(&id), sizeof(::Gadgetron::Core::MessageID));
        stream.write(reinterpret_cast<const char*>(&size), sizeof(uint32_t));
        stream.write(raw_ismrmrd_header.c_str(), raw_ismrmrd_header.size());
    }

    void write_data(std::ostream& stream, std::shared_ptr<ISMRMRD::Dataset> dataset)
    {
        auto acquisitions = dataset->getNumberOfAcquisitions();
        auto waveforms = dataset->getNumberOfWaveforms();

        ISMRMRD::Acquisition acq_tmp;
        ISMRMRD::Waveform wav_tmp;

        uint32_t i(0), j(0);

        if(waveforms>0)
        {
            dataset->readAcquisition(i, acq_tmp);
            dataset->readWaveform(j, wav_tmp);

            while(i<acquisitions && j<waveforms)
            {
                while(wav_tmp.head.time_stamp < acq_tmp.getHead().acquisition_time_stamp)
                {
                    write_waveform(stream, wav_tmp);

                    j++;

                    if(j<waveforms)
                    {
                        dataset->readWaveform(j, wav_tmp);
                    }
                    else
                    {
                        break;
                    }
                }

                while (acq_tmp.getHead().acquisition_time_stamp <= wav_tmp.head.time_stamp)
                {
                    write_acquisition(stream, acq_tmp);
                    i++;

                    if(i<acquisitions)
                    {
                        dataset->readAcquisition(i, acq_tmp);
                    }
                    else
                    {
                        break;
                    }
                }

                if(j==waveforms && i<acquisitions)
                {
                    write_acquisition(stream, acq_tmp);
                    for (uint32_t ia=i+1; ia<acquisitions; ia++)
                    {
                        dataset->readAcquisition(ia, acq_tmp);
                        write_acquisition(stream, acq_tmp);
                    }
                }

                if (i==acquisitions && j<waveforms)
                {
                    write_waveform(stream, wav_tmp);
                    for (uint32_t iw = j + 1; iw<waveforms; iw++)
                    {
                        dataset->readWaveform(iw, wav_tmp);
                        write_waveform(stream, wav_tmp);
                    }
                }
            }
        }
        else
        {
            for (i=0; i<acquisitions; i++)
            {
                dataset->readAcquisition(i, acq_tmp);
                write_acquisition(stream, acq_tmp);
            }
        }
    }

    void write_acquisition(std::ostream& stream, const ISMRMRD::Acquisition& acq)
    {
        const auto id = ::Gadgetron::Core::MessageID::GADGET_MESSAGE_ISMRMRD_ACQUISITION;
        stream.write(reinterpret_cast<const char*>(&id), sizeof(::Gadgetron::Core::MessageID));
        stream.write(reinterpret_cast<const char*>(&acq.getHead()), sizeof(ISMRMRD::AcquisitionHeader));

        unsigned long trajectory_elements = acq.getHead().trajectory_dimensions * acq.getHead().number_of_samples;
        unsigned long data_elements = acq.getHead().active_channels * acq.getHead().number_of_samples;

        if (trajectory_elements) {
             stream.write(reinterpret_cast<const char*>(&acq.getTrajPtr()[0]), sizeof(float) * trajectory_elements);
        }

        if (data_elements) {
            stream.write(reinterpret_cast<const char*>(&acq.getDataPtr()[0]), 2 * sizeof(float) * data_elements);
        }
    }

    void write_waveform(std::ostream& stream, const ISMRMRD::Waveform& wav)
    {
        const auto id = ::Gadgetron::Core::MessageID::GADGET_MESSAGE_ISMRMRD_WAVEFORM;
        stream.write(reinterpret_cast<const char*>(&id), sizeof(::Gadgetron::Core::MessageID));
        stream.write(reinterpret_cast<const char*>(&wav.head), sizeof(ISMRMRD::ISMRMRD_WaveformHeader));

        unsigned long data_elements = wav.head.channels*wav.head.number_of_samples;

        if (data_elements)
        {
            stream.write(reinterpret_cast<const char*>(wav.begin_data()), sizeof(uint32_t)*data_elements);
        }
    }

    void write_close(std::ostream& stream)
    {
        const auto id = ::Gadgetron::Core::MessageID::CLOSE;
        stream.write(reinterpret_cast<const char*>(&id), sizeof(::Gadgetron::Core::MessageID));
    }
};
