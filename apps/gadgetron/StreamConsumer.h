#include <filesystem>
#include <future>
#include <memory>
#include <sstream>
#include <thread>

#include <boost/asio.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/program_options/variables_map.hpp>

#include "Channel.h"
#include "connection/Core.h"
#include "connection/Loader.h"

#include "MessageID.h"
#include "Reader.h"
#include "readers/AcquisitionBucketReader.h"
#include "readers/AcquisitionReader.h"
#include "readers/ImageReader.h"
#include "readers/IsmrmrdImageArrayReader.h"
#include "readers/WaveformReader.h"
#include "readers/TextReader.h"
#include "writers/AcquisitionBucketWriter.h"
#include "writers/ImageWriter.h"
#include "writers/IsmrmrdImageArrayWriter.h"
#include "writers/TextWriter.h"

using namespace Gadgetron::Core;
using namespace Gadgetron::Server;

namespace
{

class ErrorThrower : public Connection::ErrorReporter
{
  public:
    void operator()(const std::string& location, const std::string& message) override {
        throw std::runtime_error(("[" + location + "] ERROR: " + message));
    }
};

std::filesystem::path find_config_path(const std::string& home_dir, const std::string& config_xml)
{
    auto config_path = std::filesystem::path(home_dir) / std::filesystem::path("share/gadgetron/config") /
                       std::filesystem::path(config_xml);

    if (!std::filesystem::is_regular_file(config_path)) {
        throw std::runtime_error("Failed to find gadgetron configuration at the expected path: " +
                                 config_path.string());
    }

    return config_path;
}

} // namespace

class StreamConsumer
{
public:
    StreamConsumer(const boost::program_options::variables_map& args, std::string storage_address)
        : args_(args), storage_address_(storage_address) {}
    ~StreamConsumer() {}

    void consume(std::istream& input_stream, std::ostream& output_stream, std::string config_xml_name)
    {
        Context::Paths paths{
            args_["home"].as<boost::filesystem::path>().string(),
            args_["dir"].as<boost::filesystem::path>().string()};

        ISMRMRD::IsmrmrdHeader hdr = consume_ismrmrd_header(input_stream, output_stream);
        auto storage_spaces = setup_storage_spaces(storage_address_, hdr);

        auto context = StreamContext(hdr, paths, args_, storage_address_, storage_spaces);
        auto loader = Connection::Loader(context);
        auto config_path = find_config_path(args_["home"].as<boost::filesystem::path>().string(), config_xml_name);

        std::ifstream file(config_path, std::ios::in | std::ios::binary);
        if (!file.is_open()) {
            throw std::runtime_error("Failed to open file at path: " + config_path.string());
        }

        auto config = Connection::parse_config(file);
        file.close();

        auto stream = loader.load(config.stream);
        auto input_channel = make_channel<MessageChannel>();
        auto output_channel = make_channel<MessageChannel>();
        std::atomic<bool> processing = true;

        auto process_future = std::async(std::launch::async, [&]() {
            try
            {
                ErrorThrower error_thrower;
                Connection::ErrorHandler error_handler(error_thrower, std::string(__FILE__));
                stream->process(std::move(input_channel.input), std::move(output_channel.output), error_handler);
                processing = false;
            }
            catch (const std::exception& exc)
            {
                {
                    // Induce a ChannelClosed exception upon readers of the channel.
                    auto destruct_me = std::move(output_channel.output);
                }
                processing = false;
                throw;
            }
        });

        auto input_future = std::async(std::launch::async, [&]() {
            try
            {
                consume_input_messages(input_stream, input_channel);
            }
            catch(std::ios_base::failure& exc)
            {
                // Induce a ChannelClosed exception upon readers of the channel.
                auto destruct_me = std::move(input_channel.output);
            }
        });

        auto output_future = std::async(std::launch::async, [&]()
        {
            process_output_messages(output_channel, output_stream);
        });

        // Clean up and propagate exceptions if they occurred
        input_future.get();
        output_future.get();
        process_future.get();
    }

  private:

    ISMRMRD::IsmrmrdHeader consume_ismrmrd_header(std::istream& input_stream, std::ostream& output_stream)
    {
        ISMRMRD::IsmrmrdHeader hdr;
        MessageID id = MessageID::CLOSE;
        input_stream.read(reinterpret_cast<char*>(&id), sizeof(MessageID));
        output_stream.write(reinterpret_cast<char*>(&id), sizeof(MessageID));

        switch(id)
        {
            case MessageID::HEADER:
            {
                uint32_t hdr_size = 0;
                input_stream.read(reinterpret_cast<char*>(&hdr_size), sizeof(uint32_t));
                output_stream.write(reinterpret_cast<char*>(&hdr_size), sizeof(uint32_t));

                if(hdr_size > 0)
                {
                    std::vector<char> data(hdr_size);

                    input_stream.read(data.data(), hdr_size);
                    output_stream.write(data.data(), hdr_size);
                    ISMRMRD::deserialize(std::string(data.data(), data.size()).c_str(), hdr);
                }
                else
                {
                    throw std::runtime_error("Expected size > 0, got: " + std::to_string(hdr_size));
                }

                break;
            }
            default:
            {
                throw std::runtime_error("Expected HEADER enumeration, got: " + std::to_string(id));
            }
        }

        return hdr;
    }

    void consume_input_messages(std::istream& input_stream, ChannelPair& input_channel)
    {
        input_stream.exceptions(std::istream::failbit | std::istream::badbit | std::istream::eofbit);

        auto acq_reader = Readers::AcquisitionReader();
        auto wav_reader = Readers::WaveformReader();
        auto img_reader = Readers::ImageReader();
        auto img_array_reader = Readers::IsmrmrdImageArrayReader();
        auto acq_bucket_reader = Readers::AcquisitionBucketReader();
        auto text_reader = Readers::TextReader();
        bool closed = false;

        while (!input_stream.eof() && !closed)
        {
            MessageID id = MessageID::ERROR;
            input_stream.read(reinterpret_cast<char*>(&id), sizeof(MessageID));

            switch(id)
            {
                case MessageID::GADGET_MESSAGE_ISMRMRD_ACQUISITION:
                {
                    input_channel.output.push_message(acq_reader.read(input_stream));
                    break;
                }
                case MessageID::GADGET_MESSAGE_ISMRMRD_WAVEFORM:
                {
                    input_channel.output.push_message(wav_reader.read(input_stream));
                    break;
                }
                case MessageID::GADGET_MESSAGE_ISMRMRD_IMAGE:
                {
                    input_channel.output.push_message(img_reader.read(input_stream));
                    break;
                }
                case MessageID::GADGET_MESSAGE_ISMRMRD_IMAGE_ARRAY:
                {
                    input_channel.output.push_message(img_array_reader.read(input_stream));
                    break;
                }
                case MessageID::GADGET_MESSAGE_BUCKET:
                {
                    input_channel.output.push_message(acq_bucket_reader.read(input_stream));
                    break;
                }
                case MessageID::TEXT:
                {
                    //Gadgetron::Core::Message msg = text_reader.read(input_stream);
                    // if (convertible_to<std::string>(msg))
                    // {
                    //     std::string str = Gadgetron::Core::force_unpack<std::string>(std::move(msg));
                    //     GDEBUG_STREAM("Receive text message : " << str);
                    // }
                    input_channel.output.push_message(text_reader.read(input_stream));
                    break;
                }
                case MessageID::ERROR:
                {
                    throw std::runtime_error("Got error while processing input stream");
                }
                case MessageID::CLOSE:
                {
                    auto destruct_me = std::move(input_channel.output);
                    closed = true;
                    break;
                }
                default:
                {
                    throw std::runtime_error("Unsupported message ID: " + std::to_string(id));
                }
            }
        }
    }

    void process_output_messages(ChannelPair& output_channel, std::ostream& output_stream)
    {
        auto writer = Writers::ImageWriter();
        auto img_array_writer = Writers::IsmrmrdImageArrayWriter();
        auto acq_bucket_writer = Writers::AcquisitionBucketWriter();
        auto text_writer = Writers::TextWriter();

        while (true)
        {
            try
            {
                auto message = output_channel.input.pop();

                if (convertible_to<Gadgetron::AcquisitionBucket>(message) )
                {
                    acq_bucket_writer.write(output_stream, std::move(message));
                }
                else if (convertible_to<Gadgetron::IsmrmrdImageArray>(message) )
                {
                    img_array_writer.write(output_stream, std::move(message));
                }
                else if (convertible_to<std::string>(message) )
                {
                    text_writer.write(output_stream, std::move(message));
                }
                else
                {
                    writer.write(output_stream, std::move(message));
                }
            }
            catch (const ChannelClosed& exc)
            {
                break;
            }
        }

        MessageID close_id = MessageID::CLOSE;
        output_stream.write(reinterpret_cast<char*>(&close_id), sizeof(MessageID));
    }

    boost::program_options::variables_map args_;
    std::string storage_address_;
};
