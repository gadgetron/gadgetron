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
#include "writers/ImageWriter.h"
#include "readers/AcquisitionReader.h"
#include "readers/WaveformReader.h"
#include "readers/ImageReader.h"
#include "Reader.h"


class ErrorThrower : public ::Gadgetron::Server::Connection::ErrorReporter
{
  public:
    void operator()(const std::string& location, const std::string& message) override {
        throw std::runtime_error(("[" + location + "] ERROR: " + message));
    }
};


namespace {
const std::string DEFAULT_WORK_DIR = "/tmp/";

std::string get_gadgetron_home_dir()
{
    auto home_dir = std::getenv("GADGETRON_HOME");

    if (!home_dir) {
        home_dir = std::getenv("CONDA_PREFIX");
        if (!home_dir) {
            throw std::runtime_error(
                "Failed to find gadgetron home directory - ensure ${GADGETRON_HOME} or ${CONDA_PREFIX} is set.");
        }
    }

    return home_dir;
}

std::filesystem::path find_config_path(const std::string& config_xml)
{
    std::filesystem::path path;
    auto home_dir = get_gadgetron_home_dir();
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
        : args_(args), storage_address_(storage_address)
        , input_channel_(::Gadgetron::Core::make_channel<::Gadgetron::Core::MessageChannel>())
        , output_channel_(::Gadgetron::Core::make_channel<::Gadgetron::Core::MessageChannel>())
        , error_producer_(ErrorThrower()), error_handler_(error_producer_, std::string(__FILE__)) {}
    ~StreamConsumer() {}

    void consume(std::istream& input_stream, std::ostream& output_stream, std::string config_xml_name)
    {
        ::Gadgetron::Core::MessageID id = ::Gadgetron::Core::MessageID::ERROR;
        input_stream.read(reinterpret_cast<char*>(&id), sizeof(::Gadgetron::Core::MessageID));
        output_stream.write(reinterpret_cast<char*>(&id), sizeof(::Gadgetron::Core::MessageID));

        ISMRMRD::IsmrmrdHeader hdr;
        switch(id)
        {
            case ::Gadgetron::Core::MessageID::HEADER:
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
                throw std::runtime_error("Expected HEADER enumeration: " + std::to_string(id));
            }
        }

        boost::filesystem::path home_dir(get_gadgetron_home_dir());
        boost::filesystem::path working_dir(DEFAULT_WORK_DIR);

        ::Gadgetron::Core::Context::Paths paths{home_dir, working_dir};

        auto storage_spaces = ::Gadgetron::Server::setup_storage_spaces(storage_address_, hdr);
        context_ =
            std::make_unique<::Gadgetron::Core::StreamContext>(hdr, paths, args_, storage_address_, storage_spaces);
        loader_ = std::make_unique<::Gadgetron::Server::Connection::Loader>(*context_);

        auto config_path = find_config_path(config_xml_name);

        std::ifstream file(config_path, std::ios::in | std::ios::binary);
        if (!file.is_open()) {
            throw std::runtime_error("Failed to open file at path: " + config_path.string());
        }

        config_ = ::Gadgetron::Server::Connection::parse_config(file);
        file.close();

        stream_ = loader_->load(config_.stream);

        processing_ = true;
        process_future_ = std::async(std::launch::async, [&]() {
            try
            {
                stream_->process(std::move(input_channel_.input),
                                std::move(output_channel_.output),
                                error_handler_);
                processing_ = false;
            }
            catch (const std::exception& exc)
            {
                { // Induce a ::Gadgetron::Core::ChannelClosed exception upon readers of the channel.
                    auto destruct_me = std::move(output_channel_.output);
                }
                processing_ = false;
                throw;
            }
        });

        // Begin writing acquisitions
        auto input_future = std::async(std::launch::async, [&]() {
            try
            {
                input_stream.exceptions(std::istream::failbit | std::istream::badbit | std::istream::eofbit);
                auto acq_reader = ::Gadgetron::Core::Readers::AcquisitionReader();
                auto wav_reader = ::Gadgetron::Core::Readers::WaveformReader();
                auto img_reader = ::Gadgetron::Core::Readers::ImageReader();
                bool closed = false;

                while (!input_stream.eof() && !closed)
                {
                    ::Gadgetron::Core::MessageID id = ::Gadgetron::Core::MessageID::ERROR;
                    input_stream.read(reinterpret_cast<char*>(&id), sizeof(::Gadgetron::Core::MessageID));

                    switch(id)
                    {
                        case ::Gadgetron::Core::MessageID::GADGET_MESSAGE_ISMRMRD_ACQUISITION:
                        {
                            input_channel_.output.push_message(acq_reader.read(input_stream));
                            break;
                        }
                        case ::Gadgetron::Core::MessageID::GADGET_MESSAGE_ISMRMRD_WAVEFORM:
                        {
                            input_channel_.output.push_message(wav_reader.read(input_stream));
                            break;
                        }
                        case ::Gadgetron::Core::MessageID::GADGET_MESSAGE_ISMRMRD_IMAGE:
                        {
                            input_channel_.output.push_message(img_reader.read(input_stream));
                            break;
                        }
                        case ::Gadgetron::Core::MessageID::ERROR:
                        case ::Gadgetron::Core::MessageID::CLOSE:
                        {
                            auto destruct_me = std::move(input_channel_.output);
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
            catch(std::ios_base::failure& exc)
            {
                // Induce a ::Gadgetron::Core::ChannelClosed exception upon readers of the channel.
                auto destruct_me = std::move(input_channel_.output);
            }
        });

        auto output_future = std::async(std::launch::async, [&]()
        {
            auto writer = ::Gadgetron::Core::Writers::ImageWriter();
            while (true)
            {
                try
                {
                    auto message = output_channel_.input.pop();
                    writer.write(output_stream, std::move(message));
                }
                catch (const ::Gadgetron::Core::ChannelClosed& exc)
                {
                    break;
                }
            }

            ::Gadgetron::Core::MessageID close_id = ::Gadgetron::Core::MessageID::CLOSE;
            output_stream.write(reinterpret_cast<char*>(&close_id), sizeof(::Gadgetron::Core::MessageID));
        });

        // Clean up and propagate exceptions if they occurred
        input_future.get();
        output_future.get();
        process_future_.get();
    }

  private:
    boost::program_options::variables_map args_;
    std::string storage_address_;

    ::Gadgetron::Core::ChannelPair input_channel_;
    ::Gadgetron::Core::ChannelPair output_channel_;
    ErrorThrower error_producer_;
    ::Gadgetron::Server::Connection::ErrorHandler error_handler_;
    std::unique_ptr<::Gadgetron::Core::StreamContext> context_;
    std::unique_ptr<::Gadgetron::Server::Connection::Loader> loader_;
    std::unique_ptr<::Gadgetron::Server::Connection::Nodes::Stream> stream_;
    ::Gadgetron::Server::Connection::Config config_;

    std::future<void> process_future_;
    std::atomic<bool> processing_ = false;
};
