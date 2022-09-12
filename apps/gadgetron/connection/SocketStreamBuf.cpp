//
// Created by dchansen on 3/25/19.
//

#include "SocketStreamBuf.h"

#include "Types.h"

#include <boost/asio.hpp>
namespace {
    using boost::asio::ip::tcp;

    std::unique_ptr<tcp::socket> connect_socket(
        const std::string& host, const std::string& service, boost::asio::io_service& context) {
        tcp::resolver resolver{ context };
        auto endpoint = *resolver.resolve(tcp::resolver::query(host, service));
        auto socket   = std::make_unique<tcp::socket>(context);
        socket->connect(endpoint);
        return std::move(socket);
    }

    class SocketStreamBuf : public std::streambuf {
    public:
        explicit SocketStreamBuf(std::unique_ptr<boost::asio::ip::tcp::socket> socket, size_t buffer_size = 1024);

    protected:
        std::streamsize xsputn(const char_type* data, std::streamsize length) override;

        int sync() override;
        int underflow() override;
        int overflow(int ch = traits_type::eof()) override;

    private:
        std::unique_ptr<boost::asio::ip::tcp::socket> socket;
        std::vector<char> input_buffer;
        std::vector<char> output_buffer;

        /* Other members */
    };

    int SocketStreamBuf::sync() {
        return this->overflow() != traits_type::eof() ? 0 : -1;
    }
    int SocketStreamBuf::underflow() {

        auto elements_read = socket->read_some(boost::asio::buffer(this->eback(), input_buffer.size()));

        this->setg(this->eback(), this->eback(), this->eback() + elements_read);
        return traits_type::to_int_type(*this->gptr());
    }
    int SocketStreamBuf::overflow(int ch) {
        if (this->pptr() != this->pbase()) {
            boost::asio::write(*socket, boost::asio::buffer(this->pbase(), std::distance(this->pbase(), this->pptr())));
            this->setp(this->pbase(), this->epptr());
        }
        if (ch != traits_type::eof()) {
            this->sputc(ch);
        }

        return 0;
    }

    std::streamsize SocketStreamBuf::xsputn(const char* data, std::streamsize length) {
        this->overflow();
        return boost::asio::write(*socket, boost::asio::buffer(data, length));
    }
    SocketStreamBuf::SocketStreamBuf(std::unique_ptr<boost::asio::ip::tcp::socket> socket, size_t buffer_size)
        : socket(std::move(socket)), input_buffer(buffer_size), output_buffer(buffer_size) {
        this->setg(input_buffer.data(), input_buffer.data() + buffer_size, input_buffer.data() + buffer_size);
        this->setp(output_buffer.data(), output_buffer.data() + buffer_size);
    }

    using namespace Gadgetron::Connection;
    class SocketStream : public std::iostream {
    public:
        explicit SocketStream(std::unique_ptr<boost::asio::ip::tcp::socket> socket)
            : std::iostream(new SocketStreamBuf(std::move(socket))) {
            buffer = std::unique_ptr<SocketStreamBuf>(static_cast<SocketStreamBuf*>(this->rdbuf()));
        }

        SocketStream(const std::string& host, const std::string& service,
            std::shared_ptr<boost::asio::io_service> io_service = std::make_shared<boost::asio::io_service>())
            : SocketStream(connect_socket(host, service, *io_service)) {
            this->io_service = io_service;
        }

        ~SocketStream() override = default;

    private:
        std::shared_ptr<boost::asio::io_service> io_service;
        std::unique_ptr<SocketStreamBuf> buffer;
    };
}


std::unique_ptr<std::iostream> Gadgetron::Connection::stream_from_socket(
    std::unique_ptr<boost::asio::ip::tcp::socket> socket) {
    return std::make_unique<SocketStream>(std::move(socket));
}

std::unique_ptr<std::iostream> Gadgetron::Connection::remote_stream(
    const std::string& host, const std::string& service) {
    return std::make_unique<SocketStream>(host, service);
}
