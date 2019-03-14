#ifndef GADGETRON_SERVER_H
#define GADGETRON_SERVER_H

#include <boost/asio.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/program_options/variables_map.hpp>

namespace Gadgetron::Server {

    class Server {

        using tcp = boost::asio::ip::tcp;

    public:
        Server(boost::asio::io_service &io_service, const boost::program_options::variables_map &args);

    private:
        void accept();
        void connection_handler(const boost::system::error_code &error);

        tcp::acceptor acceptor;
        std::unique_ptr<tcp::iostream> stream;
        const boost::program_options::variables_map &args;
    };
}




#endif //GADGETRON_SERVER_H
