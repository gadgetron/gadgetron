
#include "network_utils.h"
#include <boost/asio.hpp>

#include <iostream>
#include <stack>
#include <cmath>
#include <cstdlib>
#include "log.h"
#include "Process.h"
#include <regex>

namespace Gadgetron {

    void find_gadgetron_ip(std::string& host_name, std::vector<std::string>& ip_list)
    {
        try
        
        {
            boost::asio::io_service io_service;
            boost::asio::ip::tcp::resolver resolver(io_service);

            ip_list.clear();

            host_name = boost::asio::ip::host_name();

#ifdef WIN32
            boost::asio::ip::tcp::resolver::iterator iter = resolver.resolve({ host_name, "" });
            boost::asio::ip::tcp::resolver::iterator end;

            while (iter != end)
            {
                ip_list.push_back(iter->endpoint().address().to_string());
                iter++;
            }
#else

            namespace bp = boost::process;

            bp::ipstream stream;
            auto c = Process::child(
                    bp::search_path("ifconfig"),
                    bp::std_out > stream
            );

            auto reg = std::regex(R"(inet\s+(\S+))");
            std::string line;
            const std::string inet = "inet ";
            while(c.running() && std::getline(stream, line) && !line.empty()){
                std::smatch s;
                if (!std::regex_search(line,s,reg)) continue;

                ip_list.emplace_back(s[1]);

            }
#endif // WIN32
        }
        catch (...)
        {
            GADGET_THROW("Errors in find_gadgetron_ip(std::string& host_name, std::vector<std::string>& ip_list) ...");
        }
    }


    IP_list find_gadgetron_ip(){
        auto ip_list = IP_list{};

        try {
            find_gadgetron_ip(ip_list.host_name,ip_list.ip_list);
        } 
        catch (...){
            GERROR("Errors in find_gadgetron_ip() \n");
        }

        return ip_list;



    }
}
