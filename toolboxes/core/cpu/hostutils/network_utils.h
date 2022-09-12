#pragma once
#include <string>
#include <vector>

namespace Gadgetron {

    /// find the ip of current gadgetron server
    void find_gadgetron_ip(std::string& host_name, std::vector<std::string>& ip_list);

    struct IP_list{
        std::string host_name;
        std::vector<std::string> ip_list;
    };

    IP_list find_gadgetron_ip();
}
