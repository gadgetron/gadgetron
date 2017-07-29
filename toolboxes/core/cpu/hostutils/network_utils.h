#pragma once
#include "hostutils_export.h"
#include <string>
#include <vector>

namespace Gadgetron {

    /// find the ip of current gadgetron server
    EXPORTHOSTUTILS void find_gadgetron_ip(std::string& host_name, std::vector<std::string>& ip_list);
}
