#ifndef CLOUDBUS_IO_H
#define CLOUDBUS_IO_H

#include <exception>
#include <cstring>
#include <vector>
#include <chrono>
#include <ctime>

#include "cloudbus_export.h"

namespace Gadgetron
{
  enum CLOUDBUS_MESSAGE_TYPE
  {
    GADGETRON_CLOUDBUS_MESSAGE_MIN = 0,
    GADGETRON_CLOUDBUS_NODE_INFO = 1,
    GADGETRON_CLOUDBUS_NODE_LIST_QUERY = 2,
    GADGETRON_CLOUDBUS_NODE_LIST_REPLY = 3,
    GADGETRON_CLOUDBUS_MESSAGE_MAX
  };
  
  struct EXPORTCLOUDBUS GadgetronNodeInfo
  {
    std::string uuid;
    std::string address;
    uint32_t port;
    uint32_t rest_port;
    uint32_t compute_capability;
    uint32_t active_reconstructions;
    std::time_t last_recon;
  };

  EXPORTCLOUDBUS size_t calculate_node_info_length(GadgetronNodeInfo& n);

  EXPORTCLOUDBUS size_t serialize(GadgetronNodeInfo& n, char* buffer, size_t buf_len);
  EXPORTCLOUDBUS size_t deserialize(GadgetronNodeInfo& n, char* buffer, size_t buf_len);

  EXPORTCLOUDBUS size_t calculate_node_info_list_length(std::vector<GadgetronNodeInfo>& nl);

  EXPORTCLOUDBUS size_t serialize(std::vector<GadgetronNodeInfo>& nl, char* buffer, size_t buf_len);
  EXPORTCLOUDBUS size_t deserialize(std::vector<GadgetronNodeInfo>& nl, char* buffer, size_t buf_len);
}

#endif
