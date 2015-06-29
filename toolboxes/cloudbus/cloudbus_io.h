#ifndef CLOUDBUS_IO_H
#define CLOUDBUS_IO_H

#include <exception>
#include <cstring>
#include <vector>

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
  
  struct GadgetronNodeInfo
  {
    std::string uuid;
    std::string address;
    uint32_t port;
    uint32_t compute_capability;
  };

  size_t calculate_node_info_length(GadgetronNodeInfo& n)
  {
    size_t len = 0;
    len += 4 + n.uuid.size();
    len += 4 + n.address.size();
    len += 8;
    return len;
  }
  
  size_t serialize(GadgetronNodeInfo& n, char* buffer, size_t buf_len)
  {
    size_t pos = 0;
    
    if (buf_len < calculate_node_info_length(n)) {
      throw std::runtime_error("Provided buffer is too short for serialization");
    }
    
    *((uint32_t*)(buffer + pos)) = n.uuid.size(); pos += 4;
    memcpy ((buffer+pos), n.uuid.c_str(), n.uuid.size() ); pos += n.uuid.size();

    *((uint32_t*)(buffer + pos)) = n.address.size(); pos += 4;
    memcpy ((buffer+pos), n.address.c_str(), n.address.size() ); pos += n.address.size();

    *((uint32_t*)(buffer + pos)) = n.port; pos += 4;
    *((uint32_t*)(buffer + pos)) = n.compute_capability; pos += 4;
    
    return pos;
  }

  size_t deserialize(GadgetronNodeInfo& n, char* buffer, size_t buf_len)
  {
    size_t pos = 0;
    
    if (buf_len < 16) throw std::runtime_error("Provided buffer is too small to hold node info");

    size_t uuid_size = *((uint32_t*)(buffer+pos)); pos += 4;
    n.uuid = std::string(buffer+pos,uuid_size); pos += uuid_size;

    size_t address_size = *((uint32_t*)(buffer+pos)); pos += 4;
    n.address = std::string(buffer+pos,address_size); pos += address_size;

    n.port = *((uint32_t*)(buffer+pos)); pos += 4;
    n.compute_capability = *((uint32_t*)(buffer+pos)); pos += 4;

    return pos;
  }


  size_t calculate_node_info_list_length(std::vector<GadgetronNodeInfo>& nl)
  {
    size_t length = 0;
    for (std::vector<GadgetronNodeInfo>::iterator it = nl.begin();
	 it != nl.end(); it++)
      {
	length += calculate_node_info_length(*it);
      }
    return length;
  }
  
  size_t serialize(std::vector<GadgetronNodeInfo>& nl, char* buffer, size_t buf_len)
  {
    size_t serialized_length = calculate_node_info_list_length(nl);
    if ((serialized_length+4) > buf_len) throw std::runtime_error("Buffer too short for serializing node info list");

    *((uint32_t*)buffer) = nl.size();

    size_t pos = 4;

    for (std::vector<GadgetronNodeInfo>::iterator it = nl.begin();
	 it != nl.end(); it++)
      {
	pos += serialize(*it,buffer+pos,buf_len-pos);
      }

    return pos;
  }


  size_t deserialize(std::vector<GadgetronNodeInfo>& nl, char* buffer, size_t buf_len)
  {
    nl.clear();
    size_t pos = 0;

    uint32_t num_nodes = *((uint32_t*)buffer);
    pos += 4;
    
    for (unsigned int i = 0; i < num_nodes; i++) {
      GadgetronNodeInfo n;
      pos += deserialize(n,buffer+pos,buf_len-pos);
      nl.push_back(n);
    }
    
    return pos;
  }

}

#endif
