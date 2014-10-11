#ifndef GADGETRON_CLOUDBUS_H
#define GADGETRON_CLOUDBUS_H

#include "cloudbus_export.h"
#include <ace/Task.h>
#include <ace/INET_Addr.h>
#include <ace/SOCK_Dgram_Mcast.h>
#include <ace/OS_NS_unistd.h>
#include <ace/SOCK_Acceptor.h>

#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

#include <iostream>
#include <map>
#include <utility>
#include <time.h>
#include <vector>

#define GADGETRON_DEFAULT_MULTICAST_ADDR "224.9.9.2"
#define GADGETRON_DEFAULT_MULTICAST_PORT 4148
#define GADGETRON_NODE_INFO_MESSAGE_LENGTH 16+sizeof(uint32_t)*2 //16 bytes for uuid + 2 ints

namespace Gadgetron
{

  struct GadgetronNodeInfo
  {
    std::string uuid;
    std::string address;
    uint32_t port;
    uint32_t compute_capability;
  };
  
  class CloudBusTask : public ACE_Task<ACE_MT_SYNCH>
  {
  public:
    typedef ACE_Task<ACE_MT_SYNCH> inherited;    
    CloudBusTask(int port, const char* addr);
    CloudBusTask();
    virtual int open(void* = 0);

  protected:
    ACE_SOCK_Dgram_Mcast mcast_dgram_;
    ACE_INET_Addr mcast_addr_;
  };

  class CloudBusReceiverTask : public CloudBusTask
  {
  public:
    CloudBusReceiverTask(int port, const char* addr);
    virtual int open(void* = 0);
    virtual int close(u_long flags = 0);
    virtual int svc(void);
  };


  class CloudBusSenderTask : public CloudBusTask
  {
  public:
    CloudBusSenderTask(int port, const char* addr);
    virtual int open(void* = 0);
    virtual int svc(void);
  };

  class EXPORTCLOUDBUS CloudBus
  {
    friend class CloudBusReceiverTask;
    friend class CloudBusSenderTask;

    typedef std::map<std::string, std::pair<GadgetronNodeInfo, time_t> > map_type_;

  public:
    static CloudBus* instance();
    static void set_mcast_address(const char* addr);
    static void set_mcast_port(int port);
    static void set_query_only(bool m = true);

    void set_gadgetron_port(uint32_t port)
    {
      node_info_.port = port;
    }
    
    void set_compute_capability(uint32_t c)
    {
      node_info_.compute_capability = c;
    }

    void wait();

    void get_node_info(std::vector<GadgetronNodeInfo>& nodes);
    size_t get_number_of_nodes();

  protected:
    ///Protected constructor. 
    CloudBus(int port, const char* addr);

    void update_node(const char* a, GadgetronNodeInfo& info);
    void remove_stale_nodes();
    
    static CloudBus* instance_;
    static const char* mcast_inet_addr_;
    static int mcast_port_;
    static bool query_mode_; //Listen only

    GadgetronNodeInfo node_info_;
    map_type_ nodes_;
    
    CloudBusReceiverTask receiver_;
    CloudBusSenderTask   sender_;
    ACE_Thread_Mutex mtx_;

    boost::uuids::uuid uuid_;
  };


}

#endif
