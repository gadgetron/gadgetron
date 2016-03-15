#ifndef GADGETRON_CLOUDBUS_H
#define GADGETRON_CLOUDBUS_H

#include "cloudbus_export.h"
#include <ace/Task.h>
#include <ace/INET_Addr.h>
#include <ace/OS_NS_unistd.h>
#include <ace/SOCK_Connector.h>
#include <ace/SOCK_Stream.h>
#include <ace/SOCK_Acceptor.h>
#include <ace/Svc_Handler.h>
#include <ace/Condition_T.h>

#include "log.h"

#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

#include <iostream>
#include <map>
#include <utility>
#include <time.h>
#include <vector>

#define GADGETRON_DEFAULT_RELAY_ADDR "localhost"
#define GADGETRON_DEFAULT_RELAY_PORT 8002
#define MAXHOSTNAMELENGTH 1024

#include "cloudbus_io.h"

namespace Gadgetron
{  

  class CloudBus;
  
  class CloudBusReaderTask : public ACE_Task<ACE_MT_SYNCH>
  {
    
  public:
    typedef ACE_Task<ACE_MT_SYNCH> inherited;
    
  CloudBusReaderTask(CloudBus* cloudbus)
    : inherited()
    , cloud_bus_(cloudbus)
    {
      
    }
    
    virtual int open(void*);
    virtual int close(unsigned long flags);
    virtual int svc(void);

  protected:
    ACE_SOCK_Stream* socket_;
    CloudBus* cloud_bus_;
  };


  
  class EXPORTCLOUDBUS CloudBus : public ACE_Svc_Handler<ACE_SOCK_STREAM, ACE_MT_SYNCH>
  {
    friend CloudBusReaderTask;
    
  public:
    typedef ACE_Task<ACE_MT_SYNCH> inherited;

    virtual int open(void* = 0);
    virtual int handle_close (ACE_HANDLE handle, ACE_Reactor_Mask close_mask);

    virtual int handle_input (ACE_HANDLE fd = ACE_INVALID_HANDLE);
    virtual int svc(void);

    static CloudBus* instance();
    static void set_relay_address(const char* addr);
    static void set_relay_port(int port);
    static void set_query_only(bool m = true);
    static void set_gadgetron_port(uint32_t port);
    static void set_rest_port(uint32_t port);

    void set_lb_endpoint(std::string addr, uint32_t port);
    
    void set_compute_capability(uint32_t c);

    void send_node_info();    
    void update_node_info();
    void get_node_info(std::vector<GadgetronNodeInfo>& nodes);
    void print_nodes();
    size_t get_number_of_nodes();

    unsigned int active_reconstructions();
    unsigned int port();
    const char* uuid();
    
    void report_recon_start();
    void report_recon_end();
    
  protected:
    ///Protected constructor. 
    CloudBus(int port, const char* addr);
    CloudBus(); 

    static CloudBus* instance_;
    static const char* relay_inet_addr_;
    static int relay_port_;
    static bool query_mode_; //Listen only
    static int gadgetron_port_;
    static int rest_port_;

    bool use_lb_endpoint_;
    std::string lb_address_;
    uint32_t lb_port_;

    GadgetronNodeInfo node_info_;
    std::vector<GadgetronNodeInfo> nodes_;
    
    ACE_Thread_Mutex mtx_;
    ACE_Thread_Mutex mtx_node_list_;
    ACE_Condition<ACE_Thread_Mutex> node_list_condition_;
    
    boost::uuids::uuid uuid_;
    bool connected_;
    ACE_SOCK_Stream socket_;

    CloudBusReaderTask* reader_task_;
  };
}

#endif
