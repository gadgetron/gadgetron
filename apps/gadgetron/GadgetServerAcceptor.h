#ifndef _GADGETSERVERACCEPTOR_H
#define _GADGETSERVERACCEPTOR_H

#include "ace/SOCK_Acceptor.h"
#include "ace/Reactor.h"
#include <string>
#include <map>

namespace Gadgetron{
class GadgetServerAcceptor : public ACE_Event_Handler
{
public:
  virtual ~GadgetServerAcceptor ();

  int open (const ACE_INET_Addr &listen_addr);

  virtual ACE_HANDLE get_handle (void) const
    { return this->acceptor_.get_handle (); }

  virtual int handle_input (ACE_HANDLE fd = ACE_INVALID_HANDLE);

  virtual int handle_close (ACE_HANDLE handle,
                            ACE_Reactor_Mask close_mask);

  std::map<std::string, std::string> global_gadget_parameters_;

protected:
  ACE_SOCK_Acceptor acceptor_;
};
}
#endif //_GADGETSERVERACCEPTOR_H
