#ifndef _GADGETSERVERACCEPTOR_H
#define _GADGETSERVERACCEPTOR_H

#include "ace/SOCK_Acceptor.h"
#include "ace/Reactor.h"
#include "gadgettools_export.h"
#include <string>

namespace Gadgetron{
class EXPORTGADGETTOOLS GadgetServerAcceptor : public ACE_Event_Handler
{
public:
  virtual ~GadgetServerAcceptor ();

  int open (const ACE_INET_Addr &listen_addr);

  virtual ACE_HANDLE get_handle (void) const
    { return this->acceptor_.get_handle (); }

  virtual int handle_input (ACE_HANDLE fd = ACE_INVALID_HANDLE);

  virtual int handle_close (ACE_HANDLE handle,
                            ACE_Reactor_Mask close_mask);

  std::string working_directory_;

protected:
  ACE_SOCK_Acceptor acceptor_;
};
}
#endif //_GADGETSERVERACCEPTOR_H
