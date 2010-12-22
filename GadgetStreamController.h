#ifndef GADGETSTREAMCONTROLLER_H
#define GADGETSTREAMCONTROLLER_H

#include "ace/Log_Msg.h"
#include "ace/Reactor.h"
#include "ace/SOCK_Stream.h"
#include "ace/Stream.h"

#include "GadgetStreamConfigurator.h"

class GadgetStreamController : public ACE_Event_Handler
{
public:
  ACE_SOCK_Stream &peer (void) { return this->sock_; }


  int open (void);

  virtual ACE_HANDLE get_handle (void) const { 
    return this->sock_.get_handle (); 
  }

  virtual int handle_input (ACE_HANDLE fd = ACE_INVALID_HANDLE);
  virtual int handle_output (ACE_HANDLE fd = ACE_INVALID_HANDLE);
  virtual int handle_close (ACE_HANDLE handle,
                            ACE_Reactor_Mask close_mask);

  virtual int read_configuration();

protected:
  ACE_SOCK_Stream sock_;
  ACE_Stream<ACE_MT_SYNCH> stream_;
  bool stream_configured_;
  
};

#endif //GADGETSTREAMCONTROLLER_H
