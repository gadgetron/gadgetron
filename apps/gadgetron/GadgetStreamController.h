#ifndef GADGETSTREAMCONTROLLER_H
#define GADGETSTREAMCONTROLLER_H

#include "ace/Log_Msg.h"
#include "ace/Reactor.h"
#include "ace/SOCK_Stream.h"
#include "ace/Message_Queue.h"
#include "ace/Svc_Handler.h"
#include "ace/Reactor_Notification_Strategy.h"

#include <complex>
#include <vector>

#include "gadgetbase_export.h"
#include "GadgetronConnector.h"
#include "GadgetStreamInterface.h"


namespace Gadgetron{

class EXPORTGADGETBASE GadgetStreamController 
  : public ACE_Svc_Handler<ACE_SOCK_STREAM, ACE_MT_SYNCH>
  , public GadgetStreamInterface
{
public:
  GadgetStreamController();

  virtual ~GadgetStreamController();

  virtual int open (void);
  virtual int svc(void);


  virtual int handle_input (ACE_HANDLE fd = ACE_INVALID_HANDLE);
  virtual int handle_close (ACE_HANDLE handle,
                            ACE_Reactor_Mask close_mask);

  virtual int output_ready(ACE_Message_Block* mb);

private:
  WriterTask writer_task_;
  ACE_Reactor_Notification_Strategy notifier_;
  GadgetMessageReaderContainer readers_;
  virtual int configure(const std::string&config_file_stream);
  virtual int configure_from_file(std::string filename);
};

}
#endif //GADGETSTREAMCONTROLLER_H
