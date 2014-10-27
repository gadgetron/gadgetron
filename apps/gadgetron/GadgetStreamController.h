#ifndef GADGETSTREAMCONTROLLER_H
#define GADGETSTREAMCONTROLLER_H

#include "ace/Log_Msg.h"
#include "ace/Reactor.h"
#include "ace/SOCK_Stream.h"
#include "ace/Stream.h"
#include "ace/Message_Queue.h"
#include "ace/Svc_Handler.h"
#include "ace/Reactor_Notification_Strategy.h"

#include <complex>
#include <vector>

#include "Gadgetron.h"
//#include "Gadget.h"
#include "gadgetron_paths.h"
#include "GadgetronConnector.h"

typedef ACE_Module<ACE_MT_SYNCH> GadgetModule;

namespace Gadgetron{

class GadgetStreamController 
: public ACE_Svc_Handler<ACE_SOCK_STREAM, ACE_MT_SYNCH>
{
public:
  GadgetStreamController();

  virtual ~GadgetStreamController()
    { 
      //ACE_DEBUG( (LM_INFO, ACE_TEXT("~GadgetStreamController() called\n")) );
    }

  //ACE_SOCK_Stream &peer (void) { return this->sock_; }

  int open (void);

  /*
  virtual ACE_HANDLE get_handle (void) const { 
    return this->sock_.get_handle (); 
  }
  */

  virtual int handle_input (ACE_HANDLE fd = ACE_INVALID_HANDLE);
  //virtual int handle_output (ACE_HANDLE fd = ACE_INVALID_HANDLE);
  virtual int handle_close (ACE_HANDLE handle,
                            ACE_Reactor_Mask close_mask);

  virtual int output_ready(ACE_Message_Block* mb);

  virtual Gadget* find_gadget(std::string gadget_name);

  void set_global_gadget_parameters(const std::map<std::string, std::string>& globalGadgetPara);

private:
  ACE_Stream<ACE_MT_SYNCH> stream_;
  bool stream_configured_;
  WriterTask writer_task_;

  ACE_Reactor_Notification_Strategy notifier_;

  GadgetMessageReaderContainer readers_;
  
  std::vector<ACE_DLL_Handle*> dll_handles_;

  std::map<std::string, std::string> global_gadget_parameters_;

  virtual int configure(std::string config_xml_string);
  virtual int configure_from_file(std::string config_xml_filename);

  virtual GadgetModule * create_gadget_module(const char* DLL, const char* gadget, const char* gadget_module_name);

  template <class T>  T* load_dll_component(const char* DLL, const char* component_name);
  
  std::string gadgetron_home_;

};
}
#endif //GADGETSTREAMCONTROLLER_H
