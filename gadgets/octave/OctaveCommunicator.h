#ifndef OCTAVECOMMUNICATOR_H
#define OCTAVECOMMUNICATOR_H

#include <ace/Synch.h>
#include <ace/Mutex.h>

#include <octave/oct.h>
#include <octave/octave.h>
#include <octave/parse.h>


#include "Gadget.h"

#include <map>
#include <string>

#include "gadgetronoctavecommunicator_export.h"

class EXPORTGADGETSOCTAVECOMMUNICATOR OctaveCommunicator
{

 public:
  static OctaveCommunicator* instance(); 
  
  void register_gadget(Gadgetron::Gadget* g);
  bool message_gadget(std::string g, ACE_Message_Block* m);
  octave_value_list octave_feval (const std::string &name, const octave_value_list &args=octave_value_list(), int nargout=0);

 private:
  OctaveCommunicator();
  ~OctaveCommunicator();
  
  static OctaveCommunicator* instance_;
  ACE_Thread_Mutex mutex_;
  
  std::map<std::string, Gadgetron::Gadget*> gadget_map_;
};


#endif
