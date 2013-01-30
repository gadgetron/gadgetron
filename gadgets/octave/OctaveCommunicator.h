#ifndef OCTAVECOMMUNICATOR_H
#define OCTAVECOMMUNICATOR_H

#include "gadgetronoctavecommunicator_export.h"

#include "Gadget.h"

#include <map>
#include <string>

class EXPORTGADGETSOCTAVECOMMUNICATOR OctaveCommunicator
{

 public:
  static OctaveCommunicator* instance(); 
  
  void register_gadget(Gadget* g);

  bool message_gadget(std::string g, ACE_Message_Block* m);

 private:
  OctaveCommunicator();
  ~OctaveCommunicator();
  
  static OctaveCommunicator* instance_;
  
  std::map<std::string, Gadget*> gadget_map_;
};


#endif
