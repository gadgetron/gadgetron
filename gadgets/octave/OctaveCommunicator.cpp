#include "OctaveCommunicator.h"

#include <octave/oct.h>
#include <octave/octave.h>
#include <octave/parse.h>

#include <iostream>

OctaveCommunicator* OctaveCommunicator::instance()
{
  if (!instance_) instance_ = new OctaveCommunicator();
  return instance_;
}

OctaveCommunicator::OctaveCommunicator()
{
  const char * argvv [] = {"" /* name of program, not relevant */, "--silent"}; 
  octave_main (2, (char **) argvv, true /* embedded */);
}

OctaveCommunicator::~OctaveCommunicator()
{
}

void OctaveCommunicator::register_gadget(Gadget* g)
{
  std::cout << "Registering Gadget with Communicator = " << this << std::endl;
  gadget_map_[g->module()->name()] = g;
}

bool OctaveCommunicator::message_gadget(std::string gadget)
{
  std::cout << "Messaging Gadget with Communicator = " << this << std::endl;
  std::map<std::string, Gadget*>::iterator it = gadget_map_.find(gadget);
  if (it != gadget_map_.end()) {
    std::cout << "Messaging Gadget with ID = " << gadget << ", message is: " << it->second << std::endl;
    return true;
  } else {
    std::cout << "Gadget with ID = " << gadget << " NOT FOUND!" << std::endl;
    return false;
  }
  return false;
}

OctaveCommunicator* OctaveCommunicator::instance_ = NULL;
/*
bool message_gadget(std::string gadget)
{
  return OctaveCommunicator::instance()->message_gadget(gadget);
}
*/
