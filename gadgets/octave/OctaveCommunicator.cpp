#include "OctaveCommunicator.h"


#include <iostream>

OctaveCommunicator* OctaveCommunicator::instance()
{
  if (!instance_) instance_ = new OctaveCommunicator();
  return instance_;
}

OctaveCommunicator::OctaveCommunicator()
  : mutex_("OctaveCommunicatorMutex")
{
  const char * argvv [] = {"" /* name of program, not relevant */, "--silent"}; 
  octave_main (2, (char **) argvv, true /* embedded */);
  octave_value_list in;
  octave_value_list out;

  const char* gadgetron_home = ACE_OS::getenv("GADGETRON_HOME");
  std::string path_name = std::string(gadgetron_home) + std::string("/octave");
  
  in = octave_value (path_name.c_str());
  out = feval ("addpath", in, 1);
  
}

OctaveCommunicator::~OctaveCommunicator()
{
}

void OctaveCommunicator::register_gadget(Gadget* g)
{
  mutex_.acquire();
  gadget_map_[g->module()->name()] = g;
  mutex_.release();
}

bool OctaveCommunicator::message_gadget(std::string gadget, ACE_Message_Block* m)
{
  std::map<std::string, Gadget*>::iterator it = gadget_map_.find(gadget);

  if (it != gadget_map_.end()) {
	  if (it->second->putq(m) < 0) {
		  return false;
	  } else {
		  return true;
	  }
  } else {
    std::cout << "Gadget with ID = " << gadget << " NOT FOUND!" << std::endl;
    m->release();
    return false;
  }
  return false;
}

octave_value_list OctaveCommunicator::octave_feval (const std::string &name, const octave_value_list &args, int nargout)
{
  mutex_.acquire();
  octave_value_list out = feval(name,args,nargout);
  mutex_.release();

  return out;
}

OctaveCommunicator* OctaveCommunicator::instance_ = NULL;
