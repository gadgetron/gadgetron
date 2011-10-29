#ifndef GADGET_H
#define GADGET_H
#pragma once

#include <ace/OS_NS_stdlib.h>
#include <ace/Task.h>
#include <ace/Stream.h>
#include <ace/Module.h>
#include <ace/OS_Memory.h>
#include <ace/Svc_Handler.h>
#include <ace/SOCK_Stream.h>

#include <map>
#include <boost/shared_ptr.hpp>

#include "GadgetContainerMessage.h"
#include "GadgetronExport.h"
#include "Gadgetron.h"

class GadgetStreamController;

class Gadget : public ACE_Task<ACE_MT_SYNCH>
{

public:
  typedef ACE_Task<ACE_MT_SYNCH> inherited;

  enum
  {
    GADGET_MESSAGE_CONFIG = (ACE_Message_Block::USER_FLAGS << 1)
  };

  Gadget() 
    : inherited()
    , desired_threads_(1)
    , pass_on_undesired_data_(false)
    , controller_(0)
  {
    ACE_TRACE(( ACE_TEXT("Gadget::Gadget") ));
  }

  virtual ~Gadget() 
  {
    GADGET_DEBUG2("Shutting down Gadget (%s)\n", this->module()->name());
  }
  

  virtual int init(void)
  {
    ACE_TRACE(( ACE_TEXT("Gadget::init") ));
    return 0;
  }

  virtual int open(void* = 0) 
  {
    ACE_TRACE(( ACE_TEXT("Gadget::open") ));
    
    return this->activate( THR_NEW_LWP | THR_JOINABLE,
			   this->desired_threads() );
  }

  int put(ACE_Message_Block *m, ACE_Time_Value* timeout = 0) 
  {
    ACE_TRACE(( ACE_TEXT("Gadget::put") ));

    return this->putq(m, timeout);
  }

  virtual unsigned int desired_threads()
  {
    ACE_TRACE(( ACE_TEXT("Gadget::desired_threads (get)") ));

    return desired_threads_;
  }

  virtual void desired_threads(unsigned int t)
  {
    ACE_TRACE(( ACE_TEXT("Gadget::desired_threads (set)") ));

    desired_threads_ = t;
  }
  
  virtual void set_controller(GadgetStreamController* controller) {
    controller_ = controller;
  }

  virtual GadgetStreamController* get_controller()
  {
    return controller_;
  }

  virtual int close(unsigned long flags)
  {
    ACE_TRACE(( ACE_TEXT("Gadget::close") ));
    GADGET_DEBUG2("Gadget (%s) Close Called with flags = %d\n", this->module()->name(), flags);
    int rval = 0;
    if (flags == 1) {
      ACE_Message_Block *hangup = new ACE_Message_Block();
      hangup->msg_type( ACE_Message_Block::MB_HANGUP );
      if (this->putq(hangup) == -1) {
	hangup->release();
	GADGET_DEBUG2("Gadget (%s) failed to put hang up message on queue\n", this->module()->name());
	return GADGET_FAIL;	  
      }
      GADGET_DEBUG2("Gadget (%s) waiting for thread to finish\n", this->module()->name());
      rval = this->wait();
      GADGET_DEBUG2("Gadget (%s) thread finished\n", this->module()->name());
      controller_ = 0;
    }
    return rval;
  }

  virtual int svc(void) 
  {
    ACE_TRACE(( ACE_TEXT("Gadget::svc") ));
    
    for (ACE_Message_Block *m = 0; ;) {

      //GADGET_DEBUG2("Waiting for message in Gadget (%s)\n", this->module()->name());
      if (this->getq(m) == -1) {
	GADGET_DEBUG2("Gadget (%s) failed to get message from queue\n", this->module()->name());
	return GADGET_FAIL;	  
      }
      //GADGET_DEBUG2("Message Received in Gadget (%s)\n", this->module()->name());

      //If this is a hangup message, we are done, put the message back on the queue before breaking
      if (m->msg_type() == ACE_Message_Block::MB_HANGUP) {
	//GADGET_DEBUG2("Gadget (%s) Hangup message encountered\n", this->module()->name());
	if (this->putq(m) == -1) {
	  GADGET_DEBUG2("Gadget (%s) failed to put hang up message on queue (for other threads)\n", this->module()->name());
	  return GADGET_FAIL;	  
	}
	//GADGET_DEBUG2("Gadget (%s) breaking loop\n", this->module()->name());
	break;
      }


      //Is this config info, if so call appropriate process function
      if (m->flags() & GADGET_MESSAGE_CONFIG) {
	if (this->process_config(m) == -1) {
	  m->release();
	  this->flush();
	  GADGET_DEBUG2("Gadget (%s) process config failed\n", this->module()->name());
	  return GADGET_FAIL;

	}

	//Push this onto next gadgets queue, other gadgets may need this configuration information
	if (this->next()) {
	  if (this->next()->putq(m) == -1) {
	    m->release();
	    GADGET_DEBUG2("Gadget (%s) process config failed to put config on dowstream gadget\n", this->module()->name());
	    return GADGET_FAIL;
	  }
	}
	continue;
      }

      if (this->process(m) == -1) {
	m->release();
	this->flush();
	GADGET_DEBUG2("Gadget (%s) process failed\n", this->module()->name());
	return GADGET_FAIL;
      }

    }
    return 0;
  }

  int set_parameter(const char* name, const char* val) {

    parameters_[std::string(name)] = std::string(val);
    return 0;
  }


  int get_bool_value(const char* name) {
    return (0 == ACE_OS::strcmp(get_string_value(name)->c_str(), "true"));
  }

  int get_int_value(const char* name) {
    return ACE_OS::atoi(get_string_value(name)->c_str());
  }

  double get_double_value(const char* name) {
    return ACE_OS::atof(get_string_value(name)->c_str());
  }

  boost::shared_ptr<std::string> get_string_value(const char* name) {
    std::map<std::string,std::string>::iterator it;
    
    it = parameters_.find(std::string(name));
    
    if (it != parameters_.end()) {
      return boost::shared_ptr<std::string>(new std::string(it->second));
    }

    return boost::shared_ptr<std::string>(new std::string(""));
  }


protected:
  virtual int next_step(ACE_Message_Block *m) 
  {
    return this->put_next(m);//next()->putq(m);
  }

  virtual int process(ACE_Message_Block * m) = 0;

  virtual int process_config(ACE_Message_Block * m) {
    return 0;
  }

  unsigned int desired_threads_;
  bool pass_on_undesired_data_;
  GadgetStreamController* controller_;

 private:
  std::map<std::string, std::string> parameters_;
};


class EndGadget : public Gadget
{
protected:
  virtual int process(ACE_Message_Block *m)
  {
    ACE_TRACE(( ACE_TEXT("EndGadget::process(ACE_Message_Block* m)") ));
    m->release();
    return 0;
  }

  virtual int next_step(ACE_Message_Block *m)
  {
    ACE_TRACE(( ACE_TEXT("EndGadget::next_step(ACE_Message_Block *m)") ));
    m->release();
    return 0;
  }

  virtual int process_config(ACE_Message_Block * m) {
    m->release();
    return 0;
  }

};

template <class P1> class Gadget1 : public Gadget
{
  
protected:
  int process(ACE_Message_Block* mb)
  {
    GadgetContainerMessage<P1>* m = AsContainerMessage<P1>(mb);
    
    if (!m) {
      if (!pass_on_undesired_data_) {
	ACE_ERROR_RETURN(( LM_ERROR, ACE_TEXT("%p\n"),
			   ACE_TEXT("Gadget1::process, conversion of message block")),
			 -1);
      } else {
	return (this->next()->putq(mb));
      }

    }
    
    process(m);

    return 0;
  }

  virtual int process(GadgetContainerMessage<P1>* m) = 0;

};

template <class P1, class P2> class Gadget2 : public Gadget
{
  
protected:
  int process(ACE_Message_Block* mb)
  {

    GadgetContainerMessage<P1>* m1 = AsContainerMessage<P1>(mb);
    
    GadgetContainerMessage<P2>* m2 = 0;
    if (m1) {
      m2 = AsContainerMessage<P2>(m1->cont());
    }

    if (!m1 || !m2) {
      if (!pass_on_undesired_data_) {
	ACE_DEBUG( (LM_ERROR, ACE_TEXT("%s -> %s, (%s, %s, %@, %@), (%s, %s, %@, %@)\n"),
		    this->module()->name(),
		    ACE_TEXT("Gadget2::process, Conversion of Message Block Failed"),
		    typeid(GadgetContainerMessage<P1>*).name(),
		    typeid(m1).name(),
		    mb,
		    m1,
		    typeid(GadgetContainerMessage<P2>*).name(),
		    typeid(m2).name(),
		    mb->cont(),
		    m2));
	
	return -1;
      } else {
	return (this->next()->putq(mb));
      }
    }
    
    return this->process(m1,m2);
  }

  virtual int process(GadgetContainerMessage<P1>* m1, GadgetContainerMessage<P2>* m2) = 0;

};


/* Macros for handling dyamic linking */
#define GADGET_DECLARE(GADGET) \
  GADGETRON_LOADABLE_DECLARE(GADGET)

#define GADGET_FACTORY_DECLARE(GADGET)	\
  GADGETRON_LOADABLE_FACTORY_DECLARE(Gadget,GADGET)


#endif //GADGET_H
