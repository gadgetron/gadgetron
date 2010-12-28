#ifndef GADGET_H
#define GADGET_H

#include <ace/Task.h>
#include <ace/Stream.h>
#include <ace/Module.h>

#include "GadgetContainerMessage.h"

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
  {
    ACE_TRACE(( ACE_TEXT("Gadget::Gadget") ));
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
  
  virtual int close(unsigned long flags)
  {
    ACE_TRACE(( ACE_TEXT("Gadget::close") ));

    int rval = 0;
    if (flags == 1) {
      ACE_Message_Block *hangup = new ACE_Message_Block();
      hangup->msg_type( ACE_Message_Block::MB_HANGUP );
      if (this->putq(hangup) == -1) {
	hangup->release();
	ACE_ERROR_RETURN( (LM_ERROR,
			   ACE_TEXT("%p\n"),
			   ACE_TEXT("Gadget::close, putq")),
			  -1);
      }
      rval = this->wait();
    }
    return rval;
  }

  virtual int svc(void) 
  {
    ACE_TRACE(( ACE_TEXT("Gadget::svc") ));
    
    for (ACE_Message_Block *m = 0; ;) {

      if (this->getq(m) == -1) {
	ACE_ERROR_RETURN( (LM_ERROR, ACE_TEXT("%p\n"),
			   ACE_TEXT("Gadget::getq")),
			  -1);
      }

      //If this is a hangup message, we are done, put the message back on the queue before breaking
      if (m->msg_type() == ACE_Message_Block::MB_HANGUP) {
	if (this->putq(m) == -1) {
	  ACE_ERROR_RETURN( (LM_ERROR,
			     ACE_TEXT("%p\n"),
			     ACE_TEXT("Gadget::svc, putq")),
			    -1);
	}
	break;
      }


      //Is this config info, if so call appropriate process function
      if (m->flags() & GADGET_MESSAGE_CONFIG) {
	if (this->process_config(m) == -1) {
	  m->release();
	  ACE_ERROR_RETURN( (LM_ERROR,
			     ACE_TEXT("%p\n"),
			     ACE_TEXT("Gadget::svc, process_config")),
			    -1);
	  
	}

	//Push this onto next gadgets queue, other gadgets may need this configuration information
	if (this->next()) {
	  if (this->next()->putq(m) == -1) {
	    m->release();
	    ACE_ERROR_RETURN( (LM_ERROR,
			       ACE_TEXT("%p\n"),
			       ACE_TEXT("Gadget::svc, passing config on to next gadget")),
			      -1);
	  }
	}
	continue;
      }

      if (this->process(m) == -1) {
	m->release();
	ACE_ERROR_RETURN( (LM_ERROR,
			   ACE_TEXT("%p\n"),
			   ACE_TEXT("Gadget::svc, process")),
			  -1);
      }

    }
    return 0;
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
};


class EndGadget : public Gadget
{
protected:
  virtual int process(ACE_Message_Block *m)
  {
    ACE_TRACE(( ACE_TEXT("EndGadget::process(ACE_Message_Block* m)") ));
    
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
    GadgetContainerMessage<P1>* m = dynamic_cast< GadgetContainerMessage<P1>* >(mb);
    
    if (!m) {
      ACE_ERROR_RETURN(( LM_ERROR, ACE_TEXT("%p\n"),
			 ACE_TEXT("Gadget1::process, dynamic cast failed")),
		       -1);
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

    GadgetContainerMessage<P1>* m1 = dynamic_cast< GadgetContainerMessage<P1>* >(mb);
    
    if (!m1) {
      ACE_DEBUG( (LM_ERROR, ACE_TEXT("%s -> %s, %s, %s, %@, %@\n"),
		  this->module()->name(),
		  ACE_TEXT("Gadget2::process, dynamic cast failed (arg 0)"),
		  typeid(GadgetContainerMessage<P1>*).name(),
		  typeid(m1).name(),
		  mb,
		  m1));

      return -1;
    }
    
    if (!m1->cont()) {
      ACE_ERROR_RETURN(( LM_ERROR, ACE_TEXT("%p\n"),
			 ACE_TEXT("Gadget2::process, missing argument")),
		       -1);
    }

    GadgetContainerMessage<P2>* m2 = dynamic_cast< GadgetContainerMessage<P2>* >(m1->cont());
    
    if (!m2) {
      ACE_ERROR_RETURN(( LM_ERROR, ACE_TEXT("%s\n"),
			 this->module()->name(),
			 ACE_TEXT("Gadget2::process, dynamic cast failed (arg 1)")),
		       -1);
    }
    
    return this->process(m1,m2);
  }

  virtual int process(GadgetContainerMessage<P1>* m1, GadgetContainerMessage<P2>* m2) = 0;

};

#endif //GADGET_H
