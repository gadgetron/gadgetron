#ifndef GADGETSOCKETRECEIVER_H
#define GADGETSOCKETRECEIVER_H

#include "ace/SOCK_Stream.h"
#include "ace/Task.h"

#include <complex>

#include "GadgetMRIHeaders.h"
#include "hoNDArray.h"
#include "GadgetMessageInterface.h"


/**
   Default implementation of GadgetMessageReader for Acquisition messages

 */
class GadgetAcquisitionMessageReader : public GadgetMessageReader
{

 public:
  virtual ACE_Message_Block* read(ACE_SOCK_Stream* stream) 
  {

    GadgetContainerMessage<GadgetMessageAcquisition>* acqh =
      new GadgetContainerMessage<GadgetMessageAcquisition>();

    ssize_t recv_count = 0;

    if ((recv_count = stream->recv_n(acqh->getObjectPtr(), sizeof(GadgetMessageAcquisition))) <= 0) {
      ACE_DEBUG( (LM_ERROR, ACE_TEXT("%P, %l, GadgetAcquisitionMessageReader, failed to read ACQ Header\n")) );
      acqh->release();
      return 0;
    }

    std::vector<unsigned int> dims(1);dims[0] = acqh->getObjectPtr()->samples;
    GadgetContainerMessage< hoNDArray< std::complex<float> > >* data = 
      new GadgetContainerMessage< hoNDArray< std::complex<float> > >();

    if (!data->getObjectPtr()->create(&dims)) {
      ACE_DEBUG( (LM_ERROR, ACE_TEXT("%P, %l, GadgetAcquisitionMessageReader, failed to allocate memory\n")) );
      acqh->release();
      return 0;
    }

    acqh->cont(data);

    if ((recv_count = stream->recv_n(data->getObjectPtr()->get_data_ptr(), sizeof(float)*2*acqh->getObjectPtr()->samples)) <= 0) {
      ACE_DEBUG( (LM_ERROR, ACE_TEXT("%P, %l, GadgetAcquisitionMessageReader, failed to read data from socket\n")) );
      acqh->release();
      return 0;
    }

    return acqh;
  }
};

/**
   Default implementation of GadgetMessageReader for Image messages

 */
class GadgetImageMessageReader : public GadgetMessageReader
{

 public:
  virtual ACE_Message_Block* read(ACE_SOCK_Stream* stream) 
  {
    GadgetContainerMessage<GadgetMessageImage>* imgh = 
      new GadgetContainerMessage<GadgetMessageImage>();

    ssize_t recv_count = 0;
    if ((recv_count = stream->recv_n(imgh->getObjectPtr(), sizeof(GadgetMessageImage))) <= 0) {
      ACE_DEBUG( (LM_ERROR, ACE_TEXT("%P, %l, GadgetImageMessageReader, failed to read IMAGE Header\n")) );
      imgh->release();
      return 0;
    }

    std::vector<unsigned int> dims(3);
    dims[0] = imgh->getObjectPtr()->matrix_size[0];
    dims[1] = imgh->getObjectPtr()->matrix_size[1];
    dims[2] = imgh->getObjectPtr()->matrix_size[2];

    if (imgh->getObjectPtr()->channels > 1) {
      dims.push_back(imgh->getObjectPtr()->channels);
    } 

    GadgetContainerMessage< hoNDArray< std::complex<float> > >* data = 
      new GadgetContainerMessage< hoNDArray< std::complex<float> > >();

    if (!data->getObjectPtr()->create(&dims)) {
      ACE_DEBUG( (LM_ERROR, 
		  ACE_TEXT("%P, %l, GadgetImageMessageReader, failed to allocate memory\n")) );
      imgh->release();
      return 0;
    }

    imgh->cont(data);

    if ((recv_count = stream->recv_n(data->getObjectPtr()->get_data_ptr(), sizeof(float)*2*data->getObjectPtr()->get_number_of_elements())) <= 0) {
      ACE_DEBUG( (LM_ERROR, ACE_TEXT("%P, %l, GadgetImageMessageReader, failed to read data from socket\n")) );
      imgh->release();
      return 0;
    }
    
    return imgh;
  }
};

class GadgetSocketReceiver : public ACE_Task<ACE_MT_SYNCH>
{

 public:
  typedef ACE_Task<ACE_MT_SYNCH> inherited;

  GadgetSocketReceiver(ACE_SOCK_Stream* socket) 
    : inherited()
    , desired_threads_(1)
    , socket_(socket)
  {
    ACE_TRACE(( ACE_TEXT("GadgetSocketReceiver::GadgetSocketReceiver") ));
  }

  virtual ~GadgetSocketReceiver() 
    {
      for (unsigned int i = 0; i < readers_.size(); i++) {
	if (readers_[i]) delete readers_[i]; 
      }  
    }

  virtual int init(void)
  {
    ACE_TRACE(( ACE_TEXT("Gadget::init") ));
    return 0;
  }

  virtual int open(void* = 0) 
  {
    ACE_TRACE(( ACE_TEXT("GagetSocketReceiver::open") ));
    
    return this->activate( THR_NEW_LWP | THR_JOINABLE,
			   this->desired_threads() );
  }


  virtual int register_reader(ACE_UINT16 slot, GadgetMessageReader* reader) 
  {
    if (reader == 0) {
      return 0;
    }

    if (find_reader(slot) != 0) {
      ACE_DEBUG( (LM_ERROR, ACE_TEXT("%P, %l, GagetSocketReceiver, attempting to register reader in occupied slot\n")) );
      return -1;
    }
 
    slots_.push_back(slot);
    readers_.push_back(reader);
    return 0;
  }

  virtual unsigned int desired_threads()
  {
    ACE_TRACE(( ACE_TEXT("GagetSocketReceiver::desired_threads (get)") ));
    return desired_threads_;
  }

  virtual void desired_threads(unsigned int t)
  {
    ACE_TRACE(( ACE_TEXT("GagetSocketReceiver::desired_threads (set)") ));
    desired_threads_ = t;
  }
  
  virtual int close(unsigned long flags)
  {
    ACE_TRACE(( ACE_TEXT("GagetSocketReceiver::close") ));
    return 0;
  }

  virtual int svc(void) 
  {
    ACE_TRACE(( ACE_TEXT("GadgetSocketReceiver::svc") ));
    ssize_t recv_count = 0;
    GadgetMessageIdentifier mid;
    while (1) {   

      if ((recv_count = socket_->recv_n(&mid, sizeof(GadgetMessageIdentifier))) <= 0) {
	ACE_DEBUG( (LM_ERROR, ACE_TEXT("%P, %l, GadgetSocketReceiver, failed to read message identifier\n")) );
	return -1;
      }
      
      GadgetMessageReader* r = find_reader(mid.id);
      if (r == 0) {
	ACE_DEBUG( (LM_ERROR, ACE_TEXT("%P, %l, GadgetSocketReceiver, Unknown message id %d received\n"), mid.id) );
	return -1;
      }

      ACE_Message_Block* mb = r->read(socket_);

      if (!mb) {
	ACE_DEBUG( (LM_ERROR, ACE_TEXT("%P, %l, GadgetSocketReceiver, Failed to read message\n")) );
	return -1;
      }	else {
	if (process(mid.id, mb) < 0) {
	  ACE_DEBUG( (LM_ERROR, ACE_TEXT("%P, %l, GadgetSocketReceiver, Failed to process message\n")) );	  
	  return -1;;
	}
      }
    }

    return 0;
  }


 protected:
  unsigned int desired_threads_;
  ACE_SOCK_Stream* socket_;
  std::vector<ACE_UINT16> slots_;
  std::vector<GadgetMessageReader*> readers_;

  GadgetMessageReader* find_reader(ACE_UINT16 slot) {
    GadgetMessageReader* ret = 0;
    for (unsigned int i = 0; i < slots_.size(); i++) {
      if (slots_[i] == slot) {ret = readers_[i]; break;}
    }
    return ret;
  }

  //Default implementation. To be overwritten by sub classes.
  int process(ACE_UINT16 id, ACE_Message_Block* mb) {
    mb->release();
    return 0;
  }
};

#endif //GADGETSOCKETRECEIVER_H
