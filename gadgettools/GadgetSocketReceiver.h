#ifndef GADGETSOCKETRECEIVER_H
#define GADGETSOCKETRECEIVER_H

#include "ace/SOCK_Stream.h"
#include "ace/Task.h"

#include <complex>

#include "../gadgetheaders.h"
#include "NDArray.h"

/**
   Class for reading a specific message of a socket. 
   This is an abstract class, implementations need to be done for each message type.
 */
class GadgetMessageReader
{
 public:
  /**
     Function must be implemented to read a specific message.
   */
  virtual int read(ACE_SOCK_Stream* stream) = 0;
};

/**
   Default implementation of GadgetMessageReader for Acquisition messages

 */
class GadgetAcquisitionMessageReader : public GadgetMessageReader
{

 public:
  virtual int read(ACE_SOCK_Stream* stream) 
  {
    GadgetMessageAcquisition acqh;
    ssize_t recv_count = 0;

    if ((recv_count = stream->recv_n(&acqh, sizeof(GadgetMessageAcquisition))) <= 0) {
      ACE_DEBUG( (LM_ERROR, ACE_TEXT("%P, %l, GadgetAcquisitionMessageReader, failed to read ACQ Header\n")) );
      return -1;
    }

    std::vector<int> dims(1);dims[0] = acqh.samples;
    NDArray< std::complex<float> > data;
    if (!data.create(dims)) {
      ACE_DEBUG( (LM_ERROR, ACE_TEXT("%P, %l, GadgetAcquisitionMessageReader, failed to allocate memory\n")) );
      return -1;
    }

    if ((recv_count = stream->recv_n(data.get_data_ptr(), sizeof(float)*2*acqh.samples)) <= 0) {
      ACE_DEBUG( (LM_ERROR, ACE_TEXT("%P, %l, GadgetAcquisitionMessageReader, failed to read data from socket\n")) );
      return -1;
    }

    if (process_acquisition(&acqh,&data) < 0) {
      ACE_DEBUG( (LM_ERROR, ACE_TEXT("%P, %l, GadgetAcquisitionMessageReader, failed to processdata\n")) );
      return -1;
    }

    return 0;
  }

  virtual int process_acquisition(GadgetMessageAcquisition* acq_head, NDArray< std::complex<float> >* data)
  {
    ACE_DEBUG( (LM_DEBUG, ACE_TEXT("GadgetAcquisitionMessageReader processing ACQUISITION\n")) );
    return 0;
  }
};

/**
   Default implementation of GadgetMessageReader for Image messages

 */
class GadgetImageMessageReader : public GadgetMessageReader
{

 public:
  virtual int read(ACE_SOCK_Stream* stream) 
  {
    GadgetMessageImage imgh;
    ssize_t recv_count = 0;

    if ((recv_count = stream->recv_n(&imgh, sizeof(GadgetMessageImage))) <= 0) {
      ACE_DEBUG( (LM_ERROR, ACE_TEXT("%P, %l, GadgetImageMessageReader, failed to read IMAGE Header\n")) );
      return -1;
    }

    std::vector<int> dims(3);
    dims[0] = imgh.matrix_size[0];dims[1] = imgh.matrix_size[1];dims[2] = imgh.matrix_size[2];
    NDArray< std::complex<float> > data;
    if (!data.create(dims)) {
      ACE_DEBUG( (LM_ERROR, ACE_TEXT("%P, %l, GadgetImageMessageReader, failed to allocate memory\n")) );
      return -1;
    }

    if ((recv_count = stream->recv_n(data.get_data_ptr(), sizeof(float)*2*data.get_number_of_elements())) <= 0) {
      ACE_DEBUG( (LM_ERROR, ACE_TEXT("%P, %l, GadgetImageMessageReader, failed to read data from socket\n")) );
      return -1;
    }

    if (process_image(&imgh,&data) < 0) {
      ACE_DEBUG( (LM_ERROR, ACE_TEXT("%P, %l, GadgetImageMessageReader, failed to read data\n")) );
      return -1;
    }

    return 0;
  }

  virtual int process_image(GadgetMessageImage* img_head, NDArray< std::complex<float> >* data)
  {
    ACE_DEBUG( (LM_DEBUG, ACE_TEXT("GadgetImagenMessageReader processing IMAGE\n")) );
    return 0;
  }
};

#if 1
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

      if (r->read(socket_) < 0) {
	ACE_DEBUG( (LM_ERROR, ACE_TEXT("%P, %l, GadgetSocketReceiver, Failed to read message\n")) );
	return -1;
      }	
    }
  }


 protected:
  unsigned int desired_threads_;
  ACE_SOCK_Stream* socket_;
  std::vector<ACE_UINT16> slots_;
  std::vector<GadgetMessageReader*> readers_;

  GadgetMessageReader* find_reader(ACE_UINT16 slot) {
    GadgetMessageReader* ret = 0;
    for (unsigned int i = 0; i < slots_.size(); i++) {
      if (slots_[i] == slot) ret = readers_[i]; break;
    }
    return ret;
  }
};
#endif

#endif //GADGETSOCKETRECEIVER_H
