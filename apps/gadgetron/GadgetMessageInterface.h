#ifndef GADGETMESSAGEINTERFACE_H
#define GADGETMESSAGEINTERFACE_H

#include "GadgetContainerMessage.h"
#include "GadgetronExport.h"
#include "Gadget.h"

#include <ace/SOCK_Stream.h>
#include <ace/Basic_Types.h>
#include <map>

namespace Gadgetron
{

enum GadgetronMessageID {
  GADGET_MESSAGE_INT_ID_MIN       =   0,
  GADGET_MESSAGE_CONFIG_FILE      =   1,
  GADGET_MESSAGE_CONFIG_SCRIPT    =   2,
  GADGET_MESSAGE_PARAMETER_SCRIPT =   3,
  GADGET_MESSAGE_CLOSE            =   4,
  GADGET_MESSAGE_TEXT             =   5,
  GADGET_MESSAGE_INT_ID_MAX       = 999
};

struct GadgetMessageIdentifier
{
  ACE_UINT16 id;
};

struct GadgetMessageConfigurationFile
{
  char configuration_file[1024];
};

struct GadgetMessageScript
{
  ACE_UINT32 script_length;
};


/**
   Interface for classes capable of reading a specific message

   This is an abstract class, implementations need to be done for each message type.
 */
class GadgetMessageReader
{
 public:
	virtual ~GadgetMessageReader() {}

  /**
     Function must be implemented to read a specific message.
   */
  virtual ACE_Message_Block* read(ACE_SOCK_Stream* stream) = 0;

};

/**
   Interface for classes capable of writing for writing a specific message to a socket. 
   This is an abstract class, implementations need to be done for each message type.
 */
class GadgetMessageWriter
{
 public:
	virtual ~GadgetMessageWriter() {}

   /**
     Function must be implemented to write a specific message.
   */
  virtual int write(ACE_SOCK_Stream* stream, ACE_Message_Block* mb) = 0;
};

class GadgetMessageWriterContainer
{
 public:
  virtual ~GadgetMessageWriterContainer() {
    clear();
  }


  GadgetMessageWriter* find(ACE_UINT16 slot) {
    std::map< ACE_UINT16, GadgetMessageWriter* >::iterator it;

    it = map_.find(slot);
    GadgetMessageWriter* ret = 0;
    if (it != map_.end()) {
      ret = it->second;
    }
    return ret;
  }

  int insert ( unsigned short slot, GadgetMessageWriter* dispatcher) {
    std::map< ACE_UINT16, GadgetMessageWriter* >::iterator it;

    it = map_.find(slot);
    if (it != map_.end()) {
      delete it->second;
      it->second = dispatcher;
    } else {
      map_[slot] = dispatcher;
    }
    return GADGET_OK;
  }

  int clear()
  {
    std::map< ACE_UINT16, GadgetMessageWriter* >::iterator it;
    for (it = map_.begin(); it != map_.end(); it++) {
      delete it->second;
     }
    map_.clear();
    return 0;
  }

 protected:
  std::map<ACE_UINT16, GadgetMessageWriter*> map_;
};


class GadgetMessageReaderContainer
{
 public:
  virtual ~GadgetMessageReaderContainer() {
    clear();
  }


  GadgetMessageReader* find(ACE_UINT16 slot) {
    std::map< ACE_UINT16, GadgetMessageReader* >::iterator it;

    it = map_.find(slot);
    GadgetMessageReader* ret = 0;
    if (it != map_.end()) {
      ret = it->second;
    }
    return ret;
  }

  int insert ( unsigned short slot, GadgetMessageReader* dispatcher) {
    std::map< ACE_UINT16, GadgetMessageReader* >::iterator it;

    it = map_.find(slot);
    if (it != map_.end()) {
      delete it->second;
      it->second = dispatcher;
    } else {
      map_[slot] = dispatcher;
    }
    return GADGET_OK;
  }

  int clear()
  {
    std::map< ACE_UINT16, GadgetMessageReader* >::iterator it;

    for (it = map_.begin(); it != map_.end(); it++) {
      delete it->second;
     }
    map_.clear();
    return 0;
  }
 protected:
  std::map<ACE_UINT16, GadgetMessageReader*> map_;
};

class GadgetMessageConfigFileReader : public GadgetMessageReader
{
 public:
  virtual ACE_Message_Block* read(ACE_SOCK_STREAM* stream) {

    GadgetContainerMessage<GadgetMessageConfigurationFile>* mb1 =
      new GadgetContainerMessage<GadgetMessageConfigurationFile>();
    
    if (!mb1) {
      GDEBUG("Unable to allocate GadgetMessageConfigurationFile\n");
      return 0;
    }

    ssize_t recv_cnt = 0;
    if ((recv_cnt = stream->recv_n (mb1->getObjectPtr(), sizeof(GadgetMessageConfigurationFile))) <= 0) {
      GDEBUG("Unable to read configuration file information\n");
      mb1->release();
      return 0;
    }

    return mb1;
  }
};


class GadgetMessageScriptReader : public GadgetMessageReader
{
 public:
  virtual ACE_Message_Block* read(ACE_SOCK_STREAM* stream) {

    GadgetMessageScript ms;

    ssize_t recv_cnt = 0;
    if ((recv_cnt = stream->recv_n (&ms, sizeof(GadgetMessageScript))) <= 0) {
      GDEBUG("Unable to read configuration file information\n");
       return 0;
    }
    
    ACE_Message_Block* mb = new ACE_Message_Block(ms.script_length);

    if ((recv_cnt = stream->recv_n (mb->wr_ptr(), ms.script_length)) <= 0) {
      GERROR("Unable to read script\n");
      return 0;
    }
    mb->wr_ptr(ms.script_length);
    mb->set_flags(Gadget::GADGET_MESSAGE_CONFIG);

    return mb;
  }
};

/* Macros for handling dyamic linking */

/*#define GADGETRON_READER_DECLARE(READER) \
 GADGETRON_LOADABLE_DECLARE(READER)
*/

#define GADGETRON_READER_DECLARE(READER) 

#define GADGETRON_READER_FACTORY_DECLARE(READER)	\
  GADGETRON_LOADABLE_FACTORY_DECLARE(GadgetMessageReader, READER)

/*#define GADGETRON_WRITER_DECLARE(WRITER) \
 GADGETRON_LOADABLE_DECLARE(WRITER)
*/

#define GADGETRON_WRITER_DECLARE(WRITER) 

#define GADGETRON_WRITER_FACTORY_DECLARE(WRITER)	\
  GADGETRON_LOADABLE_FACTORY_DECLARE(GadgetMessageWriter, WRITER)

}

#endif //GADGETMESSAGEINTERFACE_H
