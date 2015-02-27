#ifndef GADGETMESSAGEINTERFACE_H
#define GADGETMESSAGEINTERFACE_H

#include "GadgetContainerMessage.h"
#include "GadgetronExport.h"
#include "Gadget.h"
#include "GadgetMessageReaderWriter.h"

#include <ace/SOCK_Stream.h>
#include <ace/Basic_Types.h>
#include <map>

namespace Gadgetron
{

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

}

#endif //GADGETMESSAGEINTERFACE_H
