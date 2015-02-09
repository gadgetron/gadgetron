#ifndef GADGETMESSAGEREADERWRITER_H
#define GADGETMESSAGEREADERWRITER_H

#include <ace/SOCK_Stream.h>
#include <ace/Basic_Types.h>

namespace Gadgetron
{

enum GadgetronMessageID {
  GADGET_MESSAGE_INT_ID_MIN       =   0,
  GADGET_MESSAGE_CONFIG_FILE      =   1,
  GADGET_MESSAGE_CONFIG_SCRIPT    =   2,
  GADGET_MESSAGE_PARAMETER_SCRIPT =   3,
  GADGET_MESSAGE_CLOSE            =   4,
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


/* Macros for handling dyamic linking */

//#define GADGETRON_READER_DECLARE(READER) \
//  GADGETRON_LOADABLE_DECLARE(READER)

#define GADGETRON_READER_DECLARE(READER) 

#define GADGETRON_READER_FACTORY_DECLARE(READER)	\
  GADGETRON_LOADABLE_FACTORY_DECLARE(GadgetMessageReader, READER)

//#define GADGETRON_WRITER_DECLARE(WRITER) \
//  GADGETRON_LOADABLE_DECLARE(WRITER)

#define GADGETRON_WRITER_DECLARE(WRITER) 

#define GADGETRON_WRITER_FACTORY_DECLARE(WRITER)	\
  GADGETRON_LOADABLE_FACTORY_DECLARE(GadgetMessageWriter, WRITER)

}

#endif
