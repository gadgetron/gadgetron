#include "GadgetIsmrmrdReadWrite.h"
#include "DependencyQueryWriter.h"
#include "GadgetContainerMessage.h"

namespace Gadgetron{

  int DependencyQueryWriter::write(ACE_SOCK_Stream* sock, ACE_Message_Block* mb)
  {
    typedef unsigned long long size_t_type;

    GadgetContainerMessage<ISMRMRD::MetaContainer>* attribmb = AsContainerMessage<ISMRMRD::MetaContainer>(mb);
    if (!attribmb)
      {
	GERROR("DependencyQueryWriter::write, invalid meta attribute message objects\n");
        return -1;
      }

    ssize_t send_cnt = 0;
    GadgetMessageIdentifier id;
    id.id = GADGET_MESSAGE_DEPENDENCY_QUERY;

    if ((send_cnt = sock->send_n (&id, sizeof(GadgetMessageIdentifier))) <= 0)
      {
	GERROR("Unable to send image message identifier\n");
	return -1;
      }

    char* buf = NULL;
    size_t_type len(0);

    try
      {
        std::stringstream str;
        ISMRMRD::serialize( *attribmb->getObjectPtr(), str);
        std::string attribContent = str.str();
        len = attribContent.length()+1;

        buf = new char[len];
        GADGET_CHECK_THROW(buf != NULL);

        memset(buf, '\0', sizeof(char)*len);
        memcpy(buf, attribContent.c_str(), len-1);
      }
    catch(...)
      {
	GERROR("Unable to serialize image meta attributes\n");
        return -1;
      }

    if ( (send_cnt = sock->send_n (&len, sizeof(size_t_type))) <= 0 )
      {
	GERROR("Unable to send image meta attributes length \n");
	if ( buf != NULL ) delete [] buf;
	return -1;
      }

    if ( (send_cnt = sock->send_n (buf, len)) <= 0 )
      {
	GERROR("Unable to send image meta attributes\n");
	if ( buf != NULL ) delete [] buf;
	return -1;
      }

    if ( buf != NULL ) delete [] buf;

    return 0;
  }

  GADGETRON_WRITER_FACTORY_DECLARE(DependencyQueryWriter)
}
