#ifndef GADGETMESSAGETEXTWRITER_H
#define GADGETMESSAGETEXTWRITER_H

#include "GadgetContainerMessage.h"
#include "GadgetMessageInterface.h"
#include "gadgetron_util_gadgets_export.h"

namespace Gadgetron{

    class EXPORTUTILGADGETS GadgetMessageTextWriter : public GadgetMessageWriter
    {

    public:
        virtual int write(ACE_SOCK_Stream* sock, ACE_Message_Block* mb)
        {
	  if (!mb) {
	    GERROR("GadgetMessageTextWriter, invalid acquisition message objects");
	    return -1;
	  }

	  GadgetMessageIdentifier id;
	  id.id = GADGET_MESSAGE_TEXT;

	  ssize_t send_cnt = 0;

	  if ((send_cnt = sock->send_n (&id, sizeof(GadgetMessageIdentifier))) <= 0) {
	    GERROR("Unable to send acquisition message identifier\n");
	    return -1;
	  }

	  uint32_t text_len = (uint32_t)mb->size();
	  if ((send_cnt = sock->send_n (&text_len, sizeof(uint32_t))) <= 0) {
	    GERROR("Unable to send length of text");
	    return -1;
	  }
	  
	  if ((send_cnt = sock->send_n (mb->rd_ptr(), text_len)) <= 0) {
	    GERROR("Unable to send text of length %d\n", text_len);
	    return -1;
	  }

	  return 0;
        }
    };
}

#endif //GADGETMESSAGETEXTWRITER_H
