#include "GadgetIsmrmrdReadWrite.h"
#include "DependencyQueryWriter.h"
#include "GadgetContainerMessage.h"

namespace Gadgetron{

int DependencyQueryWriter::write(ACE_SOCK_Stream* sock, ACE_Message_Block* mb)
{
    typedef GtImageAttribType::size_t_type size_t_type;

    GadgetContainerMessage<GtImageAttribType>* attribmb = AsContainerMessage<GtImageAttribType>(mb);
    if (!attribmb)
    {
        ACE_DEBUG( (LM_ERROR, ACE_TEXT("(%P,%l), DependencyQueryWriter::write, invalid meta attribute message objects\n")) );
        return -1;
    }

    ssize_t send_cnt = 0;
    GadgetMessageIdentifier id;
    id.id = GADGET_MESSAGE_DEPENDENCY_QUERY;

    if ((send_cnt = sock->send_n (&id, sizeof(GadgetMessageIdentifier))) <= 0)
    {
        ACE_DEBUG ((LM_ERROR,
                ACE_TEXT ("(%P|%t) Unable to send image message identifier\n")));

        return -1;
    }

    char* buf = NULL;
    size_t_type len(0);

    if ( !attribmb->getObjectPtr()->serialize(buf, len) )
    {
        ACE_DEBUG ((LM_ERROR, ACE_TEXT ("(%P|%t) Unable to serialize image meta attributes \n")));

        return -1;
    }

    if ( (send_cnt = sock->send_n (buf, len)) <= 0 )
    {
        ACE_DEBUG ((LM_ERROR,
                ACE_TEXT ("(%P|%t) Unable to send image meta attributes\n")));

        if ( buf != NULL ) delete [] buf;

        return -1;
    }

    if ( buf != NULL ) delete [] buf;

    return 0;
}

GADGETRON_WRITER_FACTORY_DECLARE(DependencyQueryWriter)
}
