#ifndef GADGETSOCKETSENDER_H
#define GADGETSOCKETSENDER_H

#include "ace/SOCK_Stream.h"
#include "ace/Task.h"

#include <complex>

#include "GadgetMRIHeaders.h"
#include "hoNDArray.h"
#include "GadgetContainerMessage.h"
#include "GadgetMessageInterface.h"

namespace Gadgetron
{

/**
Default implementation of GadgetMessageWriter for Image messages
*/

template <typename T> class GadgetImageMessageWriter : public GadgetMessageWriter
{
public:
    virtual int write(ACE_SOCK_Stream* sock, ACE_Message_Block* mb) 
    {
        GadgetContainerMessage<ISMRMRD::ImageHeader>* imagemb = 
            dynamic_cast< GadgetContainerMessage<ISMRMRD::ImageHeader>* >(mb);

        GadgetContainerMessage< hoNDArray< T > >* datamb =
            dynamic_cast< GadgetContainerMessage< hoNDArray< T > >* >(imagemb->cont());

        if (!imagemb || !datamb) {
            ACE_DEBUG( (LM_ERROR, ACE_TEXT("(%P,%l), GadgetImageMessageWriter invalid image message objects")) );
            return -1;
        }


        ssize_t send_cnt = 0;
        GadgetMessageIdentifier id;

        switch (sizeof(T)) {
        case 2: //Unsigned short
            id.id = GADGET_MESSAGE_IMAGE_REAL_USHORT;
            break;
        case 4: //Float
            id.id = GADGET_MESSAGE_IMAGE_REAL_FLOAT;
            break;
        case 8: //Complex float
            id.id = GADGET_MESSAGE_IMAGE_CPLX_FLOAT;
            break;
        default:
            ACE_DEBUG( (LM_ERROR, ACE_TEXT("(%P,%l), GadgetImageMessageWriter Wrong data size detected:")) );
            return -1;
        }

        if ((send_cnt = sock->send_n (&id, sizeof(GadgetMessageIdentifier))) <= 0) {
            ACE_DEBUG ((LM_ERROR,
                ACE_TEXT ("(%P|%t) Unable to send image message identifier\n")));

            return -1;
        }

        if ((send_cnt = sock->send_n (imagemb->getObjectPtr(), sizeof(ISMRMRD::ImageHeader))) <= 0) {
            ACE_DEBUG ((LM_ERROR,
                ACE_TEXT ("(%P|%t) Unable to send image header\n")));

            return -1;
        }

        if ((send_cnt = sock->send_n (datamb->getObjectPtr()->get_data_ptr(), sizeof(T)*datamb->getObjectPtr()->get_number_of_elements())) <= 0) {
            ACE_DEBUG ((LM_ERROR,
                ACE_TEXT ("(%P|%t) Unable to send image data\n")));

            return -1;
        }

        return 0;
    }

};

}

#endif //GADGETSOCKETSENDER_H
