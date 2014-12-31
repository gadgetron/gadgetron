#include "GadgetIsmrmrdReadWrite.h"
#include "MRIImageAttribWriter.h"
#include "GadgetContainerMessage.h"
#include "hoNDArray.h"

#include <complex>

namespace Gadgetron{

template <typename T>
int MRIImageAttribWriter<T>::write(ACE_SOCK_Stream* sock, ACE_Message_Block* mb)
{
    typedef unsigned long long size_t_type;

    GadgetContainerMessage<ISMRMRD::ImageHeader>* imagemb =
            AsContainerMessage<ISMRMRD::ImageHeader>(mb);

    if (!imagemb)
    {
        ACE_DEBUG( (LM_ERROR, ACE_TEXT("(%P,%l), MRIImageAttribWriter::write, invalid image message objects, 1\n")) );
        return -1;
    }

    GadgetContainerMessage< hoNDArray< T > >* datamb =
            AsContainerMessage< hoNDArray< T > >(imagemb->cont());

    if (!datamb)
    {
        ACE_DEBUG( (LM_ERROR, ACE_TEXT("(%P,%l), MRIImageAttribWriter::write, invalid image message objects\n")) );
        return -1;
    }

    GadgetContainerMessage<ISMRMRD::MetaContainer>* attribmb =
            AsContainerMessage<ISMRMRD::MetaContainer>(datamb->cont());

    if (!attribmb)
    {
        ACE_DEBUG( (LM_ERROR, ACE_TEXT("(%P,%l), MRIImageAttribWriter::write, invalid image attribute message objects\n")) );
        return -1;
    }

    ssize_t send_cnt = 0;
    GadgetMessageIdentifier id;
    switch (sizeof(T))
    {
    case 2: //Unsigned short
        id.id = GADGET_MESSAGE_ISMRMRD_IMAGEWITHATTRIB_REAL_USHORT;
        break;
    case 4: //Float
        id.id = GADGET_MESSAGE_ISMRMRD_IMAGEWITHATTRIB_REAL_FLOAT;
        break;
    case 8: //Complex float
        id.id = GADGET_MESSAGE_ISMRMRD_IMAGEWITHATTRIB_CPLX_FLOAT;
        break;
    default:
        ACE_DEBUG ((LM_ERROR,
                ACE_TEXT ("(%P|%t) MRIImageAttribWriter Wrong data size detected\n")));
        return GADGET_FAIL;
    }

    //Let's check if the image header is consistent with the data array size before sending:
    uint16_t RO = imagemb->getObjectPtr()->matrix_size[0];
    uint16_t E1 = imagemb->getObjectPtr()->matrix_size[1];
    uint16_t E2 = imagemb->getObjectPtr()->matrix_size[2];

    unsigned long expected_elements = RO*E1*E2;

    if (expected_elements !=  datamb->getObjectPtr()->get_number_of_elements())
    {
        GDEBUG("Number of header elements %d is inconsistent with number of elements in NDArray %d\n",expected_elements, datamb->getObjectPtr()->get_number_of_elements());
        GDEBUG("Header dimensions: %d, %d, %d\n",RO,E1,E2);
        GDEBUG("Number of array dimensions: %d:\n", datamb->getObjectPtr()->get_number_of_dimensions());
        for (size_t i = 0; i < datamb->getObjectPtr()->get_number_of_dimensions(); i++)
        {
            GDEBUG("Dimensions %d: %d\n", i, datamb->getObjectPtr()->get_size(i));
        }
        return -1;
    }

    if ((send_cnt = sock->send_n (&id, sizeof(GadgetMessageIdentifier))) <= 0)
    {
        ACE_DEBUG ((LM_ERROR,
                ACE_TEXT ("(%P|%t) Unable to send image message identifier\n")));

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
        ACE_DEBUG ((LM_ERROR, ACE_TEXT ("(%P|%t) Unable to serialize image meta attributes \n")));

        return -1;
    }

    imagemb->getObjectPtr()->attribute_string_len = len;

    if ((send_cnt = sock->send_n ( imagemb->getObjectPtr(), sizeof(ISMRMRD::ImageHeader))) <= 0)
    {
        ACE_DEBUG ((LM_ERROR,
                ACE_TEXT ("(%P|%t) Unable to send image header\n")));

        return -1;
    }

    if ( (send_cnt = sock->send_n (&len, sizeof(size_t_type))) <= 0 )
    {
        ACE_DEBUG ((LM_ERROR, ACE_TEXT ("(%P|%t) Unable to send image meta attributes length \n")));
        if ( buf != NULL ) delete [] buf;
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

    if ((send_cnt = sock->send_n (datamb->getObjectPtr()->get_data_ptr(), sizeof(T)*datamb->getObjectPtr()->get_number_of_elements())) <= 0)
    {
        ACE_DEBUG ((LM_ERROR,
                ACE_TEXT ("(%P|%t) Unable to send image data\n")));

        return -1;
    }

    return 0;
}

GADGETRON_WRITER_FACTORY_DECLARE(MRIImageAttribWriterFLOAT)
GADGETRON_WRITER_FACTORY_DECLARE(MRIImageAttribWriterUSHORT)
GADGETRON_WRITER_FACTORY_DECLARE(MRIImageAttribWriterCPLX)
}
