/** \file   MRIImageWriter.h
\brief  MRI image writer with or without meta attributes.
\author Hui Xue
*/

#ifndef MRIImageWriter_H
#define MRIImageWriter_H

#include "GadgetMessageInterface.h"
#include "GadgetMRIHeaders.h"
#include "ismrmrd/meta.h"
#include "gadgetron_mricore_export.h"

#include <ismrmrd/ismrmrd.h>
#include <complex>

namespace Gadgetron{

    class MRIImageWriter : public GadgetMessageWriter
    {
    public:
        virtual int write(ACE_SOCK_Stream* sock, ACE_Message_Block* mb);

        template <typename T>
        int write_data_attrib(ACE_SOCK_Stream* sock, GadgetContainerMessage<ISMRMRD::ImageHeader>* header, GadgetContainerMessage< hoNDArray<T> >* data)
        {
            typedef unsigned long long size_t_type;

            uint16_t RO = header->getObjectPtr()->matrix_size[0];
            uint16_t E1 = header->getObjectPtr()->matrix_size[1];
            uint16_t E2 = header->getObjectPtr()->matrix_size[2];
            uint16_t CHA = header->getObjectPtr()->channels;

            unsigned long expected_elements = RO*E1*E2*CHA;

            if (expected_elements != data->getObjectPtr()->get_number_of_elements())
            {
                GDEBUG("Number of header elements %d is inconsistent with number of elements in NDArray %d\n", expected_elements, data->getObjectPtr()->get_number_of_elements());
                GDEBUG("Header dimensions: %d, %d, %d, %d\n", RO, E1, E2, CHA);
                GDEBUG("Number of array dimensions: %d:\n", data->getObjectPtr()->get_number_of_dimensions());
                for (size_t i = 0; i < data->getObjectPtr()->get_number_of_dimensions(); i++)
                {
                    GDEBUG("Dimensions %d: %d\n", i, data->getObjectPtr()->get_size(i));
                }
                return -1;
            }

            ssize_t send_cnt = 0;
            GadgetMessageIdentifier id;
            id.id = GADGET_MESSAGE_ISMRMRD_IMAGE;

            if ((send_cnt = sock->send_n(&id, sizeof(GadgetMessageIdentifier))) <= 0)
            {
                GERROR("Unable to send image message identifier\n");
                return -1;
            }

            GadgetContainerMessage<ISMRMRD::MetaContainer>* attribmb = AsContainerMessage<ISMRMRD::MetaContainer>(data->cont());

            char* buf = NULL;
            size_t_type len(0);

            if (attribmb)
            {
                try
                {
                    std::stringstream str;
                    ISMRMRD::serialize(*attribmb->getObjectPtr(), str);
                    std::string attribContent = str.str();
                    len = attribContent.length() + 1;

                    buf = new char[len];
                    GADGET_CHECK_THROW(buf != NULL);

                    memset(buf, '\0', sizeof(char)*len);
                    memcpy(buf, attribContent.c_str(), len - 1);
                }
                catch (...)
                {
                    GERROR("Unable to serialize image meta attributes \n");
                    return -1;
                }
            }

            header->getObjectPtr()->attribute_string_len = (uint32_t)len;

            if ((send_cnt = sock->send_n(header->getObjectPtr(), sizeof(ISMRMRD::ImageHeader))) <= 0)
            {
                GERROR("Unable to send image header\n");
                return -1;
            }

            if ((send_cnt = sock->send_n(&len, sizeof(size_t_type))) <= 0)
            {
                GERROR("Unable to send image meta attributes length\n");
                if (buf != NULL) delete[] buf;
                return -1;
            }

            if (len>0)
            {
                if ((send_cnt = sock->send_n(buf, len)) <= 0)
                {
                    GERROR("Unable to send image meta attributes\n");
                    if (buf != NULL) delete[] buf;
                    return -1;
                }
            }

            if (buf != NULL) delete[] buf;

            if ((send_cnt = sock->send_n(data->getObjectPtr()->get_data_ptr(), sizeof(T)*data->getObjectPtr()->get_number_of_elements())) <= 0)
            {
                GERROR("Unable to send image data\n");
                return -1;
            }

            return 0;
        }
    };

}
#endif
