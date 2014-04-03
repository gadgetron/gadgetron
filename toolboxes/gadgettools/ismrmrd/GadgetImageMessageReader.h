#ifndef GADGETSOCKETRECEIVER_H
#define GADGETSOCKETRECEIVER_H

#include "ace/SOCK_Stream.h"
#include "ace/Task.h"

#include <complex>
#include <iostream>

#include "GadgetMRIHeaders.h"
#include "ismrmrd.h"
#include "hoNDArray.h"
#include "GadgetMessageInterface.h"

namespace Gadgetron
{

/**
Default implementation of GadgetMessageReader for Image messages
*/

template <typename T> class GadgetImageMessageReader : public GadgetMessageReader
{

public:
    virtual ACE_Message_Block* read(ACE_SOCK_Stream* stream) 
    {
        GadgetContainerMessage<ISMRMRD::ImageHeader>* imgh = 
            new GadgetContainerMessage<ISMRMRD::ImageHeader>();

        ssize_t recv_count = 0;
        if ((recv_count = stream->recv_n(imgh->getObjectPtr(), sizeof(ISMRMRD::ImageHeader))) <= 0) {
            ACE_DEBUG( (LM_ERROR, ACE_TEXT("%P, %l, GadgetImageMessageReader, failed to read IMAGE Header\n")) );
            imgh->release();
            return 0;
        }

        std::vector<size_t> dims(3);
        dims[0] = imgh->getObjectPtr()->matrix_size[0];
        dims[1] = imgh->getObjectPtr()->matrix_size[1];
        dims[2] = imgh->getObjectPtr()->matrix_size[2];

        if (imgh->getObjectPtr()->channels > 1) {
            dims.push_back(imgh->getObjectPtr()->channels);
        } 

        GadgetContainerMessage< hoNDArray< T > >* data =
            new GadgetContainerMessage< hoNDArray< T > >();

        try{ data->getObjectPtr()->create(&dims);}
        catch (std::runtime_error &err){
            GADGET_DEBUG_EXCEPTION(err,"GadgetImageMessageReader, failed to allocate memory\n");
            imgh->release();
            return 0;
        }

        imgh->cont(data);

        if ((recv_count = stream->recv_n(data->getObjectPtr()->get_data_ptr(), sizeof(T)*data->getObjectPtr()->get_number_of_elements())) <= 0) {
            ACE_DEBUG( (LM_ERROR, ACE_TEXT("%P, %l, GadgetImageMessageReader, failed to read data from socket\n")) );
            imgh->release();
            return 0;
        }

        return imgh;
    }
};

}

#endif //GADGETSOCKETRECEIVER_H
