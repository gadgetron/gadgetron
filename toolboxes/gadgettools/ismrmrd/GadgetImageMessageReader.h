#ifndef GADGETSOCKETRECEIVER_H
#define GADGETSOCKETRECEIVER_H

#include "ace/SOCK_Stream.h"
#include "ace/Task.h"

#include <complex>
#include <iostream>

#include "GadgetMRIHeaders.h"
#include "ismrmrd/ismrmrd.h"
#include "hoNDArray.h"
#include "GadgetMessageInterface.h"
#include "ismrmrd/meta.h"

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
	  GERROR("GadgetImageMessageReader, failed to read IMAGE Header\n");
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
            GEXCEPTION(err,"GadgetImageMessageReader, failed to allocate memory\n");
            imgh->release();
            return 0;
        }

        imgh->cont(data);

        if ((recv_count = stream->recv_n(data->getObjectPtr()->get_data_ptr(), sizeof(T)*data->getObjectPtr()->get_number_of_elements())) <= 0) {
	  GERROR("GadgetImageMessageReader, failed to read data from socket\n");
	  imgh->release();
	  return 0;
        }

        return imgh;
    }
};

// for images with attributes
template <typename T> class GadgetImageAttribMessageReader : public GadgetMessageReader
{
public:

    typedef unsigned long long size_t_type;

    virtual ACE_Message_Block* read(ACE_SOCK_Stream* stream) 
    {
        GadgetContainerMessage<ISMRMRD::ImageHeader>* imgh = 
            new GadgetContainerMessage<ISMRMRD::ImageHeader>();

        GadgetContainerMessage<ISMRMRD::MetaContainer>* imgAttrib = 
            new GadgetContainerMessage<ISMRMRD::MetaContainer>();

        // read in ISMRMRD image header
        ssize_t recv_count = 0;
        if ((recv_count = stream->recv_n( imgh->getObjectPtr(), sizeof(ISMRMRD::ImageHeader))) <= 0)
        {
	  GERROR("GadgetImageAttribMessageReader, failed to read IMAGE Header\n");
	  imgh->release();
	  imgAttrib->release();
	  return 0;
        }

        // read in gadgetron image meta attributes
        size_t_type len(0);
        if ( ( recv_count = stream->recv_n( &len, sizeof(size_t_type)) ) <= 0 )
        {
	  GERROR("GadgetImageAttribMessageReader, failed to read IMAGE Meta Attributes length\n");
	  imgh->release();
	  imgAttrib->release();
	  return 0;
        }

        char* buf = NULL;
        try
        {
            buf = new char[len];
            if ( buf == NULL )
            {
	      GERROR("GadgetImageAttribMessageReader, failed to allocate IMAGE Meta Attributes buffer\n");
	      imgh->release();
	      imgAttrib->release();
	      return 0;
            }

            memset(buf, '\0', len);
            memcpy(buf, &len, sizeof(size_t_type));
        }
        catch (std::runtime_error &err)
        {
            GEXCEPTION(err,"GadgetImageAttribMessageReader, failed to allocate IMAGE Meta Attributes buffer\n");
            imgh->release();
            imgAttrib->release();
            return 0;
        }

        if ( ( recv_count = stream->recv_n( buf, len) ) <= 0 )
        {
	  GERROR("GadgetImageAttribMessageReader, failed to read IMAGE Meta Attributes\n");
	  imgh->release();
	  imgAttrib->release();
	  delete [] buf;
	  return 0;
        }

        try
        {
            ISMRMRD::deserialize(buf, *imgAttrib->getObjectPtr());
        }
        catch(...)
        {
	  GERROR("GadgetImageAttribMessageReader, failed to deserialize IMAGE Meta Attributes\n");
	  imgh->release();
	  imgAttrib->release();
	  delete [] buf;
	  return 0;
        }

        delete [] buf;

        // read in image content
        std::vector<size_t> dims(3);
        dims[0] = imgh->getObjectPtr()->matrix_size[0];
        dims[1] = imgh->getObjectPtr()->matrix_size[1];
        dims[2] = imgh->getObjectPtr()->matrix_size[2];

        if (imgh->getObjectPtr()->channels > 1)
        {
            dims.push_back(imgh->getObjectPtr()->channels);
        }

        GadgetContainerMessage< hoNDArray< T > >* data = new GadgetContainerMessage< hoNDArray< T > >();

        try
        {
            data->getObjectPtr()->create(&dims);
        }
        catch (std::runtime_error &err)
        {
            GEXCEPTION(err,"GadgetImageAttribMessageReader, failed to allocate memory\n");
            imgh->release();
            imgAttrib->release();
            data->release();
            return 0;
        }

        imgh->cont(data);
        data->cont(imgAttrib);

        if ((recv_count = stream->recv_n(data->getObjectPtr()->get_data_ptr(), sizeof(T)*data->getObjectPtr()->get_number_of_elements())) <= 0)
        {
	  GERROR("GadgetImageAttribMessageReader, failed to read data from socket\n");
	  imgh->release();
	  imgAttrib->release();
	  data->release();
	  return 0;
        }

        return imgh;
    }
};

}

#endif //GADGETSOCKETRECEIVER_H
