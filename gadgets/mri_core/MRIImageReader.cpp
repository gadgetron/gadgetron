#include "MRIImageReader.h"
#include <ismrmrd/ismrmrd.h>
#include <ismrmrd/meta.h>

/*


class GadgetronClientImageMessageReader : public GadgetronClientMessageReader
{

public:
    GadgetronClientImageMessageReader(std::string filename, std::string groupname)
        : file_name_(filename)
        , group_name_(groupname)
    {

    }

    ~GadgetronClientImageMessageReader() {
    } 

    template <typename T> 
    void read_data_attrib(tcp::socket* stream, const ISMRMRD::ImageHeader& h, ISMRMRD::Image<T>& im)
    {
        im.setHead(h);

        typedef unsigned long long size_t_type;

        //Read meta attributes
        size_t_type meta_attrib_length;
        boost::asio::read(*stream, boost::asio::buffer(&meta_attrib_length, sizeof(size_t_type)));

        if (meta_attrib_length>0)
        {
            std::string meta_attrib(meta_attrib_length, 0);
            boost::asio::read(*stream, boost::asio::buffer(const_cast<char*>(meta_attrib.c_str()), meta_attrib_length));
            im.setAttributeString(meta_attrib);
        }

        //Read image data
        boost::asio::read(*stream, boost::asio::buffer(im.getDataPtr(), im.getDataSize()));
        {
            if (!dataset_) {

                {
                    mtx.lock();
                    dataset_ = boost::shared_ptr<ISMRMRD::Dataset>(new ISMRMRD::Dataset(file_name_.c_str(), group_name_.c_str(), true)); // create if necessary 
                    mtx.unlock();
                }
            }

            std::stringstream st1;
            st1 << "image_" << h.image_series_index;
            std::string image_varname = st1.str();

            {
                mtx.lock();
                //TODO should this be wrapped in a try/catch?
                dataset_->appendImage(image_varname, im);
                mtx.unlock();
            }
        }
    }

    virtual void read(tcp::socket* stream) 
    {
        //Read the image headerfrom the socket
        ISMRMRD::ImageHeader h;
        boost::asio::read(*stream, boost::asio::buffer(&h,sizeof(ISMRMRD::ImageHeader)));

        if (h.data_type == ISMRMRD::ISMRMRD_USHORT)
        {
            ISMRMRD::Image<unsigned short> im;
            this->read_data_attrib(stream, h, im);
        }
        else if (h.data_type == ISMRMRD::ISMRMRD_SHORT)
        {
            ISMRMRD::Image<short> im;
            this->read_data_attrib(stream, h, im);
        }
        else if (h.data_type == ISMRMRD::ISMRMRD_UINT)
        {
            ISMRMRD::Image<unsigned int> im;
            this->read_data_attrib(stream, h, im);
        }
        else if (h.data_type == ISMRMRD::ISMRMRD_INT)
        {
            ISMRMRD::Image<int> im;
            this->read_data_attrib(stream, h, im);
        }
        else if (h.data_type == ISMRMRD::ISMRMRD_FLOAT)
        {
            ISMRMRD::Image<float> im;
            this->read_data_attrib(stream, h, im);
        }
        else if (h.data_type == ISMRMRD::ISMRMRD_DOUBLE)
        {
            ISMRMRD::Image<double> im;
            this->read_data_attrib(stream, h, im);
        }
        else if (h.data_type == ISMRMRD::ISMRMRD_CXFLOAT)
        {
            ISMRMRD::Image< std::complex<float> > im;
            this->read_data_attrib(stream, h, im);
        }
        else if (h.data_type == ISMRMRD::ISMRMRD_CXDOUBLE)
        {
            ISMRMRD::Image< std::complex<double> > im;
            this->read_data_attrib(stream, h, im);
        }
        else
        {
            throw GadgetronClientException("Invalide image data type ... ");
        }
    }



 */

namespace Gadgetron{

  ACE_Message_Block* MRIImageReader::read(ACE_SOCK_Stream* stream)
  {

    auto h = new GadgetContainerMessage< ISMRMRD::ImageHeader >();
    
    size_t recv_count = 0;
    if ((recv_count = stream->recv_n(h->getObjectPtr(), sizeof(ISMRMRD::ImageHeader))) <= 0) {
      h->release();
      GERROR("Failed to read ISMRMRD Image Header\n");
      return 0;
    }

    typedef unsigned long long size_t_type;
    
    //Read meta attributes
    size_t_type meta_attrib_length;
    if ((recv_count = stream->recv_n(&meta_attrib_length, sizeof(size_t_type))) <= 0) {
      h->release();
      GERROR("Failed to read length of meta attributes\n");
      return 0;
    }
    
    
    GadgetContainerMessage<ISMRMRD::MetaContainer>* meta = 0;
    if (meta_attrib_length>0)
    {
      char* buffer = new char[meta_attrib_length+1];
      if ((recv_count = stream->recv_n(buffer, meta_attrib_length)) <= 0) {
	h->release();
	GERROR("Failed to read meta attributes\n");
	return 0;
      }

      meta = new GadgetContainerMessage<ISMRMRD::MetaContainer>();

      ISMRMRD::deserialize(buffer, *(meta->getObjectPtr()));
      
      delete [] buffer;
    } 

    //Read actual image data
    char* data_ptr;
    size_t data_size;
    
    size_t elements =
      h->getObjectPtr()->matrix_size[0] *
      h->getObjectPtr()->matrix_size[1] *
      h->getObjectPtr()->matrix_size[2] *
      h->getObjectPtr()->channels;

    std::vector<size_t> img_dims{h->getObjectPtr()->matrix_size[0],
	h->getObjectPtr()->matrix_size[1], h->getObjectPtr()->matrix_size[2], h->getObjectPtr()->channels};

    ACE_Message_Block* data = 0;
    try {
      if (h->getObjectPtr()->data_type == ISMRMRD::ISMRMRD_USHORT)
      {
	auto d = new GadgetContainerMessage< hoNDArray<uint16_t> >();
	d->getObjectPtr()->create(img_dims);
	data_ptr = reinterpret_cast<char*>(d->getObjectPtr()->get_data_ptr());
	data_size = elements * sizeof(uint16_t);
	h->cont(d);
      }
      else if (h->getObjectPtr()->data_type == ISMRMRD::ISMRMRD_SHORT)
      {
	auto d = new GadgetContainerMessage< hoNDArray<int16_t> >();
	d->getObjectPtr()->create(img_dims);
	data_ptr = reinterpret_cast<char*>(d->getObjectPtr()->get_data_ptr());
	data_size = elements * sizeof(int16_t);
	h->cont(d);
      }
      else if (h->getObjectPtr()->data_type == ISMRMRD::ISMRMRD_UINT)
      {
	auto d = new GadgetContainerMessage< hoNDArray<uint32_t> >();
	d->getObjectPtr()->create(img_dims);
	data_ptr = reinterpret_cast<char*>(d->getObjectPtr()->get_data_ptr());
	data_size = elements * sizeof(uint32_t);
	h->cont(d);
      }
      else if (h->getObjectPtr()->data_type == ISMRMRD::ISMRMRD_INT)
      {
	auto d = new GadgetContainerMessage< hoNDArray<int32_t> >();
	d->getObjectPtr()->create(img_dims);
	data_ptr = reinterpret_cast<char*>(d->getObjectPtr()->get_data_ptr());
	data_size = elements * sizeof(int32_t);
	h->cont(d);
      }
      else if (h->getObjectPtr()->data_type == ISMRMRD::ISMRMRD_FLOAT)
      {
	auto d = new GadgetContainerMessage< hoNDArray<float> >();
	d->getObjectPtr()->create(img_dims);
	data_ptr = reinterpret_cast<char*>(d->getObjectPtr()->get_data_ptr());
	data_size = elements * sizeof(float);
	h->cont(d);
      }
      else if (h->getObjectPtr()->data_type == ISMRMRD::ISMRMRD_DOUBLE)
      {
	auto d = new GadgetContainerMessage< hoNDArray<double> >();
	d->getObjectPtr()->create(img_dims);
	data_ptr = reinterpret_cast<char*>(d->getObjectPtr()->get_data_ptr());
	data_size = elements * sizeof(double);
	h->cont(d);
      }
      else if (h->getObjectPtr()->data_type == ISMRMRD::ISMRMRD_CXFLOAT)
      {
	auto d = new GadgetContainerMessage< hoNDArray< std::complex<float> > >();
	d->getObjectPtr()->create(img_dims);
	data_ptr = reinterpret_cast<char*>(d->getObjectPtr()->get_data_ptr());
	data_size = elements * sizeof(std::complex<float>);
	h->cont(d);

      }
      else if (h->getObjectPtr()->data_type == ISMRMRD::ISMRMRD_CXDOUBLE)
      {
	auto d = new GadgetContainerMessage< hoNDArray< std::complex<double> > >();
	d->getObjectPtr()->create(img_dims);
	data_ptr = reinterpret_cast<char*>(d->getObjectPtr()->get_data_ptr());
	data_size = elements * sizeof(std::complex<double>);
	h->cont(d);
      }
      else
      {
	throw std::runtime_error("Unknow data type");
      }
    } catch (std::runtime_error &err) {
      GEXCEPTION(err, "Unable to create image array\n");
      if (h) h->release();
      if (meta) meta->release();
    }

    //Attach the meta data if it exists
    h->cont()->cont(meta);
    
    //Now lets read the data from the socket into the array
    if ((recv_count = stream->recv_n(data_ptr, data_size)) <= 0) {
      h->release();
      GERROR("Failed to read image array data\n");
      return 0;
    }
    
    return h;
  }

  
GADGETRON_READER_FACTORY_DECLARE(MRIImageReader)

}
