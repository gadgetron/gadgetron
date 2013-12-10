#ifndef IMAGEWRITER_H
#define IMAGEWRITER_H

#include <fstream>

#include "GadgetImageMessageReader.h"

namespace Gadgetron
{

template <typename T> class ImageWriter : public GadgetImageMessageReader<T>
{

public:
	ImageWriter()
	: number_of_calls_(0)
	{}

	virtual ~ImageWriter() {};

	virtual ACE_Message_Block* read(ACE_SOCK_Stream* socket) 
	{
		// Invoke parent's read
		ACE_Message_Block* mb = GadgetImageMessageReader<T>::read(socket);

		if (!mb) {
			GADGET_DEBUG1("Read failed in parent\n");
			return 0;
		}

		GadgetContainerMessage<ISMRMRD::ImageHeader> * img_head_mb =
				dynamic_cast<GadgetContainerMessage<ISMRMRD::ImageHeader> *>(mb);

		if (!img_head_mb) {
			GADGET_DEBUG1("Failed in dynamic cast\n");
			mb->release();
			return 0;
		}

		//GADGET_DEBUG2("Received image with %d channels\n", img_head_mb->getObjectPtr()->channels);

		GadgetContainerMessage<hoNDArray< T > > * img_data_mb =
				dynamic_cast<GadgetContainerMessage<hoNDArray< T > > *>(img_head_mb->cont());

		if (!img_data_mb) {
			GADGET_DEBUG1("Failed in dynamic cast\n");
			mb->release();
			return 0;
		}

		if (this->process_image(img_head_mb->getObjectPtr(), img_data_mb->getObjectPtr()) < 0) {
			GADGET_DEBUG1("Failed to process image\n");
			mb->release();
			return 0;
		}

		return mb;
	}

	virtual int process_image(ISMRMRD::ImageHeader* img_head,
			hoNDArray< T >* data)
	{
		ACE_DEBUG( (LM_DEBUG, ACE_TEXT("Image Writer writing image\n")) );

		char filename[1024];

		switch (sizeof(T)) {

		case (8): //Complex float
    			sprintf(filename, "out_%05d.cplx", (int)number_of_calls_);
		break;
		case (4): //Real floats
				sprintf(filename, "out_%05d.real", (int)number_of_calls_);
		break;
		case (2): //Unsigned short
				sprintf(filename, "out_%05d.short", (int)number_of_calls_);
		break;
		default:
			sprintf(filename, "out_%05d.cplx", (int)number_of_calls_);
			break;
		}

		std::ofstream outfile;
		outfile.open (filename, std::ios::out|std::ios::binary);

		if (outfile.good()) {
			int ndim = 4;
			int dims[4];
			size_t elements = 1;
			dims[0] = img_head->matrix_size[0]; elements*=dims[0];
			dims[1] = img_head->matrix_size[1]; elements*=dims[1];
			dims[2] = img_head->matrix_size[2]; elements*=dims[2];
			dims[3] = img_head->channels; elements*=dims[3];

			outfile.write((char*)&ndim,sizeof(int));
			outfile.write((char*)dims,sizeof(int)*4);
			outfile.write((char*)data->get_data_ptr(),sizeof(T)*elements);
			outfile.close();
			number_of_calls_++;
		} else {
			GADGET_DEBUG1("File is not good for writing\n");
			return GADGET_FAIL;
		}

		return GADGET_OK;
	}

protected:
	size_t number_of_calls_;
};

}
#endif //IMAGE_WRITER
