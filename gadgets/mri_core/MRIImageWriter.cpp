#include "GadgetIsmrmrdReadWrite.h"
#include "MRIImageWriter.h"
#include "GadgetContainerMessage.h"
#include "hoNDArray.h"

#include <complex>

namespace Gadgetron{

    int MRIImageWriter::write(ACE_SOCK_Stream* sock, ACE_Message_Block* mb)
    {
        GadgetContainerMessage<ISMRMRD::ImageHeader>* imagemb =
            AsContainerMessage<ISMRMRD::ImageHeader>(mb);

        if (!imagemb)
        {
            GERROR("MRIImageWriter::write, invalid image message objects, 1\n");
            return -1;
        }

        uint16_t data_type = imagemb->getObjectPtr()->data_type;

        if (data_type == ISMRMRD::ISMRMRD_USHORT)
        {
            GadgetContainerMessage< hoNDArray< unsigned short > >* datamb = AsContainerMessage< hoNDArray< unsigned short > >(imagemb->cont());
            if (!datamb)
            {
                GERROR("MRIImageWriter::write, invalid image message objects\n");
                return -1;
            }

            if (this->write_data_attrib(sock, imagemb, datamb) != 0)
            {
                GERROR("MRIImageWriter::write_data_attrib failed for unsigned short ... \n");
                return -1;
            }
        }
        else if (data_type == ISMRMRD::ISMRMRD_SHORT)
        {
            GadgetContainerMessage< hoNDArray< short > >* datamb = AsContainerMessage< hoNDArray< short > >(imagemb->cont());
            if (!datamb)
            {
                GERROR("MRIImageWriter::write, invalid image message objects\n");
                return -1;
            }

            if (this->write_data_attrib(sock, imagemb, datamb) != 0)
            {
                GERROR("MRIImageWriter::write_data_attrib failed for short ... \n");
                return -1;
            }
        }
        else if (data_type == ISMRMRD::ISMRMRD_UINT)
        {
            GadgetContainerMessage< hoNDArray< unsigned int > >* datamb = AsContainerMessage< hoNDArray< unsigned int > >(imagemb->cont());
            if (!datamb)
            {
                GERROR("MRIImageWriter::write, invalid image message objects\n");
                return -1;
            }

            if (this->write_data_attrib(sock, imagemb, datamb) != 0)
            {
                GERROR("MRIImageWriter::write_data_attrib failed for unsigned int ... \n");
                return -1;
            }
        }
        else if (data_type == ISMRMRD::ISMRMRD_INT)
        {
            GadgetContainerMessage< hoNDArray< int > >* datamb = AsContainerMessage< hoNDArray< int > >(imagemb->cont());
            if (!datamb)
            {
                GERROR("MRIImageWriter::write, invalid image message objects\n");
                return -1;
            }

            if (this->write_data_attrib(sock, imagemb, datamb) != 0)
            {
                GERROR("MRIImageWriter::write_data_attrib failed for int ... \n");
                return -1;
            }
        }
        else if (data_type == ISMRMRD::ISMRMRD_FLOAT)
        {
            GadgetContainerMessage< hoNDArray< float > >* datamb = AsContainerMessage< hoNDArray< float > >(imagemb->cont());
            if (!datamb)
            {
                GERROR("MRIImageWriter::write, invalid image message objects\n");
                return -1;
            }

            if (this->write_data_attrib(sock, imagemb, datamb) != 0)
            {
                GERROR("MRIImageWriter::write_data_attrib failed for float ... \n");
                return -1;
            }
        }
        else if (data_type == ISMRMRD::ISMRMRD_DOUBLE)
        {
            GadgetContainerMessage< hoNDArray< double > >* datamb = AsContainerMessage< hoNDArray< double > >(imagemb->cont());
            if (!datamb)
            {
                GERROR("MRIImageWriter::write, invalid image message objects\n");
                return -1;
            }

            if (this->write_data_attrib(sock, imagemb, datamb) != 0)
            {
                GERROR("MRIImageWriter::write_data_attrib failed for double ... \n");
                return -1;
            }
        }
        else if (data_type == ISMRMRD::ISMRMRD_CXFLOAT)
        {
            GadgetContainerMessage< hoNDArray< std::complex<float> > >* datamb = AsContainerMessage< hoNDArray< std::complex<float> > >(imagemb->cont());
            if (!datamb)
            {
                GERROR("MRIImageWriter::write, invalid image message objects\n");
                return -1;
            }

            if (this->write_data_attrib(sock, imagemb, datamb) != 0)
            {
                GERROR("MRIImageWriter::write_data_attrib failed for std::complex<float> ... \n");
                return -1;
            }
        }
        else if (data_type == ISMRMRD::ISMRMRD_CXDOUBLE)
        {
            GadgetContainerMessage< hoNDArray< std::complex<double> > >* datamb = AsContainerMessage< hoNDArray< std::complex<double> > >(imagemb->cont());
            if (!datamb)
            {
                GERROR("MRIImageWriter::write, invalid image message objects\n");
                return -1;
            }

            if (this->write_data_attrib(sock, imagemb, datamb) != 0)
            {
                GERROR("MRIImageWriter::write_data_attrib failed for std::complex<double> ... \n");
                return -1;
            }
        }

        return 0;
    }

    GADGETRON_WRITER_FACTORY_DECLARE(MRIImageWriter)

}
