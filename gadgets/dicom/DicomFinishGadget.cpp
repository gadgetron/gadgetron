#include <vector>
#include "boost/date_time/gregorian/gregorian.hpp"

#include "DicomFinishGadget.h"
#include "ismrmrd/xml.h"

namespace Gadgetron {

    int DicomFinishGadget::process_config(ACE_Message_Block* mb)
    {
        ISMRMRD::IsmrmrdHeader h;
        deserialize(mb->rd_ptr(), h);

        xml = h;

        Gadgetron::fill_dicom_image_from_ismrmrd_header(h, dcmFile);

        ISMRMRD::MeasurementInformation meas_info = *h.measurementInformation;

        if (meas_info.seriesInstanceUIDRoot) {
            seriesIUIDRoot = *meas_info.seriesInstanceUIDRoot;
        }

        if (meas_info.initialSeriesNumber) {
            this->initialSeriesNumber = (long)*meas_info.initialSeriesNumber;
        }
        else {
            this->initialSeriesNumber = 0;
        }

        return GADGET_OK;
    }

    int DicomFinishGadget::process(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1)
    {


        // --------------------------------------------------
        ISMRMRD::ImageHeader *img = m1->getObjectPtr();

        uint16_t data_type = img->data_type;

        if (data_type == ISMRMRD::ISMRMRD_USHORT)
        {
            GadgetContainerMessage< hoNDArray< unsigned short > >* datamb = AsContainerMessage< hoNDArray< unsigned short > >(m1->cont());
            if (!datamb)
            {
                GERROR("DicomFinishGadget::process, invalid image message objects\n");
                return GADGET_FAIL;
            }

            if (this->write_data_attrib(m1, datamb) != GADGET_OK)
            {
                GERROR("DicomFinishGadget::write_data_attrib failed for unsigned short ... \n");
                return GADGET_FAIL;
            }
        }
        else if (data_type == ISMRMRD::ISMRMRD_SHORT)
        {
            GadgetContainerMessage< hoNDArray< short > >* datamb = AsContainerMessage< hoNDArray< short > >(m1->cont());
            if (!datamb)
            {
                GERROR("DicomFinishGadget::process, invalid image message objects\n");
                return GADGET_FAIL;
            }

            if (this->write_data_attrib(m1, datamb) != GADGET_OK)
            {
                GERROR("DicomFinishGadget::write_data_attrib failed for short ... \n");
                return GADGET_FAIL;
            }
        }
        else if (data_type == ISMRMRD::ISMRMRD_UINT)
        {
            GadgetContainerMessage< hoNDArray< unsigned int > >* datamb = AsContainerMessage< hoNDArray< unsigned int > >(m1->cont());
            if (!datamb)
            {
                GERROR("DicomFinishGadget::process, invalid image message objects\n");
                return GADGET_FAIL;
            }

            if (this->write_data_attrib(m1, datamb) != GADGET_OK)
            {
                GERROR("DicomFinishGadget::write_data_attrib failed for unsigned int ... \n");
                return GADGET_FAIL;
            }
        }
        else if (data_type == ISMRMRD::ISMRMRD_INT)
        {
            GadgetContainerMessage< hoNDArray< int > >* datamb = AsContainerMessage< hoNDArray< int > >(m1->cont());
            if (!datamb)
            {
                GERROR("DicomFinishGadget::process, invalid image message objects\n");
                return GADGET_FAIL;
            }

            if (this->write_data_attrib(m1, datamb) != GADGET_OK)
            {
                GERROR("DicomFinishGadget::write_data_attrib failed for int ... \n");
                return GADGET_FAIL;
            }
        }
        else if (data_type == ISMRMRD::ISMRMRD_FLOAT)
        {
            GadgetContainerMessage< hoNDArray< float > >* datamb = AsContainerMessage< hoNDArray< float > >(m1->cont());
            if (!datamb)
            {
                GERROR("DicomFinishGadget::process, invalid image message objects\n");
                return GADGET_FAIL;
            }

            if (this->write_data_attrib(m1, datamb) != GADGET_OK)
            {
                GERROR("DicomFinishGadget::write_data_attrib failed for float ... \n");
                return GADGET_FAIL;
            }
        }
        else if (data_type == ISMRMRD::ISMRMRD_DOUBLE)
        {
            GadgetContainerMessage< hoNDArray< double > >* datamb = AsContainerMessage< hoNDArray< double > >(m1->cont());
            if (!datamb)
            {
                GERROR("DicomFinishGadget::process, invalid image message objects\n");
                return GADGET_FAIL;
            }

            if (this->write_data_attrib(m1, datamb) != GADGET_OK)
            {
                GERROR("DicomFinishGadget::write_data_attrib failed for double ... \n");
                return GADGET_FAIL;
            }
        }
        else if (data_type == ISMRMRD::ISMRMRD_CXFLOAT)
        {
            GERROR("DicomFinishGadget::process, does not supprot ISMRMRD_CXFLOAT data type\n");
            return GADGET_FAIL;
        }
        else if (data_type == ISMRMRD::ISMRMRD_CXDOUBLE)
        {
            GERROR("DicomFinishGadget::process, does not supprot ISMRMRD_CXDOUBLE data type\n");
            return GADGET_FAIL;
        }

        return GADGET_OK;
    }

    GADGET_FACTORY_DECLARE(DicomFinishGadget)

} /* namespace Gadgetron */
