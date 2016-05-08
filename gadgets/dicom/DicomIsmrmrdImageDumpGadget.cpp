
#include "DicomIsmrmrdImageDumpGadget.h"
#include <iomanip>
#include <boost/filesystem.hpp>

namespace bf = boost::filesystem;

namespace Gadgetron {

    DicomIsmrmrdImageDumpGadget::DicomIsmrmrdImageDumpGadget(): BaseClass(), dcmFile()
    {
    }

    DicomIsmrmrdImageDumpGadget::~DicomIsmrmrdImageDumpGadget()
    {
    }

    int DicomIsmrmrdImageDumpGadget::process_config(ACE_Message_Block* mb)
    {
        ISMRMRD::IsmrmrdHeader h;
        deserialize(mb->rd_ptr(), h);

        xml = h;
        Gadgetron::fill_dicom_image_from_ismrmrd_header(h, dcmFile);

        // get measurement id
        std::string measurement_id = "";
        std::string ismrmrd_filename = "";

        if ( h.measurementInformation )
        {
            if ( h.measurementInformation->measurementID )
            {
                measurement_id = *h.measurementInformation->measurementID;
            }
        }

        GDEBUG("Measurement ID: %s\n", measurement_id.c_str());

        // analyze measurement id
        if(measurement_id.size()>0)
        {
            std::string mid = measurement_id;
            size_t ind = mid.find("_");
            if(ind!=std::string::npos)
            {
                device_ = mid.substr(0, ind);
                mid = mid.substr(ind+1, std::string::npos);

                ind = mid.find("_");
                if(ind!=std::string::npos)
                {
                    patient_ = mid.substr(0, ind);
                    mid = mid.substr(ind+1, std::string::npos);

                    ind = mid.find("_");
                    if(ind!=std::string::npos)
                    {
                        study_ = mid.substr(0, ind);
                        measurement_ = mid.substr(ind+1, std::string::npos);
                    }
                }
            }
        }

        // get the protocol name
        if ( h.measurementInformation )
        {
            if ( h.measurementInformation->protocolName )
            {
                protocol_name_ = *h.measurementInformation->protocolName;
            }
        }

        // if patient info was not set, assemble them
        if(!h.subjectInformation.is_present() || !h.subjectInformation.get().patientName.is_present())
        {
            if(h.acquisitionSystemInformation.is_present() && h.acquisitionSystemInformation.get().institutionName.is_present())
            {
                if ( file_prefix.value().empty() )
                {
                    patient_string_ = h.acquisitionSystemInformation.get().institutionName.get() + "_" + device_ + "_" + patient_;
                }
                else
                {
                    patient_string_ = h.acquisitionSystemInformation.get().institutionName.get() + "_" + device_ + "_" + file_prefix.value() + "_" + patient_;
                }
            }
        }

        std::string folder_prefix;

        if ( file_prefix.value().empty() )
        {
            folder_prefix = "ISMRMRD_DICOM_DUMP";
        }
        else
        {
            folder_prefix = file_prefix.value();
        }

        if (append_id.value() && measurement_id.size())
        {
            folder_prefix.append("_");
            folder_prefix.append(measurement_id);
        }

        dicom_folder_ = folder.value() + "/" + folder_prefix;
        bf::path p(dicom_folder_);

        if (!exists(p))
        {
            try
            {
                bf::create_directory(p);
            }
            catch(...)
            {
                GERROR("Error caught trying to create folder %s\n", folder.value().c_str());
                return GADGET_FAIL;
            }
        }
        else
        {
            if (!is_directory(p))
            {
                GERROR("Specified path is not a directory\n");
                return GADGET_FAIL;
            }
        }

        GDEBUG_STREAM("Dicom dump folder : " << dicom_folder_);

        return GADGET_OK;
    }

    int DicomIsmrmrdImageDumpGadget::process(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1)
    {
        if (!this->controller_) {
            GERROR("Cannot return result to controller, no controller set");
            return GADGET_FAIL;
        }

        // --------------------------------------------------
        ISMRMRD::ImageHeader *img = m1->getObjectPtr();

        uint16_t data_type = img->data_type;

        if (data_type == ISMRMRD::ISMRMRD_USHORT)
        {
            GadgetContainerMessage< hoNDArray< unsigned short > >* datamb = AsContainerMessage< hoNDArray< unsigned short > >(m1->cont());
            if (!datamb)
            {
                GERROR_STREAM("DicomIsmrmrdImageDumpGadget::process, invalid image message objects");
                return GADGET_FAIL;
            }

            if (this->dump_image(m1, datamb) != GADGET_OK)
            {
                GERROR_STREAM("DicomIsmrmrdImageDumpGadget::dump_image failed for unsigned short ... ");
            }
        }
        else if (data_type == ISMRMRD::ISMRMRD_SHORT)
        {
            GadgetContainerMessage< hoNDArray< short > >* datamb = AsContainerMessage< hoNDArray< short > >(m1->cont());
            if (!datamb)
            {
                GERROR_STREAM("DicomIsmrmrdImageDumpGadget::process, invalid image message objects");
                return GADGET_FAIL;
            }

            if (this->dump_image(m1, datamb) != GADGET_OK)
            {
                GERROR_STREAM("DicomIsmrmrdImageDumpGadget::dump_image failed for short ... ");
            }
        }
        else if (data_type == ISMRMRD::ISMRMRD_UINT)
        {
            GadgetContainerMessage< hoNDArray< unsigned int > >* datamb = AsContainerMessage< hoNDArray< unsigned int > >(m1->cont());
            if (!datamb)
            {
                GERROR_STREAM("DicomIsmrmrdImageDumpGadget::process, invalid image message objects");
                return GADGET_FAIL;
            }

            if (this->dump_image(m1, datamb) != GADGET_OK)
            {
                GERROR_STREAM("DicomIsmrmrdImageDumpGadget::dump_image failed for unsigned int ... ");
            }
        }
        else if (data_type == ISMRMRD::ISMRMRD_INT)
        {
            GadgetContainerMessage< hoNDArray< int > >* datamb = AsContainerMessage< hoNDArray< int > >(m1->cont());
            if (!datamb)
            {
                GERROR_STREAM("DicomIsmrmrdImageDumpGadget::process, invalid image message objects");
                return GADGET_FAIL;
            }

            if (this->dump_image(m1, datamb) != GADGET_OK)
            {
                GERROR_STREAM("DicomIsmrmrdImageDumpGadget::dump_image failed for int ... ");
            }
        }
        else if (data_type == ISMRMRD::ISMRMRD_FLOAT)
        {
            GadgetContainerMessage< hoNDArray< float > >* datamb = AsContainerMessage< hoNDArray< float > >(m1->cont());
            if (!datamb)
            {
                GERROR_STREAM("DicomIsmrmrdImageDumpGadget::process, invalid image message objects");
                return GADGET_FAIL;
            }

            if (this->dump_image(m1, datamb) != GADGET_OK)
            {
                GERROR_STREAM("DicomIsmrmrdImageDumpGadget::dump_image failed for float ... ");
            }
        }
        else if (data_type == ISMRMRD::ISMRMRD_DOUBLE)
        {
            GadgetContainerMessage< hoNDArray< double > >* datamb = AsContainerMessage< hoNDArray< double > >(m1->cont());
            if (!datamb)
            {
                GERROR_STREAM("DicomIsmrmrdImageDumpGadget::process, invalid image message objects");
                return GADGET_FAIL;
            }

            if (this->dump_image(m1, datamb) != GADGET_OK)
            {
                GERROR_STREAM("DicomIsmrmrdImageDumpGadget::dump_image failed for double ... ");
            }
        }
        else if (data_type == ISMRMRD::ISMRMRD_CXFLOAT)
        {
            GERROR_STREAM("DicomIsmrmrdImageDumpGadget::process, does not supprot ISMRMRD_CXFLOAT data type");
        }
        else if (data_type == ISMRMRD::ISMRMRD_CXDOUBLE)
        {
            GERROR_STREAM("DicomIsmrmrdImageDumpGadget::process, does not supprot ISMRMRD_CXDOUBLE data type");
        }
        else
        {
            GERROR_STREAM("DicomIsmrmrdImageDumpGadget::process, does not supprot data type : " << data_type);
        }

        if (this->next()->putq(m1) == -1)
        {
            m1->release();
            GDEBUG_STREAM("Unable to put ismrmrd image on next gadgets queue");
            return GADGET_FAIL;
        }

        return GADGET_OK;
    }

    std::string DicomIsmrmrdImageDumpGadget::get_date_string()
    {
        time_t rawtime;
        struct tm * timeinfo;
        time ( &rawtime );
        timeinfo = localtime ( &rawtime );

        std::stringstream str;
        str << timeinfo->tm_year+1900
            << std::setw(2) << std::setfill('0') << timeinfo->tm_mon+1
            << std::setw(2) << std::setfill('0') << timeinfo->tm_mday;

        std::string ret = str.str();

        return ret;
    }

    std::string DicomIsmrmrdImageDumpGadget::get_time_string()
    {
        time_t rawtime;
        struct tm * timeinfo;
        time ( &rawtime );
        timeinfo = localtime ( &rawtime );

        std::stringstream str;
        str << std::setw(2) << std::setfill('0') << timeinfo->tm_hour
            << std::setw(2) << std::setfill('0') << timeinfo->tm_min
            << std::setw(2) << std::setfill('0') << timeinfo->tm_sec;

        std::string ret = str.str();

        return ret;
    }

    GADGET_FACTORY_DECLARE(DicomIsmrmrdImageDumpGadget)

} /* namespace Gadgetron */
