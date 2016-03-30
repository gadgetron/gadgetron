#include "GadgetIsmrmrdReadWrite.h"
#include "IsmrmrdDumpGadget.h"
#include <iomanip>
#include <boost/filesystem.hpp>
#include <ismrmrd/ismrmrd.h>
#include <ismrmrd/xml.h>

namespace bf = boost::filesystem;

namespace Gadgetron
{
    std::string get_date_time_string()
    {
        time_t rawtime;
        struct tm * timeinfo;
        time ( &rawtime );
        timeinfo = localtime ( &rawtime );

        std::stringstream str;
        str << timeinfo->tm_year+1900
            << std::setw(2) << std::setfill('0') << timeinfo->tm_mon+1
            << std::setw(2) << std::setfill('0') << timeinfo->tm_mday
            << "-"
            << std::setw(2) << std::setfill('0') << timeinfo->tm_hour
            << std::setw(2) << std::setfill('0') << timeinfo->tm_min
            << std::setw(2) << std::setfill('0') << timeinfo->tm_sec;

        std::string ret = str.str();

        return ret;
    }

    IsmrmrdDumpGadget::IsmrmrdDumpGadget()
                    : Gadget2<ISMRMRD::AcquisitionHeader,hoNDArray< std::complex<float> > >()
    {
    }

    int IsmrmrdDumpGadget::process_config(ACE_Message_Block* mb)
    {
        
        ISMRMRD::IsmrmrdHeader ismrmrd_header;
        ISMRMRD::deserialize(mb->rd_ptr(), ismrmrd_header);

        std::string measurement_id = "";
        std::string ismrmrd_filename = "";
        
        if ( ismrmrd_header.measurementInformation ) {
            if ( ismrmrd_header.measurementInformation->measurementID ) {
                measurement_id = *ismrmrd_header.measurementInformation->measurementID;
            }
        }

        GDEBUG("Measurement ID: %s\n", measurement_id.c_str());
        
        bf::path p(folder.value());

        if (!exists(p)) {
            try {
                bf::create_directory(p);
            } catch(...) {
                GERROR("Error caught trying to create folder %s\n", folder.value().c_str());
                return GADGET_FAIL;
            }
        } else {
            if (!is_directory(p)) {
                GERROR("Specified path is not a directory\n");
                return GADGET_FAIL;
            }
        }
        
        if ( file_prefix.value().empty() )
        {
            ismrmrd_filename = "ISMRMRD_DUMP";
        } else {
            ismrmrd_filename = file_prefix.value();
        }

        if (append_id.value() && measurement_id.size()) {
            ismrmrd_filename.append("_");
            ismrmrd_filename.append(measurement_id);
        }
        
        //Generate filename
        if (append_timestamp.value())
        {
            ismrmrd_filename.append("_");
            ismrmrd_filename.append(get_date_time_string());
        }

        ismrmrd_filename.append(".h5");

        p /= ismrmrd_filename;

        ismrmrd_filename = p.string();
        ismrmrd_dataset_ = boost::shared_ptr<ISMRMRD::Dataset>(new ISMRMRD::Dataset(ismrmrd_filename.c_str(), "dataset"));

        std::string xml_config(mb->rd_ptr());

        try {
            ismrmrd_dataset_->writeHeader(xml_config);
        }
        catch (...)
        {
            GDEBUG("Failed to write XML header to HDF file\n");
            return GADGET_FAIL;
        }

        return GADGET_OK;
    }

    int IsmrmrdDumpGadget::process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1, GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2)
    {
        ISMRMRD::Acquisition ismrmrd_acq;

        ismrmrd_acq.setHead(*m1->getObjectPtr());

        memcpy((void *)ismrmrd_acq.getDataPtr(), m2->getObjectPtr()->get_data_ptr(), 
               sizeof(float)*m2->getObjectPtr()->get_number_of_elements()*2);

        if (m2->cont())
        {
            //Write trajectory
            if (ismrmrd_acq.trajectory_dimensions() == 0)
            {
                GDEBUG("Malformed dataset. Trajectory attached but trajectory dimensions == 0\n");
                return GADGET_FAIL;
            }

            GadgetContainerMessage< hoNDArray<float> >* m3 = AsContainerMessage< hoNDArray<float> >(m2->cont());

            if (!m3)
            {
                GDEBUG("Error casting trajectory data package");
                return GADGET_FAIL;
            } 

            memcpy((void *)ismrmrd_acq.getTrajPtr(), m3->getObjectPtr()->get_data_ptr(),
		   sizeof(float)*m3->getObjectPtr()->get_number_of_elements());

        }
        else
        {
            if (ismrmrd_acq.trajectory_dimensions() != 0)
            {
                GDEBUG("Malformed dataset. Trajectory dimensions not zero but no trajectory attached\n");
                return GADGET_FAIL;
            }
        }

        {
            try {
                ismrmrd_dataset_->appendAcquisition(ismrmrd_acq);
            }
            catch (...)
            {
                GDEBUG("Error appending ISMRMRD Dataset\n");
                return GADGET_FAIL;
            }
        }

        //It is enough to put the first one, since they are linked
        if (this->next()->putq(m1) == -1)
        {
            m1->release();
	    GERROR("IsmrmrdDumpGadget::process, passing data on to next gadget");
	    return -1;
        }

        return 0;
    }

    GADGET_FACTORY_DECLARE(IsmrmrdDumpGadget)
}
