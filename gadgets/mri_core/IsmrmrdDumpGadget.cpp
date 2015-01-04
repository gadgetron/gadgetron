#include "GadgetIsmrmrdReadWrite.h"
#include "IsmrmrdDumpGadget.h"

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
                    , file_prefix_("ISMRMRD_DUMP")
                    , ismrmrd_file_name_("ISMRMRD_DUMP.h5") //This will be reset during configuration
                    , append_timestamp_(true)
    {
        file_prefix_ = "ISMRMRD_DUMP";
        append_timestamp_ = true;
    }

    int IsmrmrdDumpGadget::process_config(ACE_Message_Block* mb)
    {
        file_prefix_ = *(get_string_value("file_prefix").get());
        if ( file_prefix_.empty() )
        {
            file_prefix_ = "ISMRMRD_DUMP";
        }

        append_timestamp_ = get_bool_value("append_timestamp");

        //Generate filename
        if (append_timestamp_)
        {
            ismrmrd_file_name_ = file_prefix_ + std::string("_") + get_date_time_string() + std::string(".h5");
        }
        else
        {
            ismrmrd_file_name_ = file_prefix_ + std::string(".h5");
        }

        ismrmrd_dataset_ = boost::shared_ptr<ISMRMRD::Dataset>(new ISMRMRD::Dataset(ismrmrd_file_name_.c_str(), "dataset"));

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
