#include "IsmrmrdDumpGadget.h"
#include <iomanip>
#include <boost/filesystem.hpp>
#include "network_utils.h"

namespace bf = boost::filesystem;

namespace Gadgetron
{
    std::string get_date_time_string()
    {
        time_t rawtime;
        struct tm * timeinfo;
        time(&rawtime);
        timeinfo = localtime(&rawtime);

        std::stringstream str;
        str << timeinfo->tm_year + 1900
            << std::setw(2) << std::setfill('0') << timeinfo->tm_mon + 1
            << std::setw(2) << std::setfill('0') << timeinfo->tm_mday
            << "-"
            << std::setw(2) << std::setfill('0') << timeinfo->tm_hour
            << std::setw(2) << std::setfill('0') << timeinfo->tm_min
            << std::setw(2) << std::setfill('0') << timeinfo->tm_sec;

        std::string ret = str.str();

        return ret;
    }

    IsmrmrdDumpGadget::IsmrmrdDumpGadget() : BaseClass(), first_call_(true), save_ismrmrd_data_(true)
    {
    }

    IsmrmrdDumpGadget::~IsmrmrdDumpGadget()
    {
    }

    int IsmrmrdDumpGadget::process_config(ACE_Message_Block* mb)
    {

        ISMRMRD::deserialize(mb->rd_ptr(), ismrmrd_header_);
        ismrmrd_xml_ = std::string(mb->rd_ptr());

        // if ip_no_data_saving is set, check current ip of gadgetron server
        if (!this->ip_no_data_saving.value().empty())
        {
            std::vector<std::string> ip_list = this->ip_no_data_saving.value();

            GDEBUG_STREAM("IsmrmrdDumpGadget, find pre-set ip for no-data-saving : " << ip_list.size() << " [ ");

            size_t n;
            for (n = 0; n < ip_list.size(); n++)
            {
                GDEBUG_STREAM(ip_list[n]);
            }
            GDEBUG_STREAM(" ] ");

            // get current ip of gadgetron
            std::string host_name;
            std::vector<std::string> gt_ip_list;
            try
            {
                Gadgetron::find_gadgetron_ip(host_name, gt_ip_list);
            }
            catch (...)
            {
                GERROR_STREAM("Find gadgetrion ip failed ... ");
                gt_ip_list.clear();
            }

            GDEBUG_STREAM("IsmrmrdDumpGadget, find gadgetron host name : " << host_name);

            GDEBUG_STREAM("IsmrmrdDumpGadget, find gadgetron ip : " << gt_ip_list.size() << " [ ");
            for (n = 0; n < gt_ip_list.size(); n++)
            {
                GDEBUG_STREAM(gt_ip_list[n]);
            }
            GDEBUG_STREAM(" ] ");

            size_t m;

            for (n = 0; n < ip_list.size(); n++)
            {
                for (m = 0; m < gt_ip_list.size(); m++)
                {
                    if (ip_list[n] == gt_ip_list[m])
                    {
                        this->save_ismrmrd_data_ = false;
                        GDEBUG_STREAM("IsmrmrdDumpGadget, find matching ip pair : " << ip_list[n]);
                        break;
                    }
                }

                if (!this->save_ismrmrd_data_)
                {
                    break;
                }
            }
        }

        if (this->save_ismrmrd_data_)
        {
            boost::filesystem::path p(folder.value());

            if (!exists(p))
            {
                try
                {
                    boost::filesystem::create_directory(p);
                }
                catch (...)
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
        }

        return GADGET_OK;
    }

    int IsmrmrdDumpGadget::create_ismrmrd_dataset()
    {
        try {
            std::string measurement_id = "";
            std::string ismrmrd_filename = "";

            if (ismrmrd_header_.measurementInformation)
            {
                if (ismrmrd_header_.measurementInformation->measurementID)
                {
                    measurement_id = *ismrmrd_header_.measurementInformation->measurementID;
                }
            }

            GDEBUG("Measurement ID: %s\n", measurement_id.c_str());

            boost::filesystem::path p(folder.value());

            if (file_prefix.value().empty())
            {
                // try to use the protocol name
                if (ismrmrd_header_.measurementInformation.is_present())
                {
                    if (ismrmrd_header_.measurementInformation.get().protocolName.is_present())
                    {
                        ismrmrd_filename = ismrmrd_header_.measurementInformation.get().protocolName.get();

                        for (std::string::size_type i = 0; (i = ismrmrd_filename.find(" ", i)) != std::string::npos;)
                        {
                            ismrmrd_filename.replace(i, 1, "_");
                            i += 1;
                        }

                        GDEBUG("ismrmrd_filename: %s\n", ismrmrd_filename.c_str());
                    }
                    else
                    {
                        ismrmrd_filename = "ISMRMRD_DUMP";
                    }
                }
                else
                {
                    ismrmrd_filename = "ISMRMRD_DUMP";
                }
            }
            else
            {
                ismrmrd_filename = file_prefix.value();
            }

            ismrmrd_filename.append("_");
            ismrmrd_filename.append(measurement_id);

            std::string study_date, study_time;
            if (ismrmrd_header_.studyInformation)
            {
                if (ismrmrd_header_.studyInformation->studyDate)
                {
                    study_date = *ismrmrd_header_.studyInformation->studyDate;

                    std::string d(study_date);
                    d.erase(std::remove(d.begin(), d.end(), '-'), d.end());
                    study_date = d;
                }

                if (ismrmrd_header_.studyInformation->studyTime)
                {
                    study_time = *ismrmrd_header_.studyInformation->studyTime;

                    std::string d(study_time);
                    d.erase(std::remove(d.begin(), d.end(), ':'), d.end());
                    study_time = d;
                }
            }

            if (!study_date.empty() && !study_time.empty())
            {
                ismrmrd_filename.append("_");
                ismrmrd_filename.append(study_date);

                ismrmrd_filename.append("-");
                ismrmrd_filename.append(study_time);
            }
            else
            {
                ismrmrd_filename.append("_");
                ismrmrd_filename.append(get_date_time_string());
            }

            ismrmrd_filename.append(".h5");

            p /= ismrmrd_filename;

            ismrmrd_filename = p.string();
            GDEBUG_STREAM("KSpace dump file name : " << ismrmrd_filename);

            ismrmrd_dataset_ = boost::shared_ptr<ISMRMRD::Dataset>(new ISMRMRD::Dataset(ismrmrd_filename.c_str(), "dataset", true));
        }
        catch (...)
        {
            GERROR_STREAM("Error happened in IsmrmrdDumpGadget::create_ismrmrd_dataset(...) ... ");
            return GADGET_FAIL;
        }

        return GADGET_OK;
    }

    int IsmrmrdDumpGadget::process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1)
    {
        if (first_call_)
        {
            if (this->save_ismrmrd_data_)
            {
                this->create_ismrmrd_dataset();

                try
                {
                    ismrmrd_dataset_->writeHeader(ismrmrd_xml_);
                }
                catch (...)
                {
                    GDEBUG("Failed to write XML header to HDF file\n");
                    return GADGET_FAIL;
                }

                GDEBUG_STREAM("IsmrmrdDumpGadget, save ismrmrd xml header ... ");
            }
            else
            {
                GDEBUG_STREAM("IsmrmrdDumpGadget, do NOT save ismrmrd data ... ");
            }

            first_call_ = false;
        }

        if (this->save_ismrmrd_data_ && !this->save_xml_header_only.value())
        {
            ISMRMRD::Acquisition ismrmrd_acq;
            ismrmrd_acq.setHead(*m1->getObjectPtr());

            GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2 = AsContainerMessage< hoNDArray< std::complex<float> > >(m1->cont());
            if (!m2)
            {
                GDEBUG("Error casting acquisition data package");
                return GADGET_FAIL;
            }

            memcpy((void *)ismrmrd_acq.getDataPtr(), m2->getObjectPtr()->get_data_ptr(), sizeof(float)*m2->getObjectPtr()->get_number_of_elements() * 2);

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
                try
                {
                    ismrmrd_dataset_->appendAcquisition(ismrmrd_acq);
                }
                catch (...)
                {
                    GDEBUG("Error appending ISMRMRD Dataset\n");
                    return GADGET_FAIL;
                }
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

    int IsmrmrdDumpGadget::process(GadgetContainerMessage<ISMRMRD::WaveformHeader>* m1)
    {
        if (first_call_)
        {
            if (this->save_ismrmrd_data_)
            {
                this->create_ismrmrd_dataset();

                try
                {
                    ismrmrd_dataset_->writeHeader(ismrmrd_xml_);
                }
                catch (...)
                {
                    GDEBUG("Failed to write XML header to HDF file\n");
                    return GADGET_FAIL;
                }

                GDEBUG_STREAM("IsmrmrdDumpGadget, first call is on waveform, save ismrmrd xml header ... ");
            }
            else
            {
                GDEBUG_STREAM("IsmrmrdDumpGadget, first call is on waveform, do NOT save ismrmrd data ... ");
            }

            first_call_ = false;
        }

        if (this->save_ismrmrd_data_ && !this->save_xml_header_only.value())
        {
            ISMRMRD::Waveform ismrmrd_wav;
            ismrmrd_wav.head = *m1->getObjectPtr();

            GadgetContainerMessage< hoNDArray< uint32_t > >* m2 = AsContainerMessage< hoNDArray< uint32_t > >(m1->cont());
            if (!m2)
            {
                GDEBUG("Error casting waveform data package");
                return GADGET_FAIL;
            }

            ismrmrd_wav.data = m2->getObjectPtr()->begin();

            try
            {
                ismrmrd_dataset_->appendWaveform(ismrmrd_wav);
            }
            catch (...)
            {
                GDEBUG("Error appending ISMRMRD Waveform\n");
                return GADGET_FAIL;
            }

            ismrmrd_wav.data = NULL;
        }

        // TODO: remove this check
        if(this->pass_waveform_downstream.value())
        {
            if (this->next()->putq(m1) == -1)
            {
                m1->release();
                GERROR("IsmrmrdDumpGadget::process, passing waveform on to next gadget");
                return -1;
            }
        }
        else
        {
            m1->release();
        }

        return 0;
    }

    GADGET_FACTORY_DECLARE(IsmrmrdDumpGadget)
}
