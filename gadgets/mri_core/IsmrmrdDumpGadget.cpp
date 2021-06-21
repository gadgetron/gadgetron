#include "IsmrmrdDumpGadget.h"
#include <iomanip>
#include <boost/filesystem.hpp>
#include "network_utils.h"
#include <thread>

namespace bf = boost::filesystem;


namespace Gadgetron
{
    static std::string get_date_time_string()
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
    bool IsmrmrdDumpGadget::is_ip_on_blacklist() const {
        if (ip_no_data_saving.empty()){
            return false; 
        }
        GDEBUG_STREAM("IsmrmrdDumpGadget, find pre-set ip for no-data-saving : " << ip_no_data_saving.size() << " [ ");
        for (auto ip : ip_no_data_saving)
                GDEBUG_STREAM(ip);
        GDEBUG_STREAM(" ] ");

            const auto [host_name, gt_ip_list ] = Gadgetron::find_gadgetron_ip();

            GDEBUG_STREAM("IsmrmrdDumpGadget, find gadgetron host name : " << host_name);

            GDEBUG_STREAM("IsmrmrdDumpGadget, find gadgetron ip : " << gt_ip_list.size() << " [ ");
            for (const auto & ip : gt_ip_list) 
                GDEBUG_STREAM(ip);
            GDEBUG_STREAM(" ] ");


            for (const auto & ip : gt_ip_list ){
                if (ip_no_data_saving.count(ip)){

                    GDEBUG_STREAM("IsmrmrdDumpGadget, find matching ip pair : " << ip);
                    return true;
                }
            }
            return false;
       }


    

    IsmrmrdDumpGadget::IsmrmrdDumpGadget(const Core::Context& context, const Core::GadgetProperties& props ) : Core::ChannelGadget<Core::variant<Core::Acquisition,Core::Waveform>>(context, props), save_ismrmrd_data_(!is_ip_on_blacklist()) {
        if (save_ismrmrd_data_)
        {
            if (!exists(folder ))
            {
                try
                {
                    boost::filesystem::create_directory(folder );
                }
                catch (...)
                {
                    std::stringstream stream("Error caught trying to create folder ");
                    stream << folder;
                    GADGET_THROW(stream.str());
                }
            }
            else
            {
                if (!is_directory(folder ))
                {
                    GADGET_THROW("Specified path is not a directory\n");
                }
            }
        }
    }

    ISMRMRD::Dataset IsmrmrdDumpGadget::create_ismrmrd_dataset() const 
    {
            std::string measurement_id = "";
            std::string ismrmrd_filename = "";

            if (header.measurementInformation)
            {
                if (header.measurementInformation->measurementID)
                {
                    measurement_id = *header.measurementInformation->measurementID;
                }
            }

            GDEBUG("Measurement ID: %s\n", measurement_id.c_str());


            if (file_prefix.empty())
            {
                // try to use the protocol name
                if (header.measurementInformation.is_present())
                {
                    if (header.measurementInformation.get().protocolName.is_present())
                    {
                        ismrmrd_filename = header.measurementInformation.get().protocolName.get();

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
                ismrmrd_filename = file_prefix;
            }

            ismrmrd_filename.append("_");
            ismrmrd_filename.append(measurement_id);

            std::string study_date, study_time;
            if (header.studyInformation)
            {
                if (header.studyInformation->studyDate)
                {
                    study_date = *header.studyInformation->studyDate;

                    std::string d(study_date);
                    d.erase(std::remove(d.begin(), d.end(), '-'), d.end());
                    study_date = d;
                }

                if (header.studyInformation->studyTime)
                {
                    study_time = *header.studyInformation->studyTime;

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

            auto filepath = folder / ismrmrd_filename;

            ismrmrd_filename = filepath.string();
            GDEBUG_STREAM("KSpace dump file name : " << ismrmrd_filename);

            return ISMRMRD::Dataset(ismrmrd_filename.c_str(), "dataset", true);
        }

    static void append_to_dataset(const Core::Acquisition& acq, ISMRMRD::Dataset& dataset){
        const auto& [acq_head, data, traj]  = acq;         
                   ISMRMRD::Acquisition ismrmrd_acq;
                   ismrmrd_acq.setHead(acq_head);
                   ismrmrd_acq.setData(const_cast<std::complex<float>*>(data.data()));
                   if (traj) ismrmrd_acq.setTraj(const_cast<float*>(traj->data()));
                   dataset.appendAcquisition(ismrmrd_acq);
    }
    static void append_to_dataset(const Core::Waveform& acq, ISMRMRD::Dataset& dataset){
        const auto& [wav_head, data ]  = acq;         
                   ISMRMRD::Waveform ismrmrd_wav(wav_head.number_of_samples,wav_head.channels);
                   ismrmrd_wav.head = wav_head;
                   std::copy(data.begin(),data.end(),ismrmrd_wav.data);
                   dataset.appendWaveform(ismrmrd_wav);
    }


    void IsmrmrdDumpGadget::process(Core::InputChannel<Core::variant<Core::Acquisition,Core::Waveform>>& input, Core::OutputChannel& output)
    {

        auto is_valid_type =[this](const auto& item){ return std::holds_alternative<Core::Acquisition>(item);};

        auto move_if = [](auto& input, auto& output, const auto& pred ){
            for (auto item : input ){
                if (pred(item)){
                    output.push(std::move(item));
                }
            }
        };

        if (!save_ismrmrd_data_) {
            GDEBUG_STREAM("IsmrmrdDumpGadget, do NOT save ismrmrd data ... ");
            move_if(input,output, is_valid_type );
            return;
        }



        Core::MPMCChannel<Core::variant<Core::Acquisition,Core::Waveform>> data_buffer;

        auto save_thread = std::thread([&data_buffer,this](){
            auto dataset = create_ismrmrd_dataset();

            {
                auto stream = std::stringstream();
                ISMRMRD::serialize(header,stream);
                dataset.writeHeader(stream.str());
                GDEBUG_STREAM("IsmrmrdDumpGadget, save ismrmrd xml header ... ");
            }


            try {
                for (;;) {
                    Core::visit([&dataset](const auto& item) { append_to_dataset(item, dataset); }, data_buffer.pop());
                }
            } catch (const Core::ChannelClosed& closed) {              
            }
        });

        if (save_xml_header_only){
            GDEBUG_STREAM("Only saving header");
            data_buffer.close();
            move_if(input,output, is_valid_type);
            save_thread.join();
            return;
        }
        
        for (auto item : input){
            data_buffer.push(item);
            if (is_valid_type(item))
                output.push(std::move(item));
        }
        data_buffer.close();
        save_thread.join();
    }
    GADGETRON_GADGET_EXPORT(IsmrmrdDumpGadget);



}
