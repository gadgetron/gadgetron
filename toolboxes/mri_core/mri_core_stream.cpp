
#include "mri_core_stream.h"

namespace Gadgetron 
{
    GenericReconIsmrmrdStreamer::GenericReconIsmrmrdStreamer() : verbose_(false)
    {
    }

    GenericReconIsmrmrdStreamer::GenericReconIsmrmrdStreamer(const std::map<std::string, std::string>& parameters) : verbose_(false)
    {
        this->initialize_stream_name_buffer(parameters);
    }

    GenericReconIsmrmrdStreamer::~GenericReconIsmrmrdStreamer()
    {
    }

    void GenericReconIsmrmrdStreamer::initialize_stream_name_buffer(const std::map<std::string, std::string>& parameters)
    {
        this->initialize_stream_name_buffer(parameters, GENERIC_RECON_STREAM_ISMRMRD_HEADER);
        this->initialize_stream_name_buffer(parameters, GENERIC_RECON_STREAM_UNDERSAMPLED_KSPACE);
        this->initialize_stream_name_buffer(parameters, GENERIC_RECON_STREAM_REF_KSPACE);
        this->initialize_stream_name_buffer(parameters, GENERIC_RECON_STREAM_REF_KSPACE_FOR_COILMAP);
        this->initialize_stream_name_buffer(parameters, GENERIC_RECON_STREAM_COILMAP);
        this->initialize_stream_name_buffer(parameters, GENERIC_RECON_STREAM_GFACTOR_MAP);
        this->initialize_stream_name_buffer(parameters, GENERIC_RECON_STREAM_RECONED_KSPACE);
        this->initialize_stream_name_buffer(parameters, GENERIC_RECON_STREAM_RECONED_COMPLEX_IMAGE);
        this->initialize_stream_name_buffer(parameters, GENERIC_RECON_STREAM_RECONED_COMPLEX_IMAGE_AFTER_POSTPROCESSING);
    }

    void GenericReconIsmrmrdStreamer::initialize_stream_name_buffer(const std::map<std::string, std::string>& parameters, const std::string& name)
    {
        if (parameters.find(name) != parameters.end())
        {
            buffer_names_[name].first = parameters.at(name);
            GDEBUG_CONDITION_STREAM(this->verbose_, "Buffer to store " << name << " is " << this->buffer_names_[name].first);
        }
    }

    void GenericReconIsmrmrdStreamer::close_stream_buffer()
    {
        for (auto const& x : this->buffer_names_)
        {
            GDEBUG_CONDITION_STREAM(this->verbose_, "GenericReconIsmrmrdStreamer::close_stream_buffer, stream is for " << x.first << " - " << this->buffer_names_[x.first].first);

            if ( (x.first == GENERIC_RECON_STREAM_ISMRMRD_HEADER)  
                | (x.first == GENERIC_RECON_STREAM_COILMAP) 
                | (x.first == GENERIC_RECON_STREAM_GFACTOR_MAP)
                | (x.first == GENERIC_RECON_STREAM_RECONED_COMPLEX_IMAGE)
                | (x.first == GENERIC_RECON_STREAM_RECONED_COMPLEX_IMAGE_AFTER_POSTPROCESSING)
            )
            {
                if(this->buffer_names_[x.first].second)
                {
                    std::ofstream& os = *this->buffer_names_[x.first].second;
                    if (os.is_open())
                    {
                        GDEBUG_CONDITION_STREAM(this->verbose_, "GenericReconIsmrmrdStreamer::close_stream_buffer, stream is open for " << x.first << "; put in the close message ... ");
                        ISMRMRD::OStreamView ws(os);
                        ISMRMRD::ProtocolSerializer serializer(ws);
                        serializer.close();
                        os.flush();
                    }
                    else
                    {
                        GDEBUG_CONDITION_STREAM(this->verbose_, "GenericReconIsmrmrdStreamer::close_stream_buffer, stream is not open for " << x.first << " ... ");
                    }
                }
            }
        }
    }

    void GenericReconIsmrmrdStreamer::stream_ismrmrd_header(const ISMRMRD::IsmrmrdHeader& hdr)
    {
        if (this->buffer_names_.find(GENERIC_RECON_STREAM_ISMRMRD_HEADER)!=this->buffer_names_.end())
        {
            std::string buf_name = this->buffer_names_[GENERIC_RECON_STREAM_ISMRMRD_HEADER].first;
            if (!this->buffer_names_[GENERIC_RECON_STREAM_ISMRMRD_HEADER].second)
            {
                this->buffer_names_[GENERIC_RECON_STREAM_ISMRMRD_HEADER].second = std::make_shared<std::ofstream>(std::ofstream(buf_name, std::ios::out ));
            }

            std::ofstream& os = *this->buffer_names_[GENERIC_RECON_STREAM_ISMRMRD_HEADER].second;
            if (os.is_open())
            {
                GDEBUG_STREAM("GenericReconIsmrmrdStreamer, stream the ismrmrd header to the array buffer " << buf_name);
                ISMRMRD::serialize(hdr, os);
                os.flush();
            }
            else
            {
                GERROR_STREAM("GenericReconIsmrmrdStreamer, unable to open the ismrmrd header buffer " << buf_name << " ... ");
            }
        }
        else
        {
            GWARN_CONDITION_STREAM(this->verbose_, "GenericReconIsmrmrdStreamer, the pre-set buffer names do not include " << GENERIC_RECON_STREAM_ISMRMRD_HEADER << "; the header will not be saved into the buffer ...");
        }
    }
}