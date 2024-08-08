
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
        this->initialize_stream_name_buffer(parameters, GENERIC_RECON_STREAM_WAVEFORM);
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

    std::shared_ptr<std::ofstream> GenericReconIsmrmrdStreamer::find_and_open_stream(const std::string& name, std::string& buf_name)
    {
        if (this->buffer_names_.find(name)!=this->buffer_names_.end())
        {
            buf_name = this->buffer_names_[name].first;

            if (!this->buffer_names_[name].second)
            {
                GDEBUG_STREAM("Create the stream for the first time - " << buf_name);
                this->buffer_names_[name].second = std::make_shared<std::ofstream>(std::ofstream(buf_name, std::ios::out | std::ios::binary | std::ios::app));
            }

            return this->buffer_names_[name].second;
        }
        else
        {
            GWARN_CONDITION_STREAM(this->verbose_, "The pre-set buffer names do not include " << name << " ...");
            return std::shared_ptr<std::ofstream>();
        }
    }

    void GenericReconIsmrmrdStreamer::stream_ismrmrd_header(const ISMRMRD::IsmrmrdHeader& hdr)
    {
        std::string buf_name;
        std::shared_ptr<std::ofstream> os = find_and_open_stream(GENERIC_RECON_STREAM_ISMRMRD_HEADER, buf_name);
        if (os && os->is_open())
        {
            GDEBUG_STREAM("GenericReconIsmrmrdStreamer, stream the ismrmrd header to the array buffer " << buf_name);
            ISMRMRD::OStreamView ws(*os);
            ISMRMRD::ProtocolSerializer serializer(ws);
            serializer.serialize(hdr);
            os->flush();
        }
    }

    void GenericReconIsmrmrdStreamer::stream_ismrmrd_waveform(const std::vector< ISMRMRD::Waveform>& wav)
    {
        std::string buf_name;
        std::shared_ptr<std::ofstream> os = find_and_open_stream(GENERIC_RECON_STREAM_WAVEFORM, buf_name);
        if (os && os->is_open())
        {
            GDEBUG_STREAM("GenericReconIsmrmrdStreamer, stream the waveform to buffer " << buf_name);
            ISMRMRD::OStreamView ws(*os);
            ISMRMRD::ProtocolSerializer serializer(ws);

            for (auto w : wav)
            {
                serializer.serialize(w);
            }
            os->flush();
        }
    }

    template <typename T> 
    void convert_hoNDArray_to_ismrmrd_ndarray(const hoNDArray<T>& ho_arr, ISMRMRD::NDArray<T>& arr)
    {
        size_t NDim = ho_arr.get_number_of_dimensions();
        if (NDim > ISMRMRD::ISMRMRD_NDARRAY_MAXDIM)
        {
            GWARN_STREAM("convert_hoNDArray_to_ismrmrd_ndarray, ho_arr, number of dimensions > ISMRMRD_NDARRAY_MAXDIM" << NDim << " > " << ISMRMRD::ISMRMRD_NDARRAY_MAXDIM);
            return;
        }

        std::vector<size_t> dims;
        ho_arr.get_dimensions(dims);

        arr.resize(dims);
        memcpy(arr.getDataPtr(), ho_arr.get_data_ptr(), ho_arr.get_number_of_bytes());
    }

    template void convert_hoNDArray_to_ismrmrd_ndarray(const hoNDArray<short>& ho_arr, ISMRMRD::NDArray<short>& arr);
    template void convert_hoNDArray_to_ismrmrd_ndarray(const hoNDArray<unsigned short>& ho_arr, ISMRMRD::NDArray<unsigned short>& arr);
    template void convert_hoNDArray_to_ismrmrd_ndarray(const hoNDArray<int>& ho_arr, ISMRMRD::NDArray<int>& arr);
    template void convert_hoNDArray_to_ismrmrd_ndarray(const hoNDArray<unsigned int>& ho_arr, ISMRMRD::NDArray<unsigned int>& arr);
    template void convert_hoNDArray_to_ismrmrd_ndarray(const hoNDArray<float>& ho_arr, ISMRMRD::NDArray<float>& arr);
    template void convert_hoNDArray_to_ismrmrd_ndarray(const hoNDArray<double>& ho_arr, ISMRMRD::NDArray<double>& arr);
    template void convert_hoNDArray_to_ismrmrd_ndarray(const hoNDArray< std::complex<float> >& ho_arr, ISMRMRD::NDArray< std::complex<float> >& arr);
    template void convert_hoNDArray_to_ismrmrd_ndarray(const hoNDArray< std::complex<double> >& ho_arr, ISMRMRD::NDArray< std::complex<double> >& arr);

    template <typename T> 
    void convert_ismrmrd_ndarray_to_hoNDArray(const ISMRMRD::NDArray<T>& arr, hoNDArray<T>& ho_arr)
    {
        size_t NDim = arr.getNDim();

        std::vector<size_t> dim(NDim);
        for (auto i=0; i<NDim; i++) dim[i] = arr.getDims()[i];

        ho_arr.create(dim);
        memcpy(ho_arr.get_data_ptr(), arr.getDataPtr(), ho_arr.get_number_of_bytes());
    }

    template void convert_ismrmrd_ndarray_to_hoNDArray(const ISMRMRD::NDArray<short>& arr, hoNDArray<short>& ho_arr);
    template void convert_ismrmrd_ndarray_to_hoNDArray(const ISMRMRD::NDArray<unsigned short>& arr, hoNDArray<unsigned short>& ho_arr);
    template void convert_ismrmrd_ndarray_to_hoNDArray(const ISMRMRD::NDArray<int>& arr, hoNDArray<int>& ho_arr);
    template void convert_ismrmrd_ndarray_to_hoNDArray(const ISMRMRD::NDArray<unsigned int>& arr, hoNDArray<unsigned int>& ho_arr);
    template void convert_ismrmrd_ndarray_to_hoNDArray(const ISMRMRD::NDArray<float>& arr, hoNDArray<float>& ho_arr);
    template void convert_ismrmrd_ndarray_to_hoNDArray(const ISMRMRD::NDArray<double>& arr, hoNDArray<double>& ho_arr);
    template void convert_ismrmrd_ndarray_to_hoNDArray(const ISMRMRD::NDArray< std::complex<float> >& arr, hoNDArray< std::complex<float> >& ho_arr);
    template void convert_ismrmrd_ndarray_to_hoNDArray(const ISMRMRD::NDArray< std::complex<double> >& arr, hoNDArray< std::complex<double> >& ho_arr);
}