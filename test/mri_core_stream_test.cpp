#include <gtest/gtest.h>
#include "hoNDImage.h"
#include "hoMRImage.h"
#include "hoNDArray_reductions.h"
#include "hoNDArray_elemwise.h"
#include "mri_core_stream.h"
#include <numeric>
#include <random>
#include <filesystem>
#include <iostream>

void fill_random(hoNDArray<std::complex<float>>& data, int seed = 35879)
{
    std::default_random_engine generator(seed);
    std::normal_distribution<float> distribution(2.0, 12.0);

    std::generate(data.begin(), data.end(), [&distribution, &generator]() { return distribution(generator); });
}

void fill_image_array(const hoNDArray<std::complex<float>>& imgs, hoNDArray< ISMRMRD::ImageHeader >& headers, std::vector< ISMRMRD::MetaContainer >& meta)
{
    size_t RO = imgs.get_size(0);
    size_t E1 = imgs.get_size(1);
    size_t E2 = imgs.get_size(2);
    size_t CHA = imgs.get_size(3);
    size_t N = imgs.get_size(4);
    size_t S = imgs.get_size(5);
    size_t SLC = imgs.get_size(6);

    headers.create(N, S, SLC);

    size_t n, s, slc;
    for (slc=0; slc<SLC; slc++)
    {
        for (s=0; s<S; s++)
        {
            for (n=0; n<N; n++)
            {
                headers(n, s, slc).image_index = n+s*N+slc*N*S;
                headers(n, s, slc).data_type = Gadgetron::Core::IO::ismrmrd_data_type<std::complex<float>>();
                headers(n, s, slc).matrix_size[0] = RO;
                headers(n, s, slc).matrix_size[1] = E1;
                headers(n, s, slc).matrix_size[2] = E2;
                headers(n, s, slc).channels = CHA;
                headers(n, s, slc).field_of_view[0] = RO;
                headers(n, s, slc).field_of_view[1] = E1;
                headers(n, s, slc).field_of_view[2] = E2;
            }
        }
    }

    meta.resize(N*S*SLC);
    for (auto &m : meta)
    {
        m.set("META_Field_1", 3.0);
        m.set("META_Field_2", 13.0);
        m.set("META_Field_3", 125.2);
        m.append("META_Field_3", 56.4);
    }
}

void remove_parameter_files(const std::map<std::string, std::string>& parameters) 
{
    for (auto &m : parameters)
    {
        if (std::filesystem::exists(m.second)) std::remove(m.second.c_str());
    }
}

TEST(GenericReconIsmrmrdStreamerTest, test_streamer)
{
    std::string tmp_path = std::string(std::filesystem::temp_directory_path());

    std::map<std::string, std::string> parameters;
    parameters[GENERIC_RECON_STREAM_ISMRMRD_HEADER] = tmp_path + "/recon_header.dat";
    parameters[GENERIC_RECON_STREAM_UNDERSAMPLED_KSPACE] = tmp_path + "/undersampled_kspace.dat";
    parameters[GENERIC_RECON_STREAM_REF_KSPACE] = tmp_path + "/ref_kspace.dat";
    parameters[GENERIC_RECON_STREAM_REF_KSPACE_FOR_COILMAP] = tmp_path + "/ref_kspace_for_coil_map.dat";
    parameters[GENERIC_RECON_STREAM_COILMAP] = tmp_path + "/coil_map.dat";
    parameters[GENERIC_RECON_STREAM_GFACTOR_MAP] = tmp_path + "/gfactor.dat";
    parameters[GENERIC_RECON_STREAM_RECONED_KSPACE] = tmp_path + "/reconed_kspace.dat";
    parameters[GENERIC_RECON_STREAM_RECONED_COMPLEX_IMAGE] = tmp_path + "/reconed_images.dat";
    parameters[GENERIC_RECON_STREAM_RECONED_COMPLEX_IMAGE_AFTER_POSTPROCESSING] = tmp_path + "/reconed_images_after_post_processing.dat";
    parameters[GENERIC_RECON_STREAM_WAVEFORM] = tmp_path + "/waveform.dat";

    try
    {
        GenericReconIsmrmrdStreamer gt_streamer(parameters);
        gt_streamer.verbose_ = true;

        ISMRMRD::IsmrmrdHeader hdr;
        hdr.encoding.resize(1);
        hdr.encoding[0].reconSpace.matrixSize.x = 323;
        hdr.encoding[0].reconSpace.matrixSize.y = 458;
        hdr.encoding[0].reconSpace.matrixSize.z = 7;
        gt_streamer.stream_ismrmrd_header(hdr);

        // kspace data
        hoNDArray<std::complex<float>> data;
        data.create(32, 45, 67, 8);
        fill_random(data);

        gt_streamer.stream_to_array_buffer(GENERIC_RECON_STREAM_UNDERSAMPLED_KSPACE, data);
        gt_streamer.stream_to_array_buffer(GENERIC_RECON_STREAM_RECONED_KSPACE, data);

        hoNDArray<std::complex<float>> ref;
        ref.create(64, 32, 32, 8, 2, 1);
        fill_random(ref);
        gt_streamer.stream_to_array_buffer(GENERIC_RECON_STREAM_REF_KSPACE, ref);
        gt_streamer.stream_to_array_buffer(GENERIC_RECON_STREAM_REF_KSPACE_FOR_COILMAP, ref);

        // images
        hoNDArray<std::complex<float>> imgs;

        size_t RO=32, E1=16, E2=1, CHA=4;
        size_t N=3, S=2, SLC=2;

        imgs.create(RO, E1, E2, CHA, N, S, SLC);
        fill_random(imgs);
        GDEBUG_STREAM("fill in imgs ...");

        hoNDArray< ISMRMRD::ImageHeader > headers;
        std::vector< ISMRMRD::MetaContainer > meta;

        fill_image_array(imgs, headers, meta);
        GDEBUG_STREAM("fill_image_array ...");

        gt_streamer.stream_to_ismrmrd_image_buffer(GENERIC_RECON_STREAM_RECONED_COMPLEX_IMAGE, imgs, headers, meta);
        gt_streamer.stream_to_ismrmrd_image_buffer(GENERIC_RECON_STREAM_RECONED_COMPLEX_IMAGE_AFTER_POSTPROCESSING, imgs, headers, meta);

        // waveform
        std::vector<ISMRMRD::Waveform> wavs(12);

        uint16_t num_samples=24, cha=4;
        ISMRMRD::Waveform a_wav(num_samples, cha);
        uint32_t* p_data = a_wav.begin_data();
        uint32_t* p_data_end = a_wav.end_data();

        std::default_random_engine generator(35554);
        std::normal_distribution<float> distribution(8.0, 47.0);
        std::generate(p_data, p_data_end, [&distribution, &generator]() { return distribution(generator); });

        for (auto i=0; i<wavs.size(); i++) wavs[i] = a_wav;

        gt_streamer.stream_ismrmrd_waveform(wavs);

        // close all streams
        gt_streamer.close_stream_buffer();

        // deserialize
        float v;
        {
            std::ifstream fd(parameters[GENERIC_RECON_STREAM_ISMRMRD_HEADER], std::ios::in | std::ios::binary);

            ISMRMRD::IStreamView rs(fd);
            ISMRMRD::ProtocolDeserializer deserializer(rs);

            ISMRMRD::IsmrmrdHeader hdr_deserialized;
            deserializer.deserialize(hdr_deserialized);

            ASSERT_EQ(hdr.encoding[0].reconSpace.matrixSize.x, hdr_deserialized.encoding[0].reconSpace.matrixSize.x);
            ASSERT_EQ(hdr.encoding[0].reconSpace.matrixSize.y, hdr_deserialized.encoding[0].reconSpace.matrixSize.y);
            ASSERT_EQ(hdr.encoding[0].reconSpace.matrixSize.z, hdr_deserialized.encoding[0].reconSpace.matrixSize.z);
        }
        {
            std::ifstream fd(parameters[GENERIC_RECON_STREAM_UNDERSAMPLED_KSPACE], std::ios::in | std::ios::binary);
            ASSERT_EQ(fd.is_open(), true);

            ISMRMRD::IStreamView rs(fd);
            ISMRMRD::ProtocolDeserializer deserializer(rs);

            ISMRMRD::NDArray<std::complex<float>> arr;
            deserializer.deserialize(arr);

            hoNDArray<std::complex<float>> data_deserialized, diff;
            Gadgetron::convert_ismrmrd_ndarray_to_hoNDArray(arr, data_deserialized);

            Gadgetron::subtract(data, data_deserialized, diff);
            v = Gadgetron::nrm2(diff);
            EXPECT_LE(v, 0.001);
        }
        {
            std::ifstream fd_ref(parameters[GENERIC_RECON_STREAM_REF_KSPACE], std::ios::in | std::ios::binary);

            ISMRMRD::IStreamView rs(fd_ref);
            ISMRMRD::ProtocolDeserializer deserializer(rs);

            ISMRMRD::NDArray<std::complex<float>> arr;
            deserializer.deserialize(arr);

            hoNDArray<std::complex<float>> ref_deserialized, diff;
            Gadgetron::convert_ismrmrd_ndarray_to_hoNDArray(arr, ref_deserialized);

            Gadgetron::subtract(ref, ref_deserialized, diff);
            v = Gadgetron::nrm2(diff);
            EXPECT_LE(v, 0.001);
        }
        {
            std::ifstream is(parameters[GENERIC_RECON_STREAM_RECONED_COMPLEX_IMAGE_AFTER_POSTPROCESSING].c_str(), std::ios::binary);
            ASSERT_EQ(is.is_open(), true);

            ISMRMRD::IStreamView rs(is);
            ISMRMRD::ProtocolDeserializer deserializer(rs);

            int ind = 0;
            while (deserializer.peek() != ISMRMRD::ISMRMRD_MESSAGE_CLOSE)
            {
                int n, s, slc;
                slc = ind / (N*S);
                s = (ind - slc*N*S) / N;
                n = ind - s*N - slc*N*S;
                GDEBUG_STREAM("ProtocolDeserializer for image " << ind << " - " << n << " " << s << " " << slc);

                ASSERT_EQ(deserializer.peek(), ISMRMRD::ISMRMRD_MESSAGE_IMAGE);
                ASSERT_EQ(deserializer.peek_image_data_type(), ISMRMRD::ISMRMRD_CXFLOAT);

                ISMRMRD::Image<std::complex<float> > img;
                deserializer.deserialize(img);

                hoNDArray<std::complex<float>> a_img, a_img_deserialized, diff;
                a_img.create(RO, E1, E2, CHA, &imgs(0, 0, 0, 0, n, s, slc));
                a_img_deserialized.create(RO, E1, E2, CHA, img.getDataPtr());

                Gadgetron::subtract(a_img, a_img_deserialized, diff);
                v = Gadgetron::nrm2(diff);
                EXPECT_LE(v, 0.001);

                ASSERT_EQ(headers(n, s, slc).image_index, img.getHead().image_index);
                ASSERT_EQ(headers(n, s, slc).data_type, img.getHead().data_type);
                ASSERT_EQ(headers(n, s, slc).matrix_size[0], img.getHead().matrix_size[0]);
                ASSERT_EQ(headers(n, s, slc).matrix_size[1], img.getHead().matrix_size[1]);
                ASSERT_EQ(headers(n, s, slc).matrix_size[2], img.getHead().matrix_size[2]);
                ASSERT_EQ(headers(n, s, slc).channels, img.getHead().channels);
                ASSERT_EQ(headers(n, s, slc).field_of_view[0], img.getHead().field_of_view[0]);
                ASSERT_EQ(headers(n, s, slc).field_of_view[1], img.getHead().field_of_view[1]);
                ASSERT_EQ(headers(n, s, slc).field_of_view[2], img.getHead().field_of_view[2]);

                ISMRMRD::MetaContainer a_meta;
                ISMRMRD::deserialize(img.getAttributeString(), a_meta);
                ASSERT_FLOAT_EQ(meta[ind].as_double("META_Field_1"), a_meta.as_double("META_Field_1"));
                ASSERT_FLOAT_EQ(meta[ind].as_double("META_Field_2"), a_meta.as_double("META_Field_2"));
                ASSERT_FLOAT_EQ(meta[ind].as_double("META_Field_3"), a_meta.as_double("META_Field_3"));
                ASSERT_FLOAT_EQ(meta[ind].as_double("META_Field_3", 1), a_meta.as_double("META_Field_3", 1));

                ind++;
            }
        }
        {
            std::vector<ISMRMRD::Waveform> wavs_deserialized;

            std::ifstream is(parameters[GENERIC_RECON_STREAM_WAVEFORM].c_str(), std::ios::binary);
            ASSERT_EQ(is.is_open(), true);

            ISMRMRD::IStreamView rs(is);
            ISMRMRD::ProtocolDeserializer deserializer(rs);

            while (deserializer.peek() != ISMRMRD::ISMRMRD_MESSAGE_CLOSE)
            {
                ISMRMRD::Waveform a_wav;
                deserializer.deserialize(a_wav);
                wavs_deserialized.push_back(a_wav);
            }

            ASSERT_EQ(wavs_deserialized.size(), wavs.size());

            for (auto i=0; i<wavs.size(); i++)
            {
                ASSERT_EQ(wavs_deserialized[i].head.number_of_samples, wavs[i].head.number_of_samples);
                ASSERT_EQ(wavs_deserialized[i].head.channels, wavs[i].head.channels);

                size_t N = wavs[i].size();
                for (auto k=0; k<N; k++)
                {
                    ASSERT_FLOAT_EQ(wavs_deserialized[i].data[k], wavs[i].data[k]);
                }
            }
        }
        // clean the files
        remove_parameter_files(parameters);
    }
    catch(...)
    {
        remove_parameter_files(parameters);
        throw;
    }
}
