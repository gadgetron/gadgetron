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

void fill_image_array(const hoNDArray<std::complex<float>>& imgs, hoNDArray< mrd::ImageHeader >& headers, hoNDArray< mrd::ImageMeta >& meta)
{
    size_t RO = imgs.get_size(0);
    size_t E1 = imgs.get_size(1);
    size_t E2 = imgs.get_size(2);
    size_t CHA = imgs.get_size(3);
    size_t N = imgs.get_size(4);
    size_t S = imgs.get_size(5);
    size_t SLC = imgs.get_size(6);

    headers.create(N, S, SLC);
    meta.create(N, S, SLC);

    size_t n, s, slc;
    for (slc=0; slc<SLC; slc++)
    {
        for (s=0; s<S; s++)
        {
            for (n=0; n<N; n++)
            {
                headers(n, s, slc).image_index = n+s*N+slc*N*S;
                headers(n, s, slc).field_of_view[0] = RO;
                headers(n, s, slc).field_of_view[1] = E1;
                headers(n, s, slc).field_of_view[2] = E2;

                meta(n, s, slc)["META_Field_1"] = {3.0};
                meta(n, s, slc)["META_Field_2"] = {13.0};
                meta(n, s, slc)["META_Field_3"] = {125.2, 56.4};
            }
        }
    }
}

void remove_parameter_files(const std::map<std::string, std::string>& parameters)
{
    for (auto &m : parameters)
    {
        if (std::filesystem::exists(m.second)) std::remove(m.second.c_str());
    }
}

TEST(GenericReconMrdStreamerTest, test_streamer)
{
    std::string tmp_path = std::string(std::filesystem::temp_directory_path());

    std::map<std::string, std::string> parameters;
    parameters[GENERIC_RECON_STREAM_MRD_HEADER] = tmp_path + "/recon_header.dat";
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
        GenericReconMrdStreamer gt_streamer(parameters);
        gt_streamer.verbose_ = true;

        mrd::Header hdr;
        hdr.encoding.resize(1);
        hdr.encoding[0].recon_space.matrix_size.x = 323;
        hdr.encoding[0].recon_space.matrix_size.y = 458;
        hdr.encoding[0].recon_space.matrix_size.z = 7;
        gt_streamer.stream_mrd_header(hdr);

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

        hoNDArray< mrd::ImageHeader > headers;
        hoNDArray< mrd::ImageMeta > meta;

        fill_image_array(imgs, headers, meta);
        GDEBUG_STREAM("fill_image_array ...");

        gt_streamer.stream_to_mrd_image_buffer(GENERIC_RECON_STREAM_RECONED_COMPLEX_IMAGE, imgs, headers, meta);
        gt_streamer.stream_to_mrd_image_buffer(GENERIC_RECON_STREAM_RECONED_COMPLEX_IMAGE_AFTER_POSTPROCESSING, imgs, headers, meta);

        // waveform
        std::vector<mrd::WaveformUint32> wavs(12);
        uint16_t num_samples=24, cha=4;
        mrd::WaveformUint32 a_wav;
        a_wav.data.create(cha, num_samples);
        uint32_t* p_data = a_wav.data.begin();
        uint32_t* p_data_end = a_wav.data.end();

        std::default_random_engine generator(35554);
        std::normal_distribution<float> distribution(8.0, 47.0);
        std::generate(p_data, p_data_end, [&distribution, &generator]() { return distribution(generator); });

        for (auto i=0; i<wavs.size(); i++) wavs[i] = a_wav;

        gt_streamer.stream_mrd_waveforms(wavs);

        // close all streams
        gt_streamer.close_stream_buffer();

        std::optional<mrd::Header> header_deserialized;
        mrd::StreamItem item;
        hoNDArray<std::complex<float>> diff;

        // deserialize and check MRD Header
        mrd::binary::MrdReader header_reader(parameters[GENERIC_RECON_STREAM_MRD_HEADER]);
        header_reader.ReadHeader(header_deserialized);
        while (header_reader.ReadData(item)) {}
        ASSERT_EQ(hdr.encoding[0].recon_space.matrix_size.x, header_deserialized->encoding[0].recon_space.matrix_size.x);
        ASSERT_EQ(hdr.encoding[0].recon_space.matrix_size.y, header_deserialized->encoding[0].recon_space.matrix_size.y);
        ASSERT_EQ(hdr.encoding[0].recon_space.matrix_size.z, header_deserialized->encoding[0].recon_space.matrix_size.z);

        // deserialize and check undersampled kspace
        mrd::binary::MrdReader fd_reader(parameters[GENERIC_RECON_STREAM_UNDERSAMPLED_KSPACE]);
        fd_reader.ReadHeader(header_deserialized);
        while (fd_reader.ReadData(item)) {
            auto dat = std::get<mrd::ArrayComplexFloat>(item);
            Gadgetron::subtract(data, dat, diff);
            float v = Gadgetron::nrm2(diff);
            EXPECT_LE(v, 0.001);
        }

        // deserialize and check ref kspace
        mrd::binary::MrdReader fd_ref_reader(parameters[GENERIC_RECON_STREAM_REF_KSPACE]);
        fd_ref_reader.ReadHeader(header_deserialized);
        while (fd_ref_reader.ReadData(item)) {
            auto dat = std::get<mrd::ArrayComplexFloat>(item);
            Gadgetron::subtract(ref, dat, diff);
            float v = Gadgetron::nrm2(diff);
            EXPECT_LE(v, 0.001);
        }

        // deserialize and check images
        mrd::binary::MrdReader img_reader(parameters[GENERIC_RECON_STREAM_RECONED_COMPLEX_IMAGE_AFTER_POSTPROCESSING]);
        img_reader.ReadHeader(header_deserialized);
        int ind = 0;
        while (img_reader.ReadData(item)) {
            auto img = std::get<mrd::ImageComplexFloat>(item);

            int slc = ind / (N*S);
            int s = (ind - slc*N*S) / N;
            int n = ind - s*N - slc*N*S;
            GDEBUG_STREAM("ProtocolDeserializer for image " << ind << " - " << n << " " << s << " " << slc);

            hoNDArray<std::complex<float>> a_img(RO, E1, E2, CHA, &imgs(0, 0, 0, 0, n, s, slc));

            Gadgetron::subtract(a_img, img.data, diff);
            float v = Gadgetron::nrm2(diff);
            EXPECT_LE(v, 0.001);

            ASSERT_EQ(headers(n, s, slc).image_index, img.head.image_index);
            ASSERT_EQ(headers(n, s, slc).field_of_view, img.head.field_of_view);

            ASSERT_FLOAT_EQ(std::get<double>(meta(n, s, slc)["META_Field_1"][0]), std::get<double>(img.meta["META_Field_1"][0]));
            ASSERT_FLOAT_EQ(std::get<double>(meta(n, s, slc)["META_Field_2"][0]), std::get<double>(img.meta["META_Field_2"][0]));
            ASSERT_FLOAT_EQ(std::get<double>(meta(n, s, slc)["META_Field_3"][0]), std::get<double>(img.meta["META_Field_3"][0]));
            ASSERT_FLOAT_EQ(std::get<double>(meta(n, s, slc)["META_Field_3"][1]), std::get<double>(img.meta["META_Field_3"][1]));

            ind++;
        }

        {
            std::vector<mrd::WaveformUint32> wavs_deserialized;
            mrd::binary::MrdReader wav_reader(parameters[GENERIC_RECON_STREAM_WAVEFORM]);
            wav_reader.ReadHeader(header_deserialized);
            while (wav_reader.ReadData(item)) {
                auto wav = std::get<mrd::WaveformUint32>(item);
                wavs_deserialized.push_back(wav);
            }

            ASSERT_EQ(wavs_deserialized.size(), wavs.size());

            for (auto i=0; i<wavs.size(); i++)
            {
                ASSERT_EQ(wavs_deserialized[i].NumberOfSamples(), wavs[i].NumberOfSamples());
                ASSERT_EQ(wavs_deserialized[i].Channels(), wavs[i].Channels());

                size_t N = wavs[i].data.size();
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
