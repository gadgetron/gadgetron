#include "python_toolbox.h"
#include "GadgetronTimer.h"
#include <gtest/gtest.h>
#include "hoNDArray_utils.h"
#include "hoNDArray_elemwise.h"

using namespace Gadgetron;
using testing::Types;

class python_converter_test : public ::testing::Test
{
protected:
    virtual void SetUp()
    {
        Gadgetron::initialize_python();
    }
};


TEST_F(python_converter_test, no_return_value)
{
    int a = -42;
    float b = 3.141592;
    std::string c("hello, world");
    unsigned int d(117);
    std::complex<double> e(2.12894, -1.103103);

    std::vector<size_t> vec = {1, 2, 3, 4};

    PythonFunction<> foo("builtins", "print");
    foo(a, b, c, d, e, vec);
}

TEST_F(python_converter_test, single_return_value)
{
    PythonFunction<float> atan2("math", "atan2");
    int x = 7, y = 4;
    float atan = atan2(x, y);
    EXPECT_FLOAT_EQ(atan, 1.05165);
}

TEST_F(python_converter_test, tuple_return_value)
{
    PythonFunction<float, float> divmod("builtins", "divmod");
    float w = 6.89;
    float z = 4.12;
    float fsum = 0, fdiff = 0;
    std::tie(fsum, fdiff) = divmod(w, z);
    EXPECT_FLOAT_EQ(fsum, 1);
    EXPECT_FLOAT_EQ(fdiff, 2.77);
}

TEST_F(python_converter_test, tuple_len)
{
    PythonFunction<int> tuplen("builtins", "len");
    int l = tuplen(std::make_tuple(-7, 0, 7));
    EXPECT_EQ(l, 3);
}

TEST_F(python_converter_test, std_vec_complex)
{
    {
        GILLock gl;     // this is needed
        boost::python::object main(boost::python::import("__main__"));
        boost::python::object global(main.attr("__dict__"));
        boost::python::exec("from numpy.random import random\n"
            "def rand_cplx_array(length): \n"
            "    return random(length) + 1j * random(length)\n",
            global, global);
    }

    std::vector<std::complex<double>> vec;
    PythonFunction<std::vector<std::complex<double>>> make_vec("__main__", "rand_cplx_array");
    vec = make_vec(32);
    EXPECT_EQ(vec.size(), 32);
}

TEST_F(python_converter_test, numpy_hoNDArray)
{
    PythonFunction<hoNDArray<float>> arange("numpy", "arange");
    hoNDArray<float> evens = arange(0, 100, 2, "f");
    EXPECT_FLOAT_EQ(evens(0), 0);
    EXPECT_FLOAT_EQ(evens(1), 2);
    EXPECT_FLOAT_EQ(evens(2), 4);
    EXPECT_EQ(evens.get_number_of_elements(), 50);
}

TEST_F(python_converter_test, numpy_hoNDArray_two_inputs)
{
    {
        GILLock gl;     // this is needed
        boost::python::object main(boost::python::import("__main__"));
        boost::python::object global(main.attr("__dict__"));
        boost::python::exec("import numpy as np\n"
            "def scale_array(a): \n"
            "   b = a + 100\n"
            "   return a,b\n",
            global, global);
    }

    hoNDArray<float> a, a2, b;
    a.create(32, 64);
    Gadgetron::fill(a, float(45));

    PythonFunction< hoNDArray<float>, hoNDArray<float> > scale_array_test("__main__", "scale_array");
    std::tie(a2, b) = scale_array_test(a);

    EXPECT_FLOAT_EQ(a2[0], 45);
    EXPECT_FLOAT_EQ(b[12], 145);
}

TEST_F(python_converter_test, numpy_hoNDArray_three_outputs)
{
    {
        GILLock gl;     // this is needed
        boost::python::object main(boost::python::import("__main__"));
        boost::python::object global(main.attr("__dict__"));
        boost::python::exec("import numpy as np\n"
            "def scale_array(a, b): \n"
            "   c = a + b\n"
            "   return a,b,c\n",
            global, global);
    }

    hoNDArray<float> a, b;
    a.create(32, 64);
    Gadgetron::fill(a, float(45));

    b.create(32, 64);
    Gadgetron::fill(b, float(210));

    hoNDArray<float> a2, b2, c;

    PythonFunction< hoNDArray<float>, hoNDArray<float>, hoNDArray<float> > scale_array_test("__main__", "scale_array");
    std::tie(a2, b2, c) = scale_array_test(a, b);

    EXPECT_FLOAT_EQ(a2[0], 45);
    EXPECT_FLOAT_EQ(b2[12], 210);
    EXPECT_FLOAT_EQ(c[20], 255);
}

TEST_F(python_converter_test, mrd_acquisitionheader)
{
    {
        GILLock gl;     // this is needed
        boost::python::object main(boost::python::import("__main__"));
        boost::python::object global(main.attr("__dict__"));
        boost::python::exec("def modify(head): head.measurement_uid += 1; return head",
            global, global);
    }

    mrd::AcquisitionHeader acq_head;
    acq_head.flags.SetFlags(mrd::AcquisitionFlags::kIsNoiseMeasurement);
    acq_head.measurement_uid = 41;
    acq_head.idx.repetition = 1;
    acq_head.center_sample = 64;

    PythonFunction<mrd::AcquisitionHeader> modify_acq_header("__main__", "modify");
    mrd::AcquisitionHeader acq_head2 = modify_acq_header(acq_head);
    EXPECT_TRUE(acq_head2.flags.HasFlags(mrd::AcquisitionFlags::kIsNoiseMeasurement));
    EXPECT_EQ(acq_head2.measurement_uid, 42);
    EXPECT_EQ(acq_head2.idx.repetition, 1);
    EXPECT_EQ(acq_head2.center_sample, 64);
}

TEST_F(python_converter_test, array_mrd_acquisitionheader)
{
    {
        GILLock gl;     // this is needed
        register_converter<mrd::AcquisitionHeader>();

        boost::python::object main(boost::python::import("__main__"));
        boost::python::object global(main.attr("__dict__"));
        boost::python::exec("import mrd\n"
            "def mk_acq_headers(acq_head_array):\n"
            "   acq_head_array[8, 24].measurement_uid = 120\n"
            "   return acq_head_array\n",
            global, global);
    }

    hoNDArray<mrd::AcquisitionHeader> acq_head_array;
    acq_head_array.create(30, 10);
    for (int n = 0; n<acq_head_array.get_number_of_elements(); n++) {
        acq_head_array(n).measurement_uid = 345;
        acq_head_array(n).position = { 1.0f, 2.0f, 3.0f };
        acq_head_array(n).user_int = { 0, 1, 2, 3, 4 };
    }

    PythonFunction<hoNDArray<mrd::AcquisitionHeader> > mk_acq_headers("__main__", "mk_acq_headers");
    acq_head_array = mk_acq_headers(acq_head_array);

    EXPECT_EQ(acq_head_array(23, 8).measurement_uid, 345);
    EXPECT_EQ(acq_head_array(24, 7).measurement_uid, 345);
    EXPECT_EQ(acq_head_array(24, 9).measurement_uid, 345);
    EXPECT_EQ(acq_head_array(24, 8).measurement_uid, 120);

    EXPECT_FLOAT_EQ(acq_head_array(12, 4).position[0], 1.0f);
    EXPECT_FLOAT_EQ(acq_head_array(12, 4).position[1], 2.0f);
    EXPECT_FLOAT_EQ(acq_head_array(12, 4).position[2], 3.0f);

    EXPECT_EQ(acq_head_array(12, 4).user_int[3], 3);
}

TEST_F(python_converter_test, mrd_waveform)
{
    {
        GILLock gl;     // this is needed
        boost::python::object main(boost::python::import("__main__"));
        boost::python::object global(main.attr("__dict__"));
        boost::python::exec("def modify(head): head.measurement_uid += 1; return head",
            global, global);
    }

    mrd::WaveformUint32 wav1;
    wav1.measurement_uid = 41;

    size_t channels = 1;
    size_t number_of_samples = 12;
    wav1.data.create(number_of_samples, channels);
    for (size_t n = 0; n < channels * number_of_samples; n++) {
        wav1.data(n) = n;
    }

    PythonFunction<mrd::WaveformUint32> modify_wav_header("__main__", "modify");
    mrd::WaveformUint32 wav2 = modify_wav_header(wav1);
    EXPECT_EQ(wav2.measurement_uid, 42);
    EXPECT_EQ(wav2.data(5), 5);
}

TEST_F(python_converter_test, vec_mrd_waveform)
{
    {
        GILLock gl;     // this is needed
        register_converter<mrd::WaveformUint32>();

        boost::python::object main(boost::python::import("__main__"));
        boost::python::object global(main.attr("__dict__"));
        boost::python::exec("import mrd\n"
            "def mk_vector_waveform(wav): \n"
            "   for x in wav:\n"
            "       x.measurement_uid = 121\n"
            "   return wav\n",
            global, global);
    }

    size_t channels = 1;
    size_t number_of_samples = 12;
    std::vector<mrd::WaveformUint32> wav(3);
    for (size_t n = 0; n < wav.size(); n++)
    {
        wav[n].measurement_uid = 41;
        wav[n].data.create(number_of_samples, channels);
        for (size_t k = 0; k < channels * number_of_samples; k++) {
            wav[n].data(k) = k;
        }
    }

    PythonFunction<std::vector<mrd::WaveformUint32>> modify_vector_wav("__main__", "mk_vector_waveform");
    wav = modify_vector_wav(wav);
    EXPECT_EQ(wav[0].measurement_uid, 121);
    EXPECT_EQ(wav[1].measurement_uid, 121);
    EXPECT_EQ(wav[2].measurement_uid, 121);
    EXPECT_EQ(wav[1].data(5), 5);
}

TEST_F(python_converter_test, mrd_imageheader)
{
    {
        GILLock gl;     // this is needed
        boost::python::object main(boost::python::import("__main__"));
        boost::python::object global(main.attr("__dict__"));
        boost::python::exec("def modify(head): head.measurement_uid = 42; return head",
            global, global);
    }

    mrd::ImageHeader img_head;
    img_head.measurement_uid = 0;
    img_head.image_type = mrd::ImageType::kComplex;
    img_head.field_of_view = { 1.0f, 2.0f, 3.0f };

    PythonFunction<mrd::ImageHeader> modify_img_header("__main__", "modify");
    mrd::ImageHeader img_head2 = modify_img_header(img_head);
    EXPECT_EQ(img_head2.measurement_uid, 42);
    EXPECT_EQ(img_head2.image_type, mrd::ImageType::kComplex);
    EXPECT_FLOAT_EQ(img_head2.field_of_view[0], 1.0f);
    EXPECT_FLOAT_EQ(img_head2.field_of_view[1], 2.0f);
    EXPECT_FLOAT_EQ(img_head2.field_of_view[2], 3.0f);
}

TEST_F(python_converter_test, array_mrd_imageheader)
{
    {
        GILLock gl;     // this is needed
        register_converter<mrd::ImageHeader>();

        boost::python::object main(boost::python::import("__main__"));
        boost::python::object global(main.attr("__dict__"));
        boost::python::exec("import mrd\n"
            "def mk_image_headers(img_head_array): \n"
            "    img_head_array[8, 24].measurement_uid = 120\n"
            "    return img_head_array\n",
            global, global);
    }

    hoNDArray<mrd::ImageHeader> img_head_array;
    img_head_array.create(30, 10);
    for (int n = 0; n < img_head_array.get_number_of_elements(); n++) {
        img_head_array(n).measurement_uid = 345;
    }

    PythonFunction<hoNDArray<mrd::ImageHeader> > make_image_header("__main__", "mk_image_headers");
    img_head_array = make_image_header(img_head_array);
    EXPECT_EQ(img_head_array(23, 8).measurement_uid, 345);
    EXPECT_EQ(img_head_array(24, 7).measurement_uid, 345);
    EXPECT_EQ(img_head_array(24, 9).measurement_uid, 345);
    EXPECT_EQ(img_head_array(24, 8).measurement_uid, 120);
    EXPECT_EQ(img_head_array.get_size(0), 30);
    EXPECT_EQ(img_head_array.get_size(1), 10);
}

TEST_F(python_converter_test, mrd_imagemeta)
{
    {
        GILLock gl;     // this is needed
        boost::python::object main(boost::python::import("__main__"));
        boost::python::object global(main.attr("__dict__"));
        boost::python::exec("import mrd\n"
            "def modify_meta(meta): \n"
            "   meta['TestLong'].append(mrd.ImageMetaValue.Int64(4))\n"
            "   for i in range(3):\n"
            "       meta['TestDouble'][i].value += 0.1\n"
            "   meta['TestString'][2] = mrd.ImageMetaValue.Int64(42)\n"
            "   return meta\n",
            global, global);
    }

    mrd::ImageMeta meta;
    meta["TestLong"] = { 1L, 2L, 3L };

    meta["TestDouble"] = { 1.0, 2.1, 3.2 };

    meta["TestString"] = {"This", "is", "a test!"};

    PythonFunction<mrd::ImageMeta> modify_meta("__main__", "modify_meta");
    meta = modify_meta(meta);

    EXPECT_EQ(std::get<long>(meta["TestLong"][0]), 1);
    EXPECT_EQ(std::get<long>(meta["TestLong"][1]), 2);
    EXPECT_EQ(std::get<long>(meta["TestLong"][2]), 3);
    EXPECT_EQ(std::get<long>(meta["TestLong"][3]), 4);

    EXPECT_FLOAT_EQ(std::get<double>(meta["TestDouble"][0]), 1.1);
    EXPECT_FLOAT_EQ(std::get<double>(meta["TestDouble"][1]), 2.2);
    EXPECT_FLOAT_EQ(std::get<double>(meta["TestDouble"][2]), 3.3);

    EXPECT_EQ(std::get<std::string>(meta["TestString"][0]), "This");
    EXPECT_EQ(std::get<std::string>(meta["TestString"][1]), "is");
    EXPECT_EQ(std::get<long>(meta["TestString"][2]), 42);
}

TEST_F(python_converter_test, array_mrd_imagemeta)
{
    {
        GILLock gl;     // this is needed
        register_converter<mrd::ImageMeta>();

        boost::python::object main(boost::python::import("__main__"));
        boost::python::object global(main.attr("__dict__"));
        boost::python::exec("import mrd\n"
            "def mk_array_meta(metas): \n"
            "   metas[6, 11]['TestLong'].append(mrd.ImageMetaValue.Int64(4))\n"
            "   for i in range(3):\n"
            "      metas[5, 10]['TestDouble'][i].value *= 2\n"
            "   metas[4, 9]['TestString'][2] = mrd.ImageMetaValue.Int64(42)\n"
            "   return metas\n",
            global, global);
    }

    hoNDArray<mrd::ImageMeta> meta(12, 7);

    for (int n = 0; n < meta.size(); n++)
    {
        meta[n]["TestLong"] = { 1L, 2L, 3L };
        meta[n]["TestDouble"] = { 1.0, 2.1, 3.2 };
        meta[n]["TestString"] = { "This", "is", "a test!" };
    }

    PythonFunction< hoNDArray<mrd::ImageMeta> > mk_array_meta("__main__", "mk_array_meta");
    meta = mk_array_meta(meta);

    EXPECT_EQ(meta.get_size(0), 12);
    EXPECT_EQ(meta.get_size(1), 7);

    EXPECT_EQ(std::get<long>(meta(0, 0)["TestLong"][0]), 1);
    EXPECT_EQ(std::get<long>(meta(1, 1)["TestLong"][1]), 2);
    EXPECT_EQ(std::get<long>(meta(11, 6)["TestLong"][2]), 3);
    EXPECT_EQ(std::get<long>(meta(11, 6)["TestLong"][3]), 4);

    EXPECT_FLOAT_EQ(std::get<double>(meta(7, 5)["TestDouble"][0]), 1.0);
    EXPECT_FLOAT_EQ(std::get<double>(meta(10, 5)["TestDouble"][0]), 2.0);
    EXPECT_FLOAT_EQ(std::get<double>(meta(10, 5)["TestDouble"][1]), 4.2);
    EXPECT_FLOAT_EQ(std::get<double>(meta(10, 5)["TestDouble"][2]), 6.4);

    EXPECT_EQ(std::get<std::string>(meta(8, 4)["TestString"][2]), "a test!");
    EXPECT_EQ(std::get<std::string>(meta(9, 4)["TestString"][0]), "This");
    EXPECT_EQ(std::get<std::string>(meta(9, 4)["TestString"][1]), "is");
    EXPECT_EQ(std::get<long>(meta(9, 4)["TestString"][2]), 42);
}

TEST_F(python_converter_test, mrd_image_array)
{
    {
        GILLock gl;     // this is needed
        register_converter<mrd::ImageArray>();

        boost::python::object main(boost::python::import("__main__"));
        boost::python::object global(main.attr("__dict__"));
        boost::python::exec("import mrd\n"
            "def mk_mrd_image_array(array_data): \n"
            "   array_data.data *= 2\n"
            "   array_data.headers[1, 4, 3].measurement_uid = 12345\n"
            "   array_data.meta[1, 4, 3]['TestString'][0] = mrd.ImageMetaValue.String('Gadgetron')\n"
            "   array_data.waveforms[5].measurement_uid = 4253\n"
            "   return array_data\n",
            global, global);
    }

    mrd::ImageArray array_data;
    array_data.data.create(192, 144, 1, 32, 4, 5, 2); // [RO E1 E2 CHA N S SLC]
    array_data.headers.create(4, 5, 2);
    array_data.meta.create(4, 5, 2);
    array_data.waveforms = std::vector<mrd::WaveformUint32>(10);

    size_t n;
    for (n = 0; n < array_data.data.get_number_of_elements(); n++)
    {
        array_data.data(n) = std::complex<float>(3.0, 124.2);
    }

    for (n = 0; n < array_data.headers.get_number_of_elements(); n++)
    {
        array_data.headers(n).measurement_uid = 123;
    }

    for (int n = 0; n < array_data.meta.get_number_of_elements(); n++)
    {
        array_data.meta(n)["TestLong"] = { 1L * n, 2L * n, 3L * n };
        array_data.meta(n)["TestDouble"] = { 1.0 * n, 2.1 * n, 3.2 * n };
        array_data.meta(n)["TestString"] = { "This", "is", "a test!" };
    }

    for (n = 0; n < array_data.waveforms.size(); n++)
    {
        array_data.waveforms[n].measurement_uid = 42;
        array_data.waveforms[n].data.create(12);
        for (size_t k = 0; k < 12; k++) {
            array_data.waveforms[n].data[k] = k;
        }
    }

    PythonFunction< mrd::ImageArray > mk_mrd_image_array("__main__", "mk_mrd_image_array");
    array_data = mk_mrd_image_array(array_data);

    EXPECT_FLOAT_EQ(array_data.data(65558).real(), 6.0);
    EXPECT_FLOAT_EQ(array_data.data(65558).imag(), 248.4);

    EXPECT_EQ(array_data.headers(0, 2, 0).measurement_uid, 123);
    EXPECT_EQ(array_data.headers(1, 3, 1).measurement_uid, 123);
    EXPECT_EQ(array_data.headers(2, 1, 0).measurement_uid, 123);
    EXPECT_EQ(array_data.headers(3, 4, 1).measurement_uid, 12345);

    EXPECT_EQ(array_data.waveforms[0].measurement_uid, 42);
    EXPECT_EQ(array_data.waveforms[5].measurement_uid, 4253);

    EXPECT_EQ(std::get<std::string>(array_data.meta(0, 2, 0)["TestString"][0]), "This");
    EXPECT_EQ(std::get<std::string>(array_data.meta(1, 3, 1)["TestString"][1]), "is");
    EXPECT_EQ(std::get<std::string>(array_data.meta(2, 1, 0)["TestString"][2]), "a test!");
    EXPECT_EQ(std::get<std::string>(array_data.meta(3, 4, 1)["TestString"][0]), "Gadgetron");
}

TEST_F(python_converter_test, mrd_recon_data)
{
    {
        GILLock gl;     // this is needed
        register_converter<mrd::ReconData>();

        boost::python::object main(boost::python::import("__main__"));
        boost::python::object global(main.attr("__dict__"));
        boost::python::exec("import mrd\n"
            "def mk_mrd_recon_data(recon_data): \n"
            "    recon_data.buffers[0].data.data *= 2\n"
            "    recon_data.buffers[0].data.headers[1, 4, 3].measurement_uid = 12345\n"
            "    return recon_data\n",
            global, global);
    }

    mrd::ReconData recon_data;
    recon_data.buffers.resize(1);
    recon_data.buffers[0].data.data.create(192, 144, 1, 2, 4, 5, 2); // [RO E1 E2 CHA N S SLC]
    recon_data.buffers[0].data.headers.create(4, 5, 2);

    size_t n;
    for (n = 0; n < recon_data.buffers[0].data.data.get_number_of_elements(); n++)
    {
        recon_data.buffers[0].data.data(n) = std::complex<float>(3.0, 124.2);
    }

    for (n = 0; n < recon_data.buffers[0].data.headers.get_number_of_elements(); n++)
    {
        recon_data.buffers[0].data.headers(n).measurement_uid = 123;
    }

    PythonFunction< mrd::ReconData > mk_mrd_recon_data("__main__", "mk_mrd_recon_data");
    recon_data = mk_mrd_recon_data(recon_data);

    EXPECT_FLOAT_EQ(recon_data.buffers[0].data.data(65558).real(), 6.0);
    EXPECT_FLOAT_EQ(recon_data.buffers[0].data.data(65558).imag(), 248.4);

    EXPECT_EQ(recon_data.buffers[0].data.headers(0, 2, 0).measurement_uid, 123);
    EXPECT_EQ(recon_data.buffers[0].data.headers(1, 3, 1).measurement_uid, 123);
    EXPECT_EQ(recon_data.buffers[0].data.headers(2, 1, 0).measurement_uid, 123);
    EXPECT_EQ(recon_data.buffers[0].data.headers(3, 4, 1).measurement_uid, 12345);
}
