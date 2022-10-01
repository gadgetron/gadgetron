#include "python_toolbox.h"
#include "ismrmrd/ismrmrd.h"
#include "hoNDArray_utils.h"
#include "hoNDArray_elemwise.h"
#include "log.h"

using namespace Gadgetron;

int main(int argc, char** argv)
{
    initialize_python();
  
    char* gt_home = std::getenv("GADGETRON_HOME");
    std::string path_name;
    if (gt_home != NULL)
    {
        size_t pos = std::string(gt_home).rfind("gadgetron");
        gt_home[pos - 1] = '\0';
        path_name = std::string(gt_home) + std::string("/share/gadgetron/python");

        std::string add_path_cmd = std::string("import sys;\nsys.path.insert(0, \"") + path_name + std::string("\")\n");
        GDEBUG_STREAM(add_path_cmd);

        GILLock gl;
        boost::python::exec(add_path_cmd.c_str(),
            boost::python::import("__main__").attr("__dict__"));
    }


    {
        GDEBUG_STREAM(" --------------------------------------------------------------------------------------------------");
        GDEBUG_STREAM("Call a function with no return value (print all arguments)");

        int a = -42;
        float b = 3.141592;
        std::string c("hello, world");
        unsigned int d(117);
        std::complex<double> e(2.12894, -1.103103);

        std::vector<size_t> dims;
        dims.push_back(4);
        dims.push_back(4);
        dims.push_back(4);
        hoNDArray<std::complex<float> > arr(dims);

#if PY_MAJOR_VERSION == 3
        PythonFunction<> foo("builtins", "print");
#else
        PythonFunction<> foo("__builtin__", "print");
#endif
        foo(a, b, c, d, e, arr);
    }

    {
        GDEBUG_STREAM(" --------------------------------------------------------------------------------------------------");
        GDEBUG_STREAM("Call a function with a single return value");
        PythonFunction<float> atan2("math", "atan2");
        int x = 7, y = 4;
        float atan = atan2(x, y);
        GINFO_STREAM(std::to_string(atan) << std::endl);
    }

    {
        GDEBUG_STREAM(" --------------------------------------------------------------------------------------------------");
        GDEBUG_STREAM("Call a function that returns a tuple");

#if PY_MAJOR_VERSION == 3
        PythonFunction<float, float> divmod("builtins", "divmod");
#else
        PythonFunction<float, float> divmod("__builtin__", "divmod");
#endif
        float w = 6.89;
        float z = 4.12;
        float fsum = 0, fdiff = 0;
        std::tie(fsum, fdiff) = divmod(w, z);
        GINFO_STREAM(fsum << ", " << fdiff << std::endl);
    }

    {
        GDEBUG_STREAM(" --------------------------------------------------------------------------------------------------");
        GDEBUG_STREAM("Call a function that expects an iterable argument (tuple)");
        
#if PY_MAJOR_VERSION == 3
        PythonFunction<int> tuplen("builtins", "len");
#else
        PythonFunction<int> tuplen("__builtin__", "len");
#endif
        int l = tuplen(std::make_tuple(-7, 0, 7));
        GINFO_STREAM("tuple length: " << l << std::endl);
    }

    {
        GDEBUG_STREAM(" --------------------------------------------------------------------------------------------------");
        GDEBUG_STREAM("Generate an hoNDArray of even #s using numpy");
        PythonFunction<hoNDArray<float>> arange("numpy", "arange");
        // hoNDArray<float> evens = arange(0, 100, 2, "float64");
        hoNDArray<float> evens = arange(0, 100, 2, "f");
        GINFO_STREAM("number of even numbers between 0 and 100: " << std::to_string(evens.get_number_of_elements()) << std::endl);
    }

    {
        GDEBUG_STREAM(" --------------------------------------------------------------------------------------------------");
        GDEBUG_STREAM("Test return mutltiple hoNDArray arrays");
        {
            GILLock gl;     // this is needed
            boost::python::object main(boost::python::import("__main__"));
            boost::python::object global(main.attr("__dict__"));
            boost::python::exec("import numpy as np\n"
                "def scale_array(a): \n"
                "   print(len(a))\n"
                "   b = a + 100\n"
                "   return a,b\n",
                global, global);
        }

        hoNDArray<float> a, a2, b;
        a.create(32, 64);
        Gadgetron::fill(a, float(45));

        PythonFunction< hoNDArray<float>, hoNDArray<float> > scale_array_test("__main__", "scale_array");
        std::tie(a2, b) = scale_array_test(a);
        GINFO_STREAM("a2[0] = " << a2[0] << std::endl);
        GINFO_STREAM("b[12] = " << b[12] << std::endl);
        GINFO_STREAM(std::endl);
    }

    {
        GDEBUG_STREAM(" --------------------------------------------------------------------------------------------------");
        GDEBUG_STREAM("Test return mutltiple hoNDArray arrays reutrn and multi-inputs");
        {
            GILLock gl;     // this is needed
            boost::python::object main(boost::python::import("__main__"));
            boost::python::object global(main.attr("__dict__"));
            boost::python::exec("import numpy as np\n"
                "def scale_array(a, b): \n"
                "   print(len(a))\n"
                "   print(len(b))\n"
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
        GINFO_STREAM("a2[0]  = " << a2[0] << std::endl);
        GINFO_STREAM("b2[12] = " << b[12] << std::endl);
        GINFO_STREAM("c[20]  = " << c[20] << std::endl);
        GINFO_STREAM(std::endl);
    }

    {
        {
            GILLock gl;     // this is needed
            boost::python::object main(boost::python::import("__main__"));
            boost::python::object global(main.attr("__dict__"));
            boost::python::exec("def modify(head): head.version = head.version+1; return head",
                global, global);
        }


        GDEBUG_STREAM(" --------------------------------------------------------------------------------------------------");
        GDEBUG_STREAM("Test converter for ISMRMRD::ISMRMRD_WaveformHeader");
        ISMRMRD::ISMRMRD_WaveformHeader wav_head, wav_head2;
        wav_head.version = 41;
        GINFO_STREAM("version before: " << wav_head.version << std::endl);
        PythonFunction<ISMRMRD::ISMRMRD_WaveformHeader> modify_wav_header("__main__", "modify");
        wav_head2 = modify_wav_header(wav_head);
        GINFO_STREAM("version after: " << wav_head2.version << std::endl);
    }

    {
        {
            GILLock gl;     // this is needed
            boost::python::object main(boost::python::import("__main__"));
            boost::python::object global(main.attr("__dict__"));
            boost::python::exec("def modify(head): head.version = head.version+1; return head",
                global, global);
        }


        GDEBUG_STREAM(" --------------------------------------------------------------------------------------------------");
        GDEBUG_STREAM("Test converter for ISMRMRD::Waveform");
        ISMRMRD::Waveform wav1, wav2;
        wav1.head.version = 41;
        wav1.head.channels = 1;
        wav1.head.number_of_samples = 12;
        wav1.data = new uint32_t[wav1.head.channels*wav1.head.number_of_samples];
        for (size_t n = 0; n < wav1.head.channels*wav1.head.number_of_samples; n++)
            wav1.data[n] = n;

        GINFO_STREAM("version before: " << wav1.head.version << std::endl);
        PythonFunction<ISMRMRD::Waveform> modify_wav_header("__main__", "modify");
        wav2 = modify_wav_header(wav1);
        GINFO_STREAM("version after: " << wav2.head.version << std::endl);
        GINFO_STREAM("contents: " << std::endl);
        for (size_t n = 0; n < wav1.head.channels*wav1.head.number_of_samples; n++)
            GINFO_STREAM(wav2.data[n] << " ");
        GINFO_STREAM(std::endl);
    }

    {
        GDEBUG_STREAM(" --------------------------------------------------------------------------------------------------");
        GDEBUG_STREAM("Test converter for std::vector<ISMRMRD::Waveform>");
        {
            GILLock gl;     // this is needed
            boost::python::object main(boost::python::import("__main__"));
            boost::python::object global(main.attr("__dict__"));
            boost::python::exec("import ismrmrd\n"
                "def mk_vector_waveform(wav): \n"
                "   print(len(wav))\n"
                "   for x in wav:\n"
                "       print(x.version)\n"
                "       x.version=121\n"
                "   return wav\n",
                global, global);
        }

        std::vector<ISMRMRD::Waveform> wav(3);
        for (size_t n = 0; n < wav.size(); n++)
        {
            wav[n].head.version = 41;
            wav[n].head.channels = 1;
            wav[n].head.number_of_samples = 12;
            wav[n].data = new uint32_t[wav[n].head.channels*wav[n].head.number_of_samples];
            for (size_t k = 0; k < wav[n].head.channels*wav[n].head.number_of_samples; k++)
                wav[n].data[k] = k;
        }

        PythonFunction<std::vector<ISMRMRD::Waveform>> modify_vector_wav("__main__", "mk_vector_waveform");
        wav = modify_vector_wav(wav);
        GINFO_STREAM("version after: " << wav[0].head.version << std::endl);
        GINFO_STREAM("contents: " << std::endl);
        for (size_t n = 0; n < wav[0].head.channels*wav[0].head.number_of_samples; n++)
            GINFO_STREAM(wav[0].data[n] << " ");
        GINFO_STREAM(std::endl);
    }

    {
        {
            GILLock gl;     // this is needed
            boost::python::object main(boost::python::import("__main__"));
            boost::python::object global(main.attr("__dict__"));
            boost::python::exec("def modify(head): head.version = 42; return head",
                    global, global);
        }

        GDEBUG_STREAM(" --------------------------------------------------------------------------------------------------");
        GDEBUG_STREAM("Test converter for ISMRMRD::ImageHeader");
        ISMRMRD::ImageHeader img_head, img_head2;
        img_head.version = 0;
        GINFO_STREAM("version before: " << img_head.version << std::endl);
        PythonFunction<ISMRMRD::ImageHeader> modify_img_header("__main__", "modify");
        img_head2 = modify_img_header(img_head);
        GINFO_STREAM("version after: " << img_head2.version << std::endl);

        GDEBUG_STREAM(" --------------------------------------------------------------------------------------------------");
        GDEBUG_STREAM("Test converter for ISMRMRD::AcquisitionHeader");
        ISMRMRD::AcquisitionHeader acq_head, acq_head2;
        acq_head.version = 0;
        GINFO_STREAM("version before: " << img_head.version << std::endl);
        PythonFunction<ISMRMRD::AcquisitionHeader> modify_acq_header("__main__", "modify");
        acq_head2 = modify_acq_header(acq_head);
        GINFO_STREAM("version after: " << acq_head2.version << std::endl);
    }

    {
        GDEBUG_STREAM(" --------------------------------------------------------------------------------------------------");
        GDEBUG_STREAM("Test converter for std::vector<std::complex<float>>");
        {
            GILLock gl;     // this is needed
            boost::python::object main(boost::python::import("__main__"));
            boost::python::object global(main.attr("__dict__"));
            boost::python::exec("from numpy.random import random\n"
                    "def rand_cplx_array(length): \n"
                    "    return random(length) + 1j * random(length)\n",
                    global, global);
        }

        std::vector<std::complex<double> > vec;
        PythonFunction<std::vector<std::complex<double> > > make_vec("__main__", "rand_cplx_array");
        vec = make_vec(32);
        GINFO_STREAM(vec[16] << std::endl);
    }

    {
        GDEBUG_STREAM(" --------------------------------------------------------------------------------------------------");
        GDEBUG_STREAM("Test for hoNDArray<ISMRMRD::ImageHeader>")
        {
            GILLock gl;     // this is needed
            boost::python::object main(boost::python::import("__main__"));
            boost::python::object global(main.attr("__dict__"));
            boost::python::exec("import ismrmrd\n"
                "def mk_image_headers(img_head_array): \n"
                "   img_head_array[2,4].version=120\n"
                "   print(img_head_array[0,0])\n"
                "   print(img_head_array[2,4])\n"
                "   return img_head_array\n",
                global, global);
        }

        hoNDArray<ISMRMRD::ImageHeader> img_head_array;
        img_head_array.create(30, 10);
        for (int n=0; n<img_head_array.get_number_of_elements(); n++)
            img_head_array(n).version = 345;

        GDEBUG_STREAM(" --------------------------------------------------------------------------------------------------");
        GDEBUG_STREAM("Test converter for PythonFunction<hoNDArray<ISMRMRD::ImageHeader> >");
        PythonFunction<hoNDArray<ISMRMRD::ImageHeader> > make_image_header("__main__", "mk_image_headers");
        img_head_array = make_image_header(img_head_array);
        GINFO_STREAM(img_head_array(2, 4).version << std::endl);

        GDEBUG_STREAM(" --------------------------------------------------------------------------------------------------");
        GDEBUG_STREAM("Test for ISMRMRD::MetaContainer")
        {
            GILLock gl;     // this is needed
            boost::python::object main(boost::python::import("__main__"));
            boost::python::object global(main.attr("__dict__"));
            boost::python::exec("import ismrmrd\n"
                "def mk_meta(meta): \n"
                "   mt = ismrmrd.Meta.deserialize(meta)\n"
                "   print(mt['TestLong'])\n"
                "   print(mt['TestDouble'])\n"
                "   print(mt['TestString'])\n"
                "   mt_str = ismrmrd.Meta.serialize(mt)\n"
                "   return mt_str\n",
                global, global);
        }
    }

    {
        GDEBUG_STREAM(" --------------------------------------------------------------------------------------------------");
        GDEBUG_STREAM("Test for hoNDArray<ISMRMRD::AcquisitionHeader>")
        {
            GILLock gl;     // this is needed
            boost::python::object main(boost::python::import("__main__"));
            boost::python::object global(main.attr("__dict__"));
            boost::python::exec("import ismrmrd\n"
                "def mk_acq_headers(acq_head_array): \n"
                "   acq_head_array[2,4].version=120\n"
                "   print(acq_head_array[0,0])\n"
                "   print(acq_head_array[2,4])\n"
                "   return acq_head_array\n",
                global, global);
        }

        hoNDArray<ISMRMRD::AcquisitionHeader> acq_head_array;
        acq_head_array.create(30, 10);
        for (int n = 0; n<acq_head_array.get_number_of_elements(); n++)
            acq_head_array(n).version = 345;

        GDEBUG_STREAM(" --------------------------------------------------------------------------------------------------");
        GDEBUG_STREAM("Test converter for PythonFunction<hoNDArray<ISMRMRD::AcquisitionHeader> >");
        PythonFunction<hoNDArray<ISMRMRD::AcquisitionHeader> > mk_acq_headers("__main__", "mk_acq_headers");
        acq_head_array = mk_acq_headers(acq_head_array);
        GINFO_STREAM(acq_head_array(2, 4).version << std::endl);
    }

    {
        ISMRMRD::MetaContainer meta;
        meta.set("TestLong", (long)1);
        meta.append("TestLong", (long)2);
        meta.append("TestLong", (long)3);

        meta.set("TestDouble", 1.0);
        meta.append("TestDouble", 2.1);
        meta.append("TestDouble", 3.2);

        meta.set("TestString", "This");
        meta.append("TestString", "is");
        meta.append("TestString", "a test!");

        PythonFunction<ISMRMRD::MetaContainer> make_meta("__main__", "mk_meta");
        ISMRMRD::MetaContainer meta_res = make_meta(meta);
        std::stringstream meta_res_str;
        ISMRMRD::serialize(meta_res, meta_res_str);
        GDEBUG_STREAM(meta_res_str.str());
    }

    {
        GDEBUG_STREAM(" --------------------------------------------------------------------------------------------------");
        GDEBUG_STREAM("Test converter for std::vector<ISMRMRD::MetaContainer>");
        {
            GILLock gl;     // this is needed
            boost::python::object main(boost::python::import("__main__"));
            boost::python::object global(main.attr("__dict__"));
            boost::python::exec("import ismrmrd\n"
                "def mk_vector_meta(meta_str): \n"
                "   print(len(meta_str))\n"
                "   mt = list()\n"
                "   for x in meta_str:\n"
                "       mt.append(ismrmrd.Meta.deserialize(x))\n"
                "   print(mt[0]['TestLong'])\n"
                "   print(mt[0]['TestDouble'])\n"
                "   print(mt[0]['TestString'])\n"
                "   print(mt[11]['TestLong'])\n"
                "   print(mt[11]['TestDouble'])\n"
                "   print(mt[11]['TestString'])\n"
                "   res_str = list()\n"
                "   for x in mt:\n"
                "       res_str_curr=ismrmrd.Meta.serialize(x)\n"
                "       res_str.append(res_str_curr)\n"
                "   return res_str\n",
                global, global);
        }

        std::vector<ISMRMRD::MetaContainer> meta(12);

        for (int n=0; n<meta.size(); n++)
        {
            meta[n].set("TestLong", (long)1 * n);
            meta[n].append("TestLong", (long)2 * n);
            meta[n].append("TestLong", (long)3 * n);

            meta[n].set("TestDouble", 1.0 * n);
            meta[n].append("TestDouble", 2.1 * n);
            meta[n].append("TestDouble", 3.2 * n);

            meta[n].set("TestString", "This");
            meta[n].append("TestString", "is");
            meta[n].append("TestString", "a test!");
        }

        PythonFunction< std::vector<ISMRMRD::MetaContainer> > mk_vector_meta("__main__", "mk_vector_meta");
        std::vector<ISMRMRD::MetaContainer> meta_res = mk_vector_meta(meta);

        for (int n = 0; n<meta.size(); n++)
        {
            GDEBUG_STREAM("Meta data : " << n);
            GDEBUG_STREAM("-------------------------------------------------");
            std::stringstream meta_res_str;
            ISMRMRD::serialize(meta_res[n], meta_res_str);
            GDEBUG_STREAM(meta_res_str.str());
        }
    }

    if (gt_home != NULL)
    {
        GDEBUG_STREAM(" --------------------------------------------------------------------------------------------------");
        GDEBUG_STREAM("Test converter for ISDMRMRD::IsmrmrdImageArray");

        {
            GILLock gl;     // this is needed
            boost::python::object main(boost::python::import("__main__"));
            boost::python::object global(main.attr("__dict__"));
            boost::python::exec("import ismrmrd\n"
                "def mk_ismrmrd_image_array(array_data): \n"
                "   print(array_data.data.shape)\n"
                "   print(array_data.data[128, 56, 0, 12, 3, 4, 1])\n"
                "   print(array_data.headers[3, 4, 0])\n"
                "   print(array_data.acq_headers[3, 4, 0])\n"
                "   mt = list()\n"
                "   for x in array_data.meta:\n"
                "       curr_meta = ismrmrd.Meta.deserialize(x)\n"
                "       curr_meta['TestString'][0]='Gadgetron'\n"
                "       mt.append(curr_meta)\n"
                "   array_data.headers[1, 2, 0].version=12345\n"
                "   res_str = list()\n"
                "   for x in mt:\n"
                "       res_str_curr=ismrmrd.Meta.serialize(x)\n"
                "       res_str.append(res_str_curr)\n"
                "   array_data.meta=res_str\n"
                "   return array_data\n",
                global, global);
        }

        Gadgetron::IsmrmrdImageArray array_data;
        array_data.data_.create(192, 144, 1, 32, 4, 5, 2); // [RO E1 E2 CHA N S SLC]
        array_data.headers_.create(4, 5, 2);
        array_data.meta_.resize(4*5*2);
        array_data.waveform_ = std::vector<ISMRMRD::Waveform>(10);
        array_data.acq_headers_ = hoNDArray<ISMRMRD::AcquisitionHeader>(4, 5, 2);

        size_t n;
        for (n=0; n<array_data.data_.get_number_of_elements(); n++)
        {
            array_data.data_(n) = std::complex<float>(3.0, 124.2);
        }

        memset(array_data.headers_.get_data_ptr(), 0, sizeof(ISMRMRD::ImageHeader)*8);

        for (n = 0; n<array_data.headers_.get_number_of_elements(); n++)
        {
            array_data.headers_(n).version = 123;
            (*array_data.acq_headers_)(n).version = 12300;
        }

        for (int n=0; n<4 * 5 * 2; n++)
        {
            array_data.meta_[n].set("TestLong", (long)1 * n);
            array_data.meta_[n].append("TestLong", (long)2 * n);
            array_data.meta_[n].append("TestLong", (long)3 * n);

            array_data.meta_[n].set("TestDouble", 1.0 * n);
            array_data.meta_[n].append("TestDouble", 2.1 * n);
            array_data.meta_[n].append("TestDouble", 3.2 * n);

            array_data.meta_[n].set("TestString", "This");
            array_data.meta_[n].append("TestString", "is");
            array_data.meta_[n].append("TestString", "a test!");
        }

        for (n = 0; n < 10; n++)
        {
            (*array_data.waveform_)[n].head.version = 42;
            (*array_data.waveform_)[n].head.channels = 1;
            (*array_data.waveform_)[n].head.number_of_samples = 12;
            (*array_data.waveform_)[n].data = new uint32_t[12];
            for (size_t k = 0; k<12; k++)
            {
                (*array_data.waveform_)[n].data[k] = k;
            }
        }

        PythonFunction< Gadgetron::IsmrmrdImageArray > mk_ismrmrd_image_array("__main__", "mk_ismrmrd_image_array");
        Gadgetron::IsmrmrdImageArray array_res = mk_ismrmrd_image_array(array_data);

        GDEBUG_STREAM(array_data.data_(65558));
        GDEBUG_STREAM(array_data.headers_(2, 2, 0).version);
        GDEBUG_STREAM(array_data.headers_(1, 2, 0).version);
        GDEBUG_STREAM( (*array_data.acq_headers_)(1, 2, 0).version);

        std::stringstream meta_res_str;
        ISMRMRD::serialize(array_res.meta_[6], meta_res_str);
        GDEBUG_STREAM(meta_res_str.str());

        GDEBUG_STREAM( (*array_data.waveform_)[1].head.version);
    }

    if (gt_home != NULL)
    {
        GDEBUG_STREAM(" --------------------------------------------------------------------------------------------------");
        GDEBUG_STREAM("Test converter for ISDMRMRD::IsmrmrdReconData");

        {
            GILLock gl;     // this is needed
            boost::python::object main(boost::python::import("__main__"));
            boost::python::object global(main.attr("__dict__"));
            boost::python::exec("import ismrmrd\n"
                "def mk_ismrmrd_recon_data(array_data): \n"
                "   print(array_data[0].data.data.shape)\n"
                "   print(array_data[0].data.headers[3, 4, 0])\n"
                "   array_data[0].data.headers[3, 4, 0].version=345\n"
                "   return array_data\n",
                global, global);
        }

        Gadgetron::IsmrmrdReconData array_data;
        array_data.rbit_.resize(1);
        array_data.rbit_[0].data_.data_.create(192, 144, 1, 2, 4, 5, 2); // [RO E1 E2 CHA N S SLC]
        array_data.rbit_[0].data_.headers_.create(4, 5, 2);

        size_t n;
        for (n = 0; n<array_data.rbit_[0].data_.data_.get_number_of_elements(); n++)
        {
            array_data.rbit_[0].data_.data_(n) = std::complex<float>(3.0, 124.2);
        }

        memset(array_data.rbit_[0].data_.headers_.get_data_ptr(), 0, sizeof(ISMRMRD::AcquisitionHeader) * 8);

        for (n = 0; n<array_data.rbit_[0].data_.headers_.get_number_of_elements(); n++)
        {
            array_data.rbit_[0].data_.headers_(n).version = 123;
        }

        PythonFunction< Gadgetron::IsmrmrdReconData > mk_ismrmrd_recon_data("__main__", "mk_ismrmrd_recon_data");
        Gadgetron::IsmrmrdReconData array_res = mk_ismrmrd_recon_data(array_data);

        GDEBUG_STREAM(array_data.rbit_[0].data_.data_(65558));
        GDEBUG_STREAM(array_data.rbit_[0].data_.headers_(2, 2, 0).version);
        GDEBUG_STREAM(array_data.rbit_[0].data_.headers_(1, 2, 0).version);
        GDEBUG_STREAM(array_data.rbit_[0].data_.headers_(3, 4, 0).version);
    }

    return 0;
}
