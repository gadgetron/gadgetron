#include "PythonGadget.h"

namespace Gadgetron {

    namespace {

        template<class T>
        void process_message(boost::python::object &class_, Core::OutputChannel &out, T &message) {
            boost::python::object process_fn = class_.attr("process");
            process_fn(message);

        };

        void process_message(boost::python::object &class_, Core::OutputChannel &out,
                             Core::tuple<ISMRMRD::WaveformHeader, hoNDArray<uint32_t>> &waveform) {

            auto &head = std::get<0>(waveform);
            auto &data = std::get<1>(waveform);

            if (boost::python::hasattr(class_, "process_waveform")) {
                boost::python::object process_fn = class_.attr("process_waveform");
                process_fn(head, data);
            } else {
                out.push(std::make_unique<ISMRMRD::WaveformHeader>(head),
                         std::make_unique<hoNDArray<uint32_t>>(std::move(data)));
            }
        };

        void process_message(boost::python::object &class_, Core::OutputChannel &out,
                             Core::tuple<ISMRMRD::AcquisitionHeader, hoNDArray<std::complex<float>>, Core::optional<hoNDArray<float>>> &acquisition) {

            auto &head = std::get<0>(acquisition);
            auto &data = std::get<1>(acquisition);
            auto &traj = std::get<2>(acquisition);

            boost::python::object process_fn = class_.attr("process");
            if (traj) {
                process_fn(head, data, traj);
            } else {
                process_fn(head, data);
            }

        };

        template<class T>
        void process_image(boost::python::object &class_, Core::OutputChannel &out, ISMRMRD::ImageHeader &head,
                           hoNDArray<T> &data, Core::optional<ISMRMRD::MetaContainer> &meta) {

            boost::python::object process_fn = class_.attr("process");
            if (meta) {
                std::stringstream str;
                ISMRMRD::serialize(*meta, str);
                process_fn(head, data, str.str());
            } else {
                process_fn(head, data);
            }

        }

        void process_message(boost::python::object &class_, Core::OutputChannel &out,
                             Core::tuple<ISMRMRD::ImageHeader, Python::IsmrmrdArrayTypes, Core::optional<ISMRMRD::MetaContainer>> &image) {
            auto &head = std::get<0>(image);
            auto &data = std::get<1>(image);
            auto &meta = std::get<2>(image);
            boost::apply_visitor([&](auto array) {
                process_image(class_, out, head, array, meta);
            }, data);

        }

    }

    void PythonGadget::process(Core::TypedInputChannel<Python::PythonTypes> &in, Core::OutputChannel &out) {

        if (! config_success_) return;

        auto ref = std::make_shared<GadgetReference>(out);
        boost::python::object set_next_gadget = class_.attr("set_next_gadget");


//        set_next_gadget(boost::python::object());

        for (auto message : in) {
            GILLock lock;
            try {
                boost::apply_visitor([this, &out](auto &&message) { process_message(class_, out, message); }, *message);
            }
            catch (boost::python::error_already_set const &) {
                GDEBUG("Passing data on to python module failed\n");
                std::string err = pyerr_to_string();
                GERROR(err.c_str());
                if (!error_ignored_mode) {
                    throw std::runtime_error(err);
                }
            }
        }

    }


    PythonGadget::PythonGadget(const Core::Context &context, const Core::GadgetProperties &params)
            : Core::TypedGadgetNode<Python::PythonTypes>(params) {

        initialize_python();

        // start python interpreter
        register_converters();


        GDEBUG("Python Path            : %s\n", python_path.c_str());
        GDEBUG("Python Module          : %s\n", python_module.c_str());
        GDEBUG("Python Class           : %s\n", python_class.c_str());

        boost::filesystem::path gadgetron_python_path = context.paths.gadgetron_home / "share" / "gadgetron" / "python";
        GDEBUG("Python folder          : %s\n", gadgetron_python_path.generic_string().c_str());

        for (std::string path : {python_path, gadgetron_python_path.generic_string()}) {
            add_python_path(path);
        }

        if (std::empty(python_module))
            throw std::runtime_error("Empty python module provided");
        if (std::empty(python_class))
            throw std::runtime_error("Empty python class provided");

        GILLock lock;
//        try {
            module_ = boost::python::import(python_module.c_str());
            GDEBUG_STREAM("Successfully import module : " << python_module)

            // Reload the module so changes take place at Gadgetron runtime
            boost::python::import("__main__").attr("__dict__")[python_module. c_str()] = module_;
            GDEBUG_STREAM("Successfully set module : " << python_module)

            std::string tmp = std::string("reload(") +
                              std::string(python_module.c_str()) +
                              std::string(")\n");
#if defined PYVER && PYVER == 3
            // prefix reload call for Python 3
            tmp = std::string("from importlib import reload;") + tmp;
#endif
            // GDEBUG("Reloading with command: %s\n", tmp.c_str());
            boost::python::exec(tmp.c_str(), boost::python::import("__main__").attr("__dict__"));
            GDEBUG_STREAM("Successfully reload module : " << python_module)
            // Create instance of class (passing gadget reference)
            class_ = module_.attr(python_class.c_str())();
            // Increment reference count of Python class so that both the C++
            // destructor and the interpreter can decrement its reference count
            boost::python::incref(class_.ptr());
            GDEBUG_STREAM("Successfully declare class : " << python_class)


            boost::python::object set_parameter_fn = class_.attr(
                    "set_parameter");        // Transfer all properties/parameters to Python gadget
            for (auto &property : params) {
                set_parameter_fn(property.first, property.second);
            }

            // retrieve and call python gadget's process_config method
            GDEBUG_STREAM("Call process_config ... ")

            auto sstream = std::stringstream();
            ISMRMRD::serialize(context.header, sstream);

            boost::python::object process_config_fn =
                    class_.attr("process_config");

            boost::python::object ignored =
                    process_config_fn(boost::python::object(sstream.str()));


            config_success_ = true;
//        } catch (boost::python::error_already_set const &) {
//            GERROR("Error loading python modules in Gadget ");
//            std::string err = pyerr_to_string();
//            GERROR(err.c_str());
//            if (!error_ignored_mode)
//                throw std::runtime_error("Python " + err);
//        }
        GDEBUG_STREAM("Call process_config completed without error ")

    }

    void
    PythonGadget::register_converters() {// ensure boost can convert between hoNDArrays and NumPy arrays automatically
        register_converter<hoNDArray<std::complex<float> > >();
        register_converter<hoNDArray<float> >();
        register_converter<hoNDArray<unsigned short> >();
        register_converter<hoNDArray<uint32_t> >();
        register_converter<hoNDArray<ISMRMRD::AcquisitionHeader> >();
        register_converter<hoNDArray<ISMRMRD::ImageHeader> >();

        // ensure boost can convert ISMRMRD headers automatically
        register_converter<ISMRMRD::ImageHeader>();
        register_converter<ISMRMRD::AcquisitionHeader>();
        register_converter<ISMRMRD::ISMRMRD_WaveformHeader>();
        register_converter<ISMRMRD::Waveform>();

        register_converter<std::vector<std::complex<float> > >();
        register_converter<std::vector<float> >();
        register_converter<std::vector<unsigned short> >();
        register_converter<std::vector<ISMRMRD::MetaContainer> >();
        register_converter<std::vector<ISMRMRD::Waveform> >();

        register_converter<IsmrmrdReconData>();
        register_converter<IsmrmrdImageArray>();
    }


    GADGETRON_GADGET_EXPORT(PythonGadget)
}

