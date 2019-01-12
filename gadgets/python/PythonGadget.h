#pragma once


#include "Node.h"
#include "hoNDArray.h"
#include "GadgetReference.h"
#include "gadgetronpython_export.h"
#include "python_toolbox.h"


#include <ismrmrd/ismrmrd.h>
#include <ismrmrd/meta.h>
#include <boost/python.hpp>
#include <boost/filesystem.hpp>

namespace Gadgetron {

    namespace Python {


        template<class... BASETYPES>
        using PythonTypePattern = Core::variant<IsmrmrdReconData,
                IsmrmrdImageArray,
                Core::tuple<ISMRMRD::WaveformHeader, hoNDArray<uint32_t>>,
                Core::tuple<ISMRMRD::AcquisitionHeader, Core::optional<hoNDArray<float>>, hoNDArray<std::complex<float>>>,
                Core::tuple<ISMRMRD::ImageHeader, hoNDArray<BASETYPES>, Core::optional<ISMRMRD::MetaContainer>>...
        >;

        using PythonTypes = PythonTypePattern<float, std::complex<float>, double, std::complex<double>, uint16_t, int16_t, uint32_t, int32_t>;
    }
    /// This PythonGadget is the gateway for c++ to call python
    class EXPORTGADGETSPYTHON PythonGadget : public Core::TypedGadgetNode<Python::PythonTypes> {
    public:

        PythonGadget(const Core::Context &context, const Core::GadgetProperties &params);

        NODE_PROPERTY(error_ignored_mode, bool, "If true failure of this python gadget will not stop the entire chain",false);


    protected:

        bool config_success_;
        NODE_PROPERTY(python_module, std::string, "Python module containing the Python Gadget class to be loaded", "");
        NODE_PROPERTY(python_class, std::string, "Python class to load from python module", "");
        NODE_PROPERTY(python_path, std::string, "Path(s) to add to the to the Python search path", "");

    private:
        boost::python::object module_;
        boost::python::object class_;

        int process_image(GadgetContainerMessage<ISMRMRD::ImageHeader> *hmi);

        std::map<std::string, std::string> parameters_python_;

        void register_converters();

        void process(Core::TypedInputChannel<Python::PythonTypes> &in, Core::OutputChannel &out) override;
    };
}
