#include "Gadget.h"
#include "GadgetReference.h"
#include "GadgetContainerMessage.h"
#include "hoNDArray.h"
#include <ismrmrd/ismrmrd.h>
#include <ismrmrd/meta.h>

/* #include <boost/preprocessor/stringize.hpp> */
#include <boost/python.hpp>

namespace Gadgetron {


    void GadgetReference::return_recondata(boost::python::object rec) {
        output.push(IsmrmrdReconData(boost::python::extract<IsmrmrdReconData>(rec)()));

    }

    void GadgetReference::return_ismrmrd_image_array(boost::python::object rec) {
        output.push(IsmrmrdImageArray(boost::python::extract<IsmrmrdImageArray>(rec)()));
    }

    template<class TH, class TD>
    void GadgetReference::return_data(TH header, boost::python::object arr, const char *meta) {

        hoNDArray<TD> data = boost::python::extract<hoNDArray<TD>>(arr);

        if (meta) {
            auto m3 = ISMRMRD::MetaContainer{};
            ISMRMRD::deserialize(meta, m3);
            output.push(std::move(header), std::move(data), std::move(m3));
            return;
        }

        output.push(std::move(header), std::move(data));
    }

    void GadgetReference::return_acquisition(ISMRMRD::AcquisitionHeader acq, boost::python::object arr) {
        return_data<ISMRMRD::AcquisitionHeader, std::complex<float> >(acq, arr, 0);
    }

    void GadgetReference::return_image_cplx(ISMRMRD::ImageHeader img, boost::python::object arr) {
        return_data<ISMRMRD::ImageHeader, std::complex<float> >(img, arr, 0);
    }

    void
    GadgetReference::return_image_cplx_attr(ISMRMRD::ImageHeader img, boost::python::object arr, const char *meta) {
        return_data<ISMRMRD::ImageHeader, std::complex<float> >(img, arr, meta);
    }


    void GadgetReference::return_image_float(ISMRMRD::ImageHeader img, boost::python::object arr) {
        return_data<ISMRMRD::ImageHeader, float>(img, arr, 0);
    }

    void
    GadgetReference::return_image_float_attr(ISMRMRD::ImageHeader img, boost::python::object arr, const char *meta) {
        return_data<ISMRMRD::ImageHeader, float>(img, arr, meta);
    }

    void GadgetReference::return_image_ushort(ISMRMRD::ImageHeader img, boost::python::object arr) {
        return_data<ISMRMRD::ImageHeader, unsigned short>(img, arr, 0);
    }

    void
    GadgetReference::return_image_ushort_attr(ISMRMRD::ImageHeader img, boost::python::object arr, const char *meta) {
        return_data<ISMRMRD::ImageHeader, unsigned short>(img, arr, meta);
    }

    GadgetReference::GadgetReference(Core::OutputChannel &channel) : output(channel) {

    }

}
