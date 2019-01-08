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
        output.push(std::make_unique<IsmrmrdReconData>(boost::python::extract<IsmrmrdReconData>(rec)()));

    }

    void GadgetReference::return_ismrmrd_image_array(boost::python::object rec) {
        output.push(std::make_unique<IsmrmrdImageArray>(boost::python::extract<IsmrmrdImageArray>(rec)()));
    }

    template<class TH, class TD>
    void GadgetReference::return_data(TH header, boost::python::object arr, const char *meta) {

        auto m1 = std::make_unique<TH>(std::move(header));

        auto m2 = std::make_unique<hoNDArray<TD>>(boost::python::extract<hoNDArray<TD>>(arr));

        if (meta) {
            auto m3 = std::make_unique<ISMRMRD::MetaContainer>();
            ISMRMRD::deserialize(meta, *m3);
            output.push(std::move(m1), std::move(m2), std::move(m3));
            return;
        }

        output.push(std::move(m1), std::move(m2));
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
