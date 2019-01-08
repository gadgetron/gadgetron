#pragma once

#include "Gadget.h"
#include "gadgetronpython_export.h"
#include <Channel.h>
#include <ismrmrd/ismrmrd.h>
#include <boost/python.hpp>
#include "mri_core_data.h"

namespace Gadgetron{

  class EXPORTGADGETSPYTHON GadgetReference
  {

  public:
    explicit GadgetReference(Core::OutputChannel& channel);
    ~GadgetReference() = default;


    void return_acquisition(ISMRMRD::AcquisitionHeader acq, boost::python::object arr);
    void return_recondata(boost::python::object arr);
    void return_ismrmrd_image_array(boost::python::object rec);
    void return_image_cplx(ISMRMRD::ImageHeader img, boost::python::object arr);
    void return_image_cplx_attr(ISMRMRD::ImageHeader img, boost::python::object arr, const char* meta);
    void return_image_float(ISMRMRD::ImageHeader img, boost::python::object arr);
    void return_image_float_attr(ISMRMRD::ImageHeader img, boost::python::object arr, const char* meta);
    void return_image_ushort(ISMRMRD::ImageHeader img, boost::python::object arr);
    void return_image_ushort_attr(ISMRMRD::ImageHeader img, boost::python::object arr, const char* meta);

  private:

      template<class TH, class TD> void return_data(TH header, boost::python::object arr, const char* meta = 0);
   Core::OutputChannel& output;
  };
}
