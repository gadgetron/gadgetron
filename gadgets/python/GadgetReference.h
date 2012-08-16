#pragma once

#include "Gadget.h"
#include "../core/GadgetMRIHeaders.h"


#include <boost/python.hpp>
#include <boost/python/numeric.hpp>
#include <boost/python/tuple.hpp>
#include "gadgetronpython_export.h"
#include "ismrmrd.h"

class EXPORTGADGETSPYTHON GadgetReference
{

 public:
  GadgetReference();
  ~GadgetReference();
  
  int set_gadget(Gadget* g)
  {
    gadget_ = g;
    return 0;
  }

  template<class T> int return_data(T header, boost::python::numeric::array arr);
  int return_acquisition(ISMRMRD::AcquisitionHeader acq, boost::python::numeric::array arr);
  int return_image(GadgetMessageImage img, boost::python::numeric::array arr);

 protected:
  Gadget* gadget_;

};
