#ifndef GADGETINSTRUMENTATIONSTREAMCONTROLLER_H
#define GADGETINSTRUMENTATIONSTREAMCONTROLLER_H

#include "GadgetStreamInterface.h"
#include "Gadget.h"
#include <ismrmrd/ismrmrd.h>
#include "python_toolbox.h"
#include <boost/python.hpp>

namespace Gadgetron{



class GadgetInstrumentationStreamController 
  : public GadgetStreamInterface
{
public:
  GadgetInstrumentationStreamController();
  int open();
  int close();
  int prepend_gadget(const char* gadgetname, 
		    const char* dllname, 
		    const char* classname);

  virtual ~GadgetInstrumentationStreamController();

  template<class T> int put_data(T header, boost::python::object arr);
  int put_acquisition(ISMRMRD::AcquisitionHeader acq, boost::python::object arr);
  int put_image(ISMRMRD::ImageHeader img, boost::python::object arr);
  int set_python_gadget(boost::python::object g)
  {
    python_gadget_ = g;
    boost::python::incref(python_gadget_.ptr());
    return GADGET_OK;
  }

  virtual int output_ready(ACE_Message_Block* mb);
  void set_parameter(const char* gadgetname, const char* parameter, const char* value);

 protected:
  boost::python::object python_gadget_;
  template <class T1, class T2> int return_data(ACE_Message_Block* mb);

};

class GadgetInstrumentationStreamControllerWrapper
{
 public:
  GadgetInstrumentationStreamControllerWrapper() 
    {
      // ensure boost can convert between hoNDArrays and NumPy arrays automatically
      register_converter<hoNDArray<std::complex<float> > >();
      // ensure boost can convert ISMRMRD headers automatically
      register_converter<ISMRMRD::ImageHeader>();
      register_converter<ISMRMRD::AcquisitionHeader>();

      cntrl_ = new GadgetInstrumentationStreamController;
    }

  ~GadgetInstrumentationStreamControllerWrapper()
    {
      delete cntrl_;
    }

  int prepend_gadget(const char* gadgetname, 
		    const char* dllname, 
		    const char* classname)
  {
    return cntrl_->prepend_gadget(gadgetname,dllname,classname);
  }

  int put_acquisition(ISMRMRD::AcquisitionHeader acq, boost::python::object arr)
  {
    return cntrl_->put_acquisition(acq, arr);
  }


  int put_image(ISMRMRD::ImageHeader img, boost::python::object arr)
  {
    return cntrl_->put_image(img,arr);
  }

  int close()
  {
    return cntrl_->close();
  }

  int set_python_gadget(boost::python::object g)
  {
    return cntrl_->set_python_gadget(g);
  }

  void set_parameter(const char* gadgetname, const char* parameter, const char* value)
  {
    cntrl_->set_parameter(gadgetname, parameter, value);
  }

 protected:
  GadgetInstrumentationStreamController* cntrl_;

};



}
#endif //GADGETINSTRUMENTATIONSTREAMCONTROLLER_H
