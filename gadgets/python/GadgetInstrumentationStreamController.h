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
  int put_config(const char* config);
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


 class StreamClosingTask : public ACE_Task<ACE_MT_SYNCH>
 {
   typedef ACE_Task<ACE_MT_SYNCH> inherited;

 public:
   StreamClosingTask(GadgetInstrumentationStreamController* cntrl)
     : inherited()
     , mtx_("StreamControllerMTX")
     , cntrl_(0)
     , closed_(false)
     {
       cntrl_ = cntrl;
     }

   virtual int open(void* = 0)
   {
     return this->activate( THR_NEW_LWP | THR_JOINABLE, 1);
   }
	
   virtual int svc(void)
   {
     if (cntrl_) {
       cntrl_->close();
     }
     mtx_.acquire();
     closed_ = true;
     mtx_.release();
   }

   bool is_closed()
   {
     bool ret;
     mtx_.acquire();
     ret = closed_;
     mtx_.release();
     return ret;
   }
     

 protected:
   ACE_Thread_Mutex mtx_;
   GadgetInstrumentationStreamController* cntrl_;
   bool closed_;
 };

class GadgetInstrumentationStreamControllerWrapper
{
 public:
  GadgetInstrumentationStreamControllerWrapper() 
    : closer_(0)
    {
      // ensure boost can convert between hoNDArrays and NumPy arrays automatically
      register_converter<hoNDArray<std::complex<float> > >();
      register_converter<hoNDArray< float > >();
      // ensure boost can convert ISMRMRD headers automatically
      register_converter<ISMRMRD::ImageHeader>();
      register_converter<ISMRMRD::AcquisitionHeader>();

      cntrl_ = new GadgetInstrumentationStreamController;
    }

  ~GadgetInstrumentationStreamControllerWrapper()
    {
      delete cntrl_;
      if (closer_) delete closer_;
    }

  int prepend_gadget(const char* gadgetname, 
		    const char* dllname, 
		    const char* classname)
  {
    return cntrl_->prepend_gadget(gadgetname,dllname,classname);
  }

  int put_config(const char* config)
  {
    return cntrl_->put_config(config);
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
    //We will spawn a seperate thread task to close the stream to avoid blocking Python while closing
    if (!closer_) {
      closer_ = new StreamClosingTask(cntrl_);
      closer_->open();
    }
    return 0;
  }

  bool is_closed()
  {
    if (!closer_) {
      return false;
    }
    return closer_->is_closed();
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
  StreamClosingTask* closer_;
  

};

}
#endif //GADGETINSTRUMENTATIONSTREAMCONTROLLER_H
