#pragma once 

#include <octave/oct.h>
#include <octave/octave.h>
#include <octave/parse.h>
#include <octave/ov-struct.h>

#include "gadgetronoctave_export.h"
#include "Gadget.h"
#include "Gadgetron.h"
#include "hoNDArray.h"
#include "ismrmrd.h"

#include "OctaveCommunicator.h"
#include "GadgetStreamController.h"

#include <stdio.h>
#include <stdlib.h>

#include <complex>

template <class T> class OctaveGadget :
public Gadget2<T, hoNDArray< std::complex<float> > >
{
 public:
  //GADGET_DECLARE(OctaveGadget);
  //virtual ~OctaveGadget();

 protected:

  int process_config(ACE_Message_Block* mb)
  {

    path_        = this->get_string_value("path");
    reffunc_     = this->get_string_value("gadget_reference_function");
    datafunc_    = this->get_string_value("input_function");
    configfunc_  = this->get_string_value("config_function");

    GADGET_DEBUG2("OCTAVE Ref Function    : %s\n", reffunc_.get()->c_str());
    GADGET_DEBUG2("OCTAVE Data Function   : %s\n", datafunc_.get()->c_str());
    GADGET_DEBUG2("OCTAVE Config Function : %s\n", configfunc_.get()->c_str());

    OctaveCommunicator::instance()->register_gadget(this);
    OctaveCommunicator::instance()->register_gadget(this->controller_->find_gadget(this->next()->module()->name()));

    octave_value_list in = octave_value (path_->c_str());
    octave_value_list out = feval ("addpath", in, 1);

    const char* gadgetron_home = ACE_OS::getenv("GADGETRON_HOME");
    std::string path_name = std::string(gadgetron_home) + std::string("/octave");

    in = octave_value (path_name.c_str());
    out = feval ("addpath", in, 1);

    in(0) = octave_value(this->module()->name());
    in(1) = octave_value(this->next()->module()->name());
    out = feval(reffunc_->c_str(), in, 2);

    in(0) = octave_value(std::string(mb->rd_ptr(),mb->length()));
    out = feval(configfunc_->c_str(), in, 1);

    return GADGET_OK;
  }

protected:
  boost::shared_ptr<std::string> path_;
  boost::shared_ptr<std::string> reffunc_;
  boost::shared_ptr<std::string> datafunc_;
  boost::shared_ptr<std::string> configfunc_;


};



class EXPORTGADGETSOCTAVE AcquisitionOctaveGadget :
public OctaveGadget<ISMRMRD::AcquisitionHeader>
{
 public:
  GADGET_DECLARE(AcquisitionPythonGadget);
  
  int process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1,
  	      GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2);

};

class EXPORTGADGETSOCTAVE ImageOctaveGadget :
public OctaveGadget<ISMRMRD::ImageHeader>
{
 public:
  GADGET_DECLARE(ImagePythonGadget);

  int process(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1,
  	      GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2);

};
