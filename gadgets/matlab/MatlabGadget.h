#pragma once 


#include "gadgetronmatlab_export.h"
#include "Gadget.h"
#include "Gadgetron.h"
#include "hoNDArray.h"
#include "ismrmrd.h"

#include "MatlabCommunicator.h"
#include "GadgetStreamController.h"

#include <stdio.h>
#include <stdlib.h>

#include <complex>

template <class T> class MatlabGadget :
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

    GADGET_DEBUG2("MATLAB Ref Function    : %s\n", reffunc_.get()->c_str());
    GADGET_DEBUG2("MATLAB Data Function   : %s\n", datafunc_.get()->c_str());
    GADGET_DEBUG2("MATLAB Config Function : %s\n", configfunc_.get()->c_str());

    //MatlabCommunicator::instance()->register_gadget(this);
    //MatlabCommunicator::instance()->register_gadget(this->controller_->find_gadget(this->next()->module()->name()));
    MatlabCommunicator::instance();

    return GADGET_OK;
  }

protected:
  boost::shared_ptr<std::string> path_;
  boost::shared_ptr<std::string> reffunc_;
  boost::shared_ptr<std::string> datafunc_;
  boost::shared_ptr<std::string> configfunc_;


};



class EXPORTGADGETSMATLAB AcquisitionMatlabGadget :
public MatlabGadget<ISMRMRD::AcquisitionHeader>
{
 public:
  GADGET_DECLARE(AcquisitionMatlabGadget);
  
  int process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1,
  	      GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2);

};

class EXPORTGADGETSMATLAB ImageMatlabGadget :
	public MatlabGadget<ISMRMRD::ImageHeader>
{
 public:
  GADGET_DECLARE(ImageMatlabGadget);

  int process(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1,
  	      GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2);

};
