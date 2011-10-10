#pragma once 

#include "gadgetron_export.h"
#include "Gadget.h"
#include "hoNDArray.h"
#include "../core/GadgetMRIHeaders.h"
#include "GadgetReference.h"

#include <boost/python.hpp>

#include <complex>

class EXPORTGADGETSCORE AcquisitionPythonGadget : 
public Gadget2<GadgetMessageAcquisition,hoNDArray< std::complex<float> > >
{
 public:
  GADGET_DECLARE(AcquisitionPythonGadget);

  ~AcquisitionPythonGadget();
  
 protected:
  virtual int process(GadgetContainerMessage<GadgetMessageAcquisition>* m1,
		      GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2);

  virtual int process_config(ACE_Message_Block* mb);

 private:
  boost::python::object python_module_;
  boost::python::object python_set_gadget_reference_function_;
  boost::python::object python_input_function_;
  GadgetReference gadget_reference_;
};

