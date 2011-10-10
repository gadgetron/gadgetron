#include "Gadgetron.h"
#include "Gadget.h"
#include "GadgetReference.h"
#include "GadgetContainerMessage.h"
#include "hoNDArray.h"

#include <boost/python.hpp>
#include <numpy/arrayobject.h>

#include <complex>


GadgetReference::GadgetReference()
  : gadget_(0)
{
  //_import_array();
}

GadgetReference::~GadgetReference()
{

}

template<class T>
int GadgetReference::return_data(T header, boost::python::numeric::array arr)
{

  int ndims = PyArray_NDIM(arr.ptr());
  npy_intp* dims = PyArray_DIMS(arr.ptr());
  std::vector<unsigned int> dimensions(ndims);
  for (int i = 0; i < ndims; i++) dimensions[i] = static_cast<unsigned int>(dims[i]);

  GadgetContainerMessage< T >*         m1 = new GadgetContainerMessage< T >;
  memcpy(m1->getObjectPtr(), &header, sizeof(T));

  GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2 = new GadgetContainerMessage< hoNDArray< std::complex<float> > >;
  m1->cont(m2);

  if (!m2->getObjectPtr()->create(&dimensions)) {
    GADGET_DEBUG1("Failed to create data storage for data returning from Python");
    return GADGET_FAIL;
    
  }

  memcpy(m2->getObjectPtr()->get_data_ptr(), PyArray_DATA(arr.ptr()), m2->getObjectPtr()->get_number_of_elements()*2*sizeof(float));

  if (gadget_) {
    return gadget_->next()->putq(m1);
  } else {
    GADGET_DEBUG1("Data received from python, but no Gadget registered for output\n");
    m1->release();
    return GADGET_OK;
  }

  return GADGET_OK;

}

int GadgetReference::return_acquisition(GadgetMessageAcquisition acq, boost::python::numeric::array arr)
{
  return return_data<GadgetMessageAcquisition>(acq, arr);
}

int GadgetReference::return_image(GadgetMessageImage img, boost::python::numeric::array arr)
{
  return return_data<GadgetMessageImage>(img, arr);
}

