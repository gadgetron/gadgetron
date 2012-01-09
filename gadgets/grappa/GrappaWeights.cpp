#include "Gadgetron.h"
#include "GrappaWeights.h"
#include "hoNDArray_fileio.h"

template <class T> int GrappaWeights<T>::
update(hoNDArray< std::complex<T> >* new_weights)
{
  /*
  ACE_Guard<ACE_Thread_Mutex> guard(mutex_);
  if (!guard.locked()) {
    return -1;
  }
  */

  mutex_.acquire();

  if (!weights_.dimensions_equal(new_weights)) {
    if (!weights_.create(new_weights->get_dimensions().get())) {
      return -2;
    }
  }

  memcpy(weights_.get_data_ptr(), new_weights->get_data_ptr(),
	 weights_.get_number_of_elements()*sizeof(T)*2);

  weights_are_valid_ = true;
  mutex_.release();
  cond_.broadcast();

  return 0;
}

template<class T> int GrappaWeights<T>::
apply(hoNDArray< std::complex<T> >* data_in,
      hoNDArray< std::complex<T> >* data_out,
      T scale)
{
  /*
  ACE_Guard<ACE_Thread_Mutex> guard(mutex_);
  if (!guard.locked()) {
    return -1;
  }
  */

  mutex_.acquire();
  if (!weights_are_valid_) {
	  GADGET_DEBUG1("Releasing Mutex to Wait for result\n");
	  mutex_.release();
	  cond_.wait();
	  mutex_.acquire();
 }


  if (weights_.get_number_of_elements()%data_in->get_number_of_elements()) {
    return -3;
  }

  unsigned int sets = weights_.get_number_of_elements()/data_in->get_number_of_elements();
  
  if (sets < 1) {
    return -4;
  }

  if (data_out->get_size(data_out->get_number_of_dimensions()-1) != sets) {
    return -5;
  }

  unsigned long image_elements = data_out->get_number_of_elements()/sets;
  unsigned int coils = weights_.get_number_of_elements()/(sets*image_elements);
  
  if (weights_.get_number_of_elements() != (image_elements*coils*sets)) {
    return -6;
  }

  if (data_in->get_number_of_elements() != (image_elements*coils)) {
    return -7;
  }

  if (data_out->get_number_of_elements() != (image_elements*sets)) {
    return -8;
  }

  std::complex<T>* weights_ptr = weights_.get_data_ptr();
  std::complex<T>* in_ptr = data_in->get_data_ptr();
  std::complex<T>* out_ptr = data_out->get_data_ptr();

  for (unsigned int i = 0; i < image_elements*sets; i++) {
    out_ptr[i] = 0;
  }

  for (unsigned int s = 0; s < sets; s++) {
    for (unsigned int p = 0; p < image_elements; p++) {
      for (unsigned int c = 0; c < coils; c++) {
	out_ptr[s*image_elements + p] += 
	  weights_ptr[s*image_elements*coils + c*image_elements + p] * 
	  in_ptr[c*image_elements + p]*scale;
      }
    }
  }

  mutex_.release();
  return 0;
}

//Template instanciation
template class EXPORTGADGETSGRAPPA GrappaWeights<float>;
template class EXPORTGADGETSGRAPPA GrappaWeights<double>;
