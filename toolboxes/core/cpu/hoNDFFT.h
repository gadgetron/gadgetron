/** \file hoNDFFT.h
    \brief Wrappers for FFTW for ndarrays of type std::complex.
*/

#ifndef hoNDFFT_H
#define hoNDFFT_H

#include "hoNDArray.h"
#include "cpucore_export.h"

#include <boost/thread/mutex.hpp>
#include <iostream>
#include <fftw3.h>
#include <complex>

namespace Gadgetron{

  /** 
      Generic class for Fast Fourier Transforms using FFTW on the hoNDArray class.
      This class is a singleton because the planning and memory allocation routines of FFTW are NOT threadsafe.
      The class' template type is a REAL, ie. float or double.

      Access using e.g.
      FFT<float>::instance()
  */
  template <typename T> class EXPORTCPUCORE hoNDFFT
  {

  public:
    static hoNDFFT<T>* instance(); 

    void fft(hoNDArray< std::complex<T> >* input, unsigned int dim_to_transform)
    {
      //-1 refers to the sign of the transform, -1 for FFTW_FORWARD
      fft_int(input,dim_to_transform,-1);
    }

    void ifft(hoNDArray< std::complex<T> >* input, unsigned int dim_to_transform)
    {
      //1 refers to the sign of the transform, +1 for FFTW_BACKWARD
      fft_int(input,dim_to_transform,1);
    }

    void fft(hoNDArray< std::complex<T> >* input)
    {
      for (int i = 0; i < input->get_number_of_dimensions(); i++) {
	//-1 refers to the sign of the transform, -1 for FFTW_FORWARD
	fft_int(input,i,-1);
      }
    }

    void ifft(hoNDArray< std::complex<T> >* input)
    {
      for (int i = 0; i < input->get_number_of_dimensions(); i++) {
	//1 refers to the sign of the transform, +1 for FFTW_BACKWARD
	fft_int(input,i,1);
      }
    }

  protected:

    //We are making these protected since this class is a singleton

    hoNDFFT() {
      set_function_pointers();
    }

    virtual ~hoNDFFT() { fftw_cleanup_ptr_(); }

    void fft_int(hoNDArray< std::complex<T> >* input, unsigned int dim_to_transform, int sign);

    void set_function_pointers();

    int   (*fftw_import_wisdom_from_file_ptr_)(FILE*);
    void  (*fftw_export_wisdom_to_file_ptr_)(FILE*);
    void  (*fftw_cleanup_ptr_)(void);
    void* (*fftw_malloc_ptr_)(size_t);
    void  (*fftw_free_ptr_)(void* p);
    void  (*fftw_execute_ptr_)(void*);
    void* (*fftw_plan_dft_1d_ptr_)(int, void*, void*, int, unsigned);
    void  (*fftw_destroy_plan_ptr_)(void*);

    static hoNDFFT<T>* instance_;
    boost::mutex mutex_;
  };
}

#endif //hoNDFFT_H
