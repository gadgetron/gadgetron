#ifndef FFT_H
#define FFT_H

#include <fftw3.h>

#include "NDArray.h"

template <class T> class FFT
{

public:
  FFT() { set_function_pointers(); }
  
  virtual ~FFT() { fftw_cleanup_ptr_(); }


  void fft(NDArray< std::complex<T> >* input, int dim_to_transform)
  {
    
  }
  
  
  void ifft(NDArray< std::complex<T> >* input, int dim_to_transform)
  {
    
  }
  
  void fft(NDArray< std::complex<T> >* input)
  {
    
  }
  
  void ifft(NDArray< std::complex<T> >* input)
  {
    
  }

protected:
	void fft_int(NDArray< std::complex<T> >& input, int dim_to_transform, int sign);

	void set_function_pointers();

	int   (*fftw_import_wisdom_from_file_ptr_)(FILE*);
	void  (*fftw_export_wisdom_to_file_ptr_)(FILE*);
	void  (*fftw_cleanup_ptr_)(void);
	void* (*fftw_malloc_ptr_)(size_t);
	void  (*fftw_free_ptr_)(void* p);
	void  (*fftw_execute_ptr_)(void*);
	void* (*fftw_plan_dft_1d_ptr_)(int, void*, void*, int, unsigned);
	void  (*fftw_destroy_plan_ptr_)(void*);
};

template<> void FFT<float>::set_function_pointers()
{
    fftw_import_wisdom_from_file_ptr_ = &fftwf_import_wisdom_from_file;
    fftw_export_wisdom_to_file_ptr_ = &fftwf_export_wisdom_to_file;
    fftw_cleanup_ptr_ = &fftwf_cleanup;
    fftw_malloc_ptr_ = &fftwf_malloc;
    fftw_free_ptr_ = &fftwf_free;
    fftw_execute_ptr_ = (void (*)(void*))(&fftwf_execute);
    fftw_plan_dft_1d_ptr_ = (void* (*)(int, void*, void*, int, unsigned))(&fftwf_plan_dft_1d);
    fftw_destroy_plan_ptr_ = (void (*)(void*))(&fftwf_destroy_plan);
}

template<> void FFT<double>::set_function_pointers()
{
    fftw_import_wisdom_from_file_ptr_ = &fftw_import_wisdom_from_file;
    fftw_export_wisdom_to_file_ptr_ = &fftw_export_wisdom_to_file;
    fftw_cleanup_ptr_ = &fftw_cleanup;
    fftw_malloc_ptr_ = &fftw_malloc;
    fftw_free_ptr_ = &fftw_free;
    fftw_execute_ptr_ = (void (*)(void*))(&fftw_execute);
    fftw_plan_dft_1d_ptr_ = (void* (*)(int, void*, void*, int, unsigned))(&fftw_plan_dft_1d);
    fftw_destroy_plan_ptr_ = (void (*)(void*))(&fftw_destroy_plan);
}

#endif //FFT_H
