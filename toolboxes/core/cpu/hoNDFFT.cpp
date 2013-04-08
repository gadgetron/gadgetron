/*
 * hoNDFFT.cpp
 *
 *  Created on: Nov 29, 2011
 *      Author: hansenms
 */

#include "hoNDFFT.h"

namespace Gadgetron{

  template<typename T> hoNDFFT<T>* hoNDFFT<T>::instance()
  {
    if (!instance_) instance_ = new hoNDFFT<T>();
    return instance_;
  }
  
  template<class T> hoNDFFT<T>* hoNDFFT<T>::instance_ = NULL;

  template<class T> void hoNDFFT<T>::fft_int(hoNDArray< std::complex<T> >* input, unsigned int dim_to_transform, int sign)
  {
    if (sign != -1 && sign != 1) return;
    if (dim_to_transform >= input->get_number_of_dimensions()) return;

    int stride     = 1;           //Distance between points in transform
    int dist       = 1;           //Distance between vectors
    int trafos     = 1;           //Transformations per chunk
    int chunks     = 1;           //Number of chunks
    int chunk_size = 1;           //Points per chunk
    int length     = 1;           //Length of each transform
    int total_dist = 1;

    T scale = 0.0;

    void* fft_plan        = 0;
    T*    fft_storage     = 0;

    T* fft_buffer = 0;
    T* data_ptr = 0;

    //Set sizes
    length = input->get_size(dim_to_transform);

    if (sign == 1)
      {
	scale = 1.0/length;
      }
    else
      {
	scale = 1.0;
      }

    if (dim_to_transform != 0)
      {
	for (unsigned int i = 0; i < dim_to_transform; i++)
	  {
	    chunk_size *= input->get_size(i);
	  }
	stride = chunk_size;
	trafos = chunk_size;
	chunk_size *= length;

	for (unsigned int i = dim_to_transform+1; i < input->get_number_of_dimensions(); i++)
	  {
	    chunks *= input->get_size(i);
	  }
      }
    else
      {
	for (unsigned int i = 1; i < input->get_number_of_dimensions(); i++)
	  {
	    trafos *= input->get_size(i);
	  }
	chunk_size = trafos*length;

	dist = length;
      }

    //*2 real and imag
    chunk_size *= 2;
    dist *= 2;
    total_dist = trafos*dist;


    //Allocate storage and make plan
    {
      mutex_.lock();
      fft_storage = (T*)fftw_malloc_ptr_(sizeof(T)*length*2);
      if (fft_storage == 0)
	{
	  std::cout << "Failed to allocate buffer for FFT" << std::endl;
	  return;
	}
      fft_buffer = (T*)fft_storage;

      unsigned planner_flags = FFTW_MEASURE | FFTW_DESTROY_INPUT;

      fft_plan = fftw_plan_dft_1d_ptr_(length, fft_storage, fft_storage, sign, planner_flags);

      if (fft_plan == 0)
	{
	  fftw_free_ptr_(fft_storage);
	  std::cout << "Failed to create plan for FFT" << std::endl;
	  return;
	}
      mutex_.unlock();
    }

    //Grab address of data
    data_ptr = reinterpret_cast<T*>(input->get_data_ptr());

    register int idx1_max = chunks*chunk_size;
    register int idx1, idx2;       //Index variables
    register int idx2_limit;
    register int middle_point = ((length+1)>>1)<<1;
    register int length2 = length<<1;
    register int stride2 = stride<<1;

    for (idx1 = 0; idx1 < idx1_max; idx1+=chunk_size) //Loop over all chunks
      {
	idx2_limit = idx1+total_dist;
	for (idx2 = idx1; idx2 < idx2_limit; idx2+=dist) //Loop over all transformations
	  {
	    ///Copy data to buffer.
	    {
	      register int j, idx3 = idx2;
	      for (j = middle_point; j < length2; idx3+=stride2)
		{
		  fft_buffer[j++] = data_ptr[idx3  ];
		  fft_buffer[j++] = data_ptr[idx3+1];
		}
	      for (j = 0; j < middle_point; idx3+=stride2)
		{
		  fft_buffer[j++] = data_ptr[idx3  ];
		  fft_buffer[j++] = data_ptr[idx3+1];
		}
	    }

	    fftw_execute_ptr_(fft_plan);

	    {
	      register int j, idx3 = idx2;

	      for (j = middle_point; j < length2; idx3+=stride2)
		{
		  data_ptr[idx3  ] = fft_buffer[j++]*scale;
		  data_ptr[idx3+1] = fft_buffer[j++]*scale;
		}
	      for (j = 0; j < middle_point; idx3+=stride2)
		{
		  data_ptr[idx3  ] = fft_buffer[j++]*scale;
		  data_ptr[idx3+1] = fft_buffer[j++]*scale;
		}
	    }

	  } //Loop over transformations
      } //Loop over chunks

    //clean up
    {
      mutex_.lock();
      if (fft_plan != 0)
	{
	  fftw_destroy_plan_ptr_(fft_plan);
	}

      if (fft_storage != 0)
	{
	  fftw_free_ptr_(fft_storage);
	}
      mutex_.unlock();
    }
  }
  
  template<> void hoNDFFT<float>::set_function_pointers()
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

  template<> void hoNDFFT<double>::set_function_pointers()
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

  // 
  // Instantiation
  //
  
  template class hoNDFFT<float>;
  template class hoNDFFT<double>;
}
