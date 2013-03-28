/*
  CUDA/GPU implementation of the NFFT.

  -----------

  Accelerating the Non-equispaced Fast Fourier Transform on Commodity Graphics Hardware.
  T.S. Sørensen, T. Schaeffter, K.Ø. Noe, M.S. Hansen. 
  IEEE Transactions on Medical Imaging 2008; 27(4):538-547.

  Real-time Reconstruction of Sensitivity Encoded Radial Magnetic Resonance Imaging Using a Graphics Processing Unit.
  T.S. Sørensen, D. Atkinson, T. Schaeffter, M.S. Hansen.
  IEEE Transactions on Medical Imaging 2009; 28(12):1974-1985. 
*/

#pragma once

#include "cuNDArray.h"
#include "vector_td.h"
#include "complext.h"
#include "gpunfft_export.h"

#include <thrust/device_vector.h>
#include <boost/shared_ptr.hpp>

template<class REAL, unsigned int D, bool ATOMICS> struct _convolve_NFFT_NC2C;

namespace Gadgetron{

  //  ------------------------------
  //  --- NFFT class declaration ---
  //  ------------------------------
  //
  //     REAL:  Desired precision : float or double
  //        D:  Dimensionality : { 1,2,3,4 }
  //  ATOMICS:  NC2C convolution using atomic device memory transactions : { true, false }
  //

  template< class REAL, unsigned int D, bool ATOMICS = false > class EXPORTGPUNFFT NFFT_plan
  {
  
  public: // Main interface
    
    // Constructors
    NFFT_plan();
    NFFT_plan( typename uintd<D>::Type matrix_size, typename uintd<D>::Type matrix_size_os,
	       REAL W, int device = -1 );

    // Destructor
    virtual ~NFFT_plan();

    // Clear internal storage in plan
    enum NFFT_wipe_mode { NFFT_WIPE_ALL, NFFT_WIPE_PREPROCESSING };
    void wipe( NFFT_wipe_mode mode );

    // Replan 
    void setup( typename uintd<D>::Type matrix_size, typename uintd<D>::Type matrix_size_os,
		REAL W, int device = -1 );
  
    // Preproces trajectory ( Cartesian to non-Cartesian / non-Cartesian to Cartesian / both )
    enum NFFT_prep_mode { NFFT_PREP_C2NC, NFFT_PREP_NC2C, NFFT_PREP_ALL };
    void preprocess( cuNDArray<typename reald<REAL,D>::Type> *trajectory, NFFT_prep_mode mode );
    
    // Execute NFFT (forwards/backwards; modes: Cartesian to non-Cartesian or non-Cartesian to Cartesian)
    enum NFFT_comp_mode { NFFT_FORWARDS_C2NC, NFFT_FORWARDS_NC2C, 
			  NFFT_BACKWARDS_C2NC, NFFT_BACKWARDS_NC2C };
    void compute( cuNDArray<complext<REAL> > *in, cuNDArray<complext<REAL> > *out,
		  cuNDArray<REAL> *dcw, NFFT_comp_mode mode );

    // Execute NFFT iteration (from Cartesian image space to non-Cartesian Fourier space and back to Cartesian image space)
    void mult_MH_M( cuNDArray<complext<REAL> > *in, cuNDArray<complext<REAL> > *out,
		    cuNDArray<REAL> *dcw, std::vector<unsigned int> halfway_dims );
  
  public: // Utilities
  
    // NFFT convolution (Cartesian to non-Cartesian or non-Cartesian to Cartesian)
    enum NFFT_conv_mode { NFFT_CONV_C2NC, NFFT_CONV_NC2C };
    void convolve( cuNDArray<complext<REAL> > *in, cuNDArray<complext<REAL> > *out, cuNDArray<REAL> *dcw,
		   NFFT_conv_mode mode, bool accumulate = false );
    
    // NFFT FFT
    enum NFFT_fft_mode { NFFT_FORWARDS, NFFT_BACKWARDS };
    void fft( cuNDArray<complext<REAL> > *data, NFFT_fft_mode mode, bool do_scale = true );
  
    // NFFT deapodization
    void deapodize( cuNDArray<complext<REAL> > *image );

  public: // Setup queries

    typename uintd<D>::Type get_matrix_size();
    typename uintd<D>::Type get_matrix_size_os();
    REAL get_W();
    unsigned int get_device();
  
  public: 

    // Custom operators new/delete for windows memory handling across dll boundaries
    void* operator new (size_t bytes) { return ::new char[bytes]; }
    void operator delete (void *ptr) { delete [] static_cast <char *> (ptr); } 
    void * operator new(size_t s, void * p) { return p; }

    friend struct _convolve_NFFT_NC2C<REAL,D,ATOMICS>;
  
  private: // Internal to the implementation

    // Validate setup / arguments
    enum NFFT_components { _NFFT_CONV_C2NC = 1, _NFFT_CONV_NC2C = 2, _NFFT_FFT = 4, _NFFT_DEAPODIZATION = 8 };
    void check_consistency( cuNDArray<complext<REAL> > *samples, cuNDArray<complext<REAL> > *image,
			    cuNDArray<REAL> *dcw, unsigned char components );

    // Shared barebones constructor
    void barebones();
    
    // Compute beta control parameter for Kaiser-Bessel kernel
    void compute_beta();

    // Compute deapodization filter
    void compute_deapodization_filter();

    // Dedicated computes
    void compute_NFFT_C2NC( cuNDArray<complext<REAL> > *in, cuNDArray<complext<REAL> > *out );
    void compute_NFFT_NC2C( cuNDArray<complext<REAL> > *in, cuNDArray<complext<REAL> > *out );
    void compute_NFFTH_NC2C( cuNDArray<complext<REAL> > *in, cuNDArray<complext<REAL> > *out );
    void compute_NFFTH_C2NC( cuNDArray<complext<REAL> > *in, cuNDArray<complext<REAL> > *out );

    // Dedicated convolutions
    void convolve_NFFT_C2NC( cuNDArray<complext<REAL> > *in, cuNDArray<complext<REAL> > *out, bool accumulate );
    void convolve_NFFT_NC2C( cuNDArray<complext<REAL> > *in, cuNDArray<complext<REAL> > *out, bool accumulate );
  
    // Internal utility
    void image_wrap( cuNDArray<complext<REAL> > *in, cuNDArray<complext<REAL> > *out, bool accumulate );

  private:
    
    typename uintd<D>::Type matrix_size;          // Matrix size
    typename uintd<D>::Type matrix_size_os;       // Oversampled matrix size
    typename uintd<D>::Type matrix_size_wrap;     // Wrap size at border

    REAL alpha;                                   // Oversampling factor
    REAL beta;                                    // Kaiser-Bessel convolution kernel control parameter
    REAL W;                                       // Kernel width in oversampled grid

    unsigned int number_of_samples;               // Number of samples per frame per coil
    unsigned int number_of_frames;                // Number of frames per reconstruction
    
    int device;                                   // Associated device id

    //
    // Internal data structures for convolution and deapodization
    //

    boost::shared_ptr< cuNDArray<complext<REAL> > > deapodization_filter;
   
    thrust::device_vector< typename reald<REAL,D>::Type > *trajectory_positions;
    thrust::device_vector<unsigned int> *tuples_last;
    thrust::device_vector<unsigned int> *bucket_begin, *bucket_end;

    //
    // State variables
    //

    bool preprocessed_C2NC, preprocessed_NC2C;
    bool initialized;
  };

  //Empty class to cause compile errors if you try to use NFFT with double and atomics
  template< unsigned int D> class EXPORTGPUNFFT NFFT_plan<double,D,true>
  {};
}
