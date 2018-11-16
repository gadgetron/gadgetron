/** \file cuNFFT.h
    \brief Cuda implementation of the non-Cartesian FFT

    Reference information on the CUDA/GPU implementation of the NFFT can be found in the papers
    
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

  /** \class cuNFFT_plan
      \brief Cuda implementation of the non-Cartesian FFT

      ------------------------------
      --- NFFT class declaration ---
      ------------------------------      
      REAL:  desired precision : float or double
      D:  dimensionality : { 1,2,3,4 }
      ATOMICS: use atomic device memory transactions : { true, false }
      
      For the tested hardware the implementation using atomic operations is slower as its non-atomic counterpart.
      However, using atomic operations has the advantage of not requiring any pre-processing.
      As the preprocessing step can be quite costly in terms of memory usage,
      the atomic mode can be necessary for very large images or for 3D/4D volumes.
      Notice: currently no devices support atomics operations in double precision.
  */
  template< class REAL, unsigned int D, bool ATOMICS = false > class EXPORTGPUNFFT cuNFFT_plan
  {
  
  public: // Main interface
    
    /** 
        Default constructor
    */
    cuNFFT_plan();

    /**
       Constructor defining the required NFFT parameters.
       \param matrix_size the matrix size to use for the NFFT. Define as a multiple of 32.
       \param matrix_size_os intermediate oversampled matrix size. Define as a multiple of 32.
       The ratio between matrix_size_os and matrix_size define the oversampling ratio for the NFFT implementation.
       Use an oversampling ratio between 1 and 2. The higher ratio the better quality results, 
       however at the cost of increased execution times. 
       \param W the concolution window size used in the NFFT implementation. 
       The larger W the better quality at the cost of increased runtime.
       \param device the device (GPU id) to use for the NFFT computation. 
       The default value of -1 indicates that the currently active device is used.
    */
    cuNFFT_plan( typename uint64d<D>::Type matrix_size, typename uint64d<D>::Type matrix_size_os,
                 REAL W, int device = -1 );

    /**
       Destructor
    */
    virtual ~cuNFFT_plan();

    /** 
        Enum to specify the desired mode for cleaning up when using the wipe() method.
    */
    enum NFFT_wipe_mode { 
      NFFT_WIPE_ALL, /**< delete all internal memory. */
      NFFT_WIPE_PREPROCESSING /**< delete internal memory holding the preprocessing data structures. */
    };

    /** 
        Clear internal storage
        \param mode enum defining the wipe mode
    */
    void wipe( NFFT_wipe_mode mode );

    /** 
        Setup the plan. Please see the constructor taking similar arguments for a parameter description.
    */
    void setup( typename uint64d<D>::Type matrix_size, typename uint64d<D>::Type matrix_size_os,
                REAL W, int device = -1 );

    /**
       Enum to specify the preprocessing mode.
    */
    enum NFFT_prep_mode { 
      NFFT_PREP_C2NC, /**< preprocess to perform a Cartesian to non-Cartesian NFFT. */
      NFFT_PREP_NC2C, /**< preprocess to perform a non-Cartesian to Cartesian NFFT. */
      NFFT_PREP_ALL /**< preprocess to perform NFFTs in both directions. */
    };

    /**
       Perform NFFT preprocessing for a given trajectory.
       \param trajectory the NFFT non-Cartesian trajectory normalized to the range [-1/2;1/2]. 
       \param mode enum specifying the preprocessing mode
    */
    void preprocess( cuNDArray<typename reald<REAL,D>::Type> *trajectory, NFFT_prep_mode mode );

    /**
       Enum defining the desired NFFT operation
    */
    enum NFFT_comp_mode { 
      NFFT_FORWARDS_C2NC, /**< forwards NFFT Cartesian to non-Cartesian. */
      NFFT_FORWARDS_NC2C, /**< forwards NFFT non-Cartesian to Cartesian. */
      NFFT_BACKWARDS_C2NC, /**< backwards NFFT Cartesian to non-Cartesian. */
      NFFT_BACKWARDS_NC2C /**< backwards NFFT non-Cartesian to Cartesian. */
    };

    /**
       Execute the NFFT.
       \param[in] in the input array.
       \param[out] out the output array.
       \param[in] dcw optional density compensation weights weighing the input samples according to the sampling density. 
       If an 0x0-pointer is provided no density compensation is used.
       \param mode enum specifying the mode of operation.
    */
    void compute( cuNDArray<complext<REAL> > *in, cuNDArray<complext<REAL> > *out,
                  cuNDArray<REAL> *dcw, NFFT_comp_mode mode );

    /**
       Execute an NFFT iteraion (from Cartesian image space to non-Cartesian Fourier space and back to Cartesian image space).
       \param[in] in the input array.
       \param[out] out the output array.
       \param[in] dcw optional density compensation weights weighing the input samples according to the sampling density. 
       If an 0x0-pointer is provided no density compensation is used.
       \param[in] halfway_dims specifies the dimensions of the intermediate Fourier space (codomain).
    */
    void mult_MH_M( cuNDArray<complext<REAL> > *in, cuNDArray<complext<REAL> > *out,
                    cuNDArray<REAL> *dcw, std::vector<size_t> halfway_dims );
  
  public: // Utilities
  
    /**
       Enum specifying the direction of the NFFT standalone convolution
    */
    enum NFFT_conv_mode { 
      NFFT_CONV_C2NC, /**< convolution: Cartesian to non-Cartesian. */
      NFFT_CONV_NC2C /**< convolution: non-Cartesian to Cartesian. */
    };
    
    /**
       Perform "standalone" convolution
       \param[in] in the input array.
       \param[out] out the output array.
       \param[in] dcw optional density compensation weights.
       \param[in] mode enum specifying the mode of the convolution
       \param[in] accumulate specifies whether the result is added to the output (accumulation) or if the output is overwritten.
    */
    void convolve( cuNDArray<complext<REAL> > *in, cuNDArray<complext<REAL> > *out, cuNDArray<REAL> *dcw,
                   NFFT_conv_mode mode, bool accumulate = false );
    
    /**
       Enum specifying the direction of the NFFT standalone FFT.
    */
    enum NFFT_fft_mode { 
      NFFT_FORWARDS, /**< forwards FFT. */
      NFFT_BACKWARDS /**< backwards FFT. */
    };

    /**
       Cartesian FFT. For completeness, just invokes the cuNDFFT class.
       \param[in,out] data the data for the inplace FFT.
       \param mode enum specifying the direction of the FFT.
       \param do_scale boolean specifying whether FFT normalization is desired.
    */
    void fft( cuNDArray<complext<REAL> > *data, NFFT_fft_mode mode, bool do_scale = true );
  
    /**
       NFFT deapodization.
       \param[in,out] image the image to be deapodized (inplace).
    */
    void deapodize( cuNDArray<complext<REAL> > *image, bool fourier_domain=false);

  public: // Setup queries
    
    /**
       Get the matrix size.
    */
    inline typename uint64d<D>::Type get_matrix_size(){
      return matrix_size;
    }

    /**
       Get the oversampled matrix size.
    */
    inline typename uint64d<D>::Type get_matrix_size_os(){
      return matrix_size_os;
    }

    /**
       Get the convolution kernel size
    */
    inline REAL get_W(){
      return W;
    }
    
    /**
       Get the assigned device id
    */
    inline unsigned int get_device(){
      return device;
    }
    
    /**
       Query of the plan has been setup
    */
    inline bool is_setup(){
      return initialized;
    }
    
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
    boost::shared_ptr<cuNDArray<complext<REAL> > > compute_deapodization_filter(bool FFTed = false);

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
    
    typename uint64d<D>::Type matrix_size;          // Matrix size
    typename uint64d<D>::Type matrix_size_os;       // Oversampled matrix size
    typename uint64d<D>::Type matrix_size_wrap;     // Wrap size at border

    typename reald<REAL,D>::Type alpha;           // Oversampling factor (for each dimension)
    typename reald<REAL,D>::Type beta;            // Kaiser-Bessel convolution kernel control parameter

    REAL W;                                       // Kernel width in oversampled grid

    unsigned int number_of_samples;               // Number of samples per frame per coil
    unsigned int number_of_frames;                // Number of frames per reconstruction
    
    int device;                                   // Associated device id

    //
    // Internal data structures for convolution and deapodization
    //

    boost::shared_ptr< cuNDArray<complext<REAL> > > deapodization_filter; //Inverse fourier transformed deapodization filter

    boost::shared_ptr< cuNDArray<complext<REAL> > > deapodization_filterFFT; //Fourier transformed deapodization filter
   
    thrust::device_vector< typename reald<REAL,D>::Type > *trajectory_positions;
    thrust::device_vector<unsigned int> *tuples_last;
    thrust::device_vector<unsigned int> *bucket_begin, *bucket_end;

    //
    // State variables
    //

    bool preprocessed_C2NC, preprocessed_NC2C;
    bool initialized;
  };

  // Pure virtual class to cause compile errors if you try to use NFFT with double and atomics
  // - since this is not supported on the device
  template< unsigned int D> class EXPORTGPUNFFT cuNFFT_plan<double,D,true>{ 
    virtual void atomics_not_supported_for_type_double() = 0; };
}
