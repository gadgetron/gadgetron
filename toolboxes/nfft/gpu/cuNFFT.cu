/*
  CUDA implementation of the NFFT.

  -----------

  Accelerating the Non-equispaced Fast Fourier Transform on Commodity Graphics Hardware.
  T.S. Sørensen, T. Schaeffter, K.Ø. Noe, M.S. Hansen. 
  IEEE Transactions on Medical Imaging 2008; 27(4):538-547.

  Real-time Reconstruction of Sensitivity Encoded Radial Magnetic Resonance Imaging Using a Graphics Processing Unit.
  T.S. Sørensen, D. Atkinson, T. Schaeffter, M.S. Hansen.
  IEEE Transactions on Medical Imaging 2009; 28(12): 1974-1985. 
*/

// Includes - Thrust
#include <thrust/scan.h>
#include <thrust/sort.h>
#include <thrust/binary_search.h>
#include <thrust/extrema.h>
// Includes - Gadgetron
#include "cuNFFT.h"
#include "cuNDFFT.h"
#include "cuNDArray_operators.h"
#include "cuNDArray_elemwise.h"
#include "cuNDArray_blas.h"
#include "cuNDArray_utils.h"
#include "vector_td_utilities.h"
#include "vector_td_io.h"
#include "cudaDeviceManager.h"
#include "check_CUDA.h"

// Includes - CUDA
#include <cuda_runtime.h>
#include <math_constants.h>
#include <cufft.h>


// Includes - stdlibs
#include <stdio.h>
#include <assert.h>
#include <limits.h>
#include <math.h>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <numeric>
#include <driver_types.h>


#include "NFFT.hpp"
//using namespace std;
using std::vector;
using namespace thrust;
using namespace Gadgetron;

// Kernel configuration  
#define NFFT_MAX_COILS_COMPUTE_1x    8
#define NFFT_MAX_COILS_COMPUTE_2x   16
#define NFFT_THREADS_PER_KERNEL    192

// Reference to shared memory
extern __shared__ char _shared_mem[];

// Includes containing the NFFT convolution implementation
#include "KaiserBessel_kernel.h"
#include "NFFT_C2NC_conv_kernel.cu"
#include "NFFT_NC2C_conv_kernel.cu"
#include "NFFT_NC2C_atomic_conv_kernel.cu"
#include "NFFT_preprocess_kernel.cu"
#include "NFFT_sparseMatrix_kernel.cu"





//
// Public class methods
//


template<class REAL, unsigned int D, ConvolutionType CONV>
Gadgetron::cuNFFT_impl<REAL, D, CONV>::cuNFFT_impl(const vector_td<size_t,D>& matrix_size,
                                                   const  vector_td<size_t,D>& matrix_size_os, REAL W, int _device) :
                                                   NFFT_plan<cuNDArray,REAL,D>(matrix_size,matrix_size_os,W){
    // Minimal initialization
    barebones();
    //
    // Check if the this->device is valid
    //

    if (_device < 0) {
        if (cudaGetDevice(&this->device) != cudaSuccess) {
            throw cuda_error("Error: cuNFFT_impl::setup: unable to determine this->device properties.");
        }
    } else
        this->device = _device;

    // The convolution does not work properly for very small convolution kernel widths
    // (experimentally observed limit)

    if (W < REAL(1.8)) {
        throw std::runtime_error("Error: the convolution kernel width for the cuNFFT plan is too small.");
    }

    typename uint64d<D>::Type vec_warp_size((size_t) (cudaDeviceManager::Instance()->warp_size(this->device)));

    //
    // Check input against certain requirements
    //

    if (sum(matrix_size % vec_warp_size) || sum(matrix_size_os % vec_warp_size)) {
        //GDEBUG_STREAM("Matrix size: " << matrix_size << std::endl);
        //GDEBUG_STREAM("Matrix size os: " << this->matrix_size_os << std::endl);
        //GDEBUG_STREAM("Warp size: " << vec_warp_size << std::endl);
        throw std::runtime_error("Error: Illegal matrix size for the cuNFFT plan (not a multiple of the warp size)");
    }

    //
    // Setup private variables
    //


    REAL W_half = REAL(0.5) * W;
    vector_td<REAL, D>
    W_vec(W_half);

    this->matrix_size_wrap = vector_td<size_t, D>(ceil(W_vec));
    this->matrix_size_wrap <<= 1;

    alpha = vector_td<REAL, D>(matrix_size_os) / vector_td<REAL, D>(matrix_size);

    typename reald<REAL, D>::Type ones(REAL(1));
    if (weak_less(alpha, ones)) {
        throw std::runtime_error("Error: cuNFFT : Illegal oversampling ratio suggested");
    }


    // Compute Kaiser-Bessel beta
    compute_beta();

    int device_no_old;
    if (cudaGetDevice(&device_no_old) != cudaSuccess) {
        throw cuda_error("Error: cuNFFT_impl::setup: unable to get this->device no");
    }
    if (this->device != device_no_old && cudaSetDevice(this->device) != cudaSuccess) {
        throw cuda_error("Error: cuNFFT_impl::setup: unable to set this->device");
    }

    if (this->device != device_no_old && cudaSetDevice(device_no_old) != cudaSuccess) {
        throw cuda_error("Error: cuNFFT_impl::setup: unable to restore this->device");
    }
    // Setup plan
}

template<class REAL, unsigned int D, ConvolutionType CONV>
Gadgetron::cuNFFT_impl<REAL, D, CONV>::~cuNFFT_impl() {
}


template<class REAL, unsigned int D, ConvolutionType CONV>
void Gadgetron::cuNFFT_impl<REAL, D, CONV>::preprocess(const cuNDArray<typename reald<REAL, D>::Type>& trajectory,
                                                       NFFT_prep_mode mode) {

    NFFT_plan<cuNDArray,REAL,D>::preprocess(trajectory,mode);
    cuNDArray<REAL> trajectory_view (std::vector<size_t>{trajectory.get_number_of_elements()*D},(REAL*)trajectory.get_data_ptr());

    // Make sure that the trajectory values are within range [-1/2;1/2]
    thrust::pair<thrust::device_ptr<REAL>, thrust::device_ptr<REAL> > mm_pair =
            thrust::minmax_element(trajectory_view.begin(), trajectory_view.end());

    if (*mm_pair.first < REAL(-0.5) || *mm_pair.second > REAL(0.5)) {
        std::stringstream ss;
        ss << "Error: cuNFFT::preprocess : trajectory [" << *mm_pair.first << "; " << *mm_pair.second
           << "] out of range [-1/2;1/2]";
        throw std::runtime_error(ss.str());
    }

    trajectory_positions = device_vector<vector_td<REAL, D> >(trajectory.get_number_of_elements());

    CHECK_FOR_CUDA_ERROR();

    vector_td<REAL, D>
    matrix_size_os_real = vector_td<REAL, D>(this->matrix_size_os);
    vector_td<REAL, D>
    matrix_size_os_plus_wrap_real = vector_td<REAL, D>((this->matrix_size_os + this->matrix_size_wrap) >> 1);

    // convert input trajectory in [-1/2;1/2] to [0;this->matrix_size_os]
    thrust::transform(trajectory.begin(), trajectory.end(), trajectory_positions.begin(),
                      trajectory_scale<REAL, D>(matrix_size_os_real, matrix_size_os_plus_wrap_real));

    CHECK_FOR_CUDA_ERROR();

    convNC2C.prepare(this, trajectory_positions);


    preprocessed_C2NC = true;

    if (mode != NFFT_prep_mode::C2NC)
        preprocessed_NC2C = true;

}

template<class REAL, unsigned int D>
void
cuNFFT::convolverNC2C<REAL, D, ConvolutionType::STANDARD>::prepare(cuNFFT_impl<REAL, D, ConvolutionType::STANDARD> *plan,
                                                               const thrust::device_vector<vector_td<REAL, D>> &trajectory_int) {// allocate storage for and compute temporary prefix-sum variable (#cells influenced per sample)

    unsigned int number_of_frames = plan->number_of_frames;
    unsigned int number_of_samples = plan->number_of_samples;
    auto matrix_size_os = plan->matrix_size_os;
    auto matrix_size_wrap = plan->matrix_size_wrap;
    auto W = plan->W;


    device_vector<unsigned int> c_p_s(trajectory_int.size());
    device_vector<unsigned int> c_p_s_ps(trajectory_int.size());
    CHECK_FOR_CUDA_ERROR();

    REAL half_W = REAL(0.5) * W;
    plus<unsigned int> binary_op;
    transform(trajectory_int.begin(), trajectory_int.end(), c_p_s.begin(),
              compute_num_cells_per_sample<REAL, D>(half_W));
    inclusive_scan(c_p_s.begin(), c_p_s.end(), c_p_s_ps.begin(), binary_op); // prefix sum

    // Build the vector of (grid_idx, sample_idx) tuples. Actually kept in two seperate vectors.
    unsigned int num_pairs = c_p_s_ps.back();
    c_p_s.clear();

    device_vector<unsigned int> tuples_first = device_vector<unsigned int>(num_pairs);
    tuples_last = device_vector<unsigned int>(num_pairs);

    CHECK_FOR_CUDA_ERROR();

    // Fill tuple vector
    write_pairs<REAL, D>(vector_td<unsigned int, D>(matrix_size_os), vector_td<unsigned int, D>(matrix_size_wrap),
                         number_of_samples, number_of_frames, W,
                         raw_pointer_cast(&trajectory_int[0]), raw_pointer_cast(&c_p_s_ps[0]),
                         raw_pointer_cast(&tuples_first[0]), raw_pointer_cast(&tuples_last[0]));
    c_p_s_ps.clear();

    // Sort by grid indices
    sort_by_key(tuples_first.begin(), tuples_first.end(), tuples_last.begin());

    // each bucket_begin[i] indexes the first element of bucket i's list of points
// each bucket_end[i] indexes one past the last element of bucket i's list of points

    bucket_begin = device_vector<unsigned int>(number_of_frames * prod(matrix_size_os + matrix_size_wrap));
    bucket_end = device_vector<unsigned int>(number_of_frames * prod(matrix_size_os + matrix_size_wrap));

    CHECK_FOR_CUDA_ERROR();

    // find the beginning of each bucket's list of points
    counting_iterator<unsigned int> search_begin(0);
    lower_bound(tuples_first.begin(), tuples_first.end(), search_begin, search_begin + number_of_frames * prod(
            matrix_size_os + matrix_size_wrap), bucket_begin.begin());

    // find the end of each bucket's list of points
    upper_bound(tuples_first.begin(), tuples_first.end(), search_begin, search_begin + number_of_frames * prod(
            matrix_size_os + matrix_size_wrap), bucket_end.begin());

}



template<class REAL, unsigned int D, ConvolutionType CONV>
void
Gadgetron::cuNFFT_impl<REAL, D, CONV>::convolve(const cuNDArray <complext<REAL>>& in, cuNDArray <complext<REAL>> & out,
                                                 NFFT_conv_mode mode, bool accumulate) {

    {
        const cuNDArray <complext<REAL>> *samples, *image;

        if (mode == NFFT_conv_mode::C2NC) {
            image = &in;
            samples = &out;
        } else {
            image = &out;
            samples = &in;
        }

        check_consistency(samples, image, 0x0 );
    }

    typename uint64d<D>::Type image_dims = from_std_vector<size_t, D>
            (*(((mode == NFFT_conv_mode::C2NC) ? in : out).get_dimensions()));
    bool oversampled_image = (image_dims == this->matrix_size_os);

    if (!oversampled_image) {
        throw std::runtime_error("Error: cuNFFT_impl::convolve: ERROR: oversampled image not provided as input.");
    }

    switch (mode) {
        case NFFT_conv_mode::C2NC:
            convC2NC.convolve_C2NC(this,&in, &out, accumulate);
            break;
        case NFFT_conv_mode::NC2C:
                convNC2C.convolve_NC2C(this, &in, &out, accumulate);
            break;
    }

}

template<class REAL, unsigned int D, ConvolutionType CONV>
void
Gadgetron::cuNFFT_impl<REAL, D, CONV>::fft(cuNDArray <complext<REAL>>& data, NFFT_fft_mode mode, bool do_scale) {

    typename uint64d<D>::Type _dims_to_transform = counting_vec<size_t, D>();
    vector<size_t> dims_to_transform = to_std_vector(_dims_to_transform);

    if (mode == NFFT_fft_mode::FORWARDS) {
        cuNDFFT<REAL>::instance()->fft(&data, &dims_to_transform, do_scale);
    } else {
        cuNDFFT<REAL>::instance()->ifft(&data, &dims_to_transform, do_scale);
    }

}

template<class REAL, unsigned int D, ConvolutionType CONV>
void
Gadgetron::cuNFFT_impl<REAL, D, CONV>::deapodize(cuNDArray <complext<REAL>>& image, bool fourier_domain) {

    typename uint64d<D>::Type image_dims = from_std_vector<size_t, D>(*image.get_dimensions());
    bool oversampled_image = (image_dims == this->matrix_size_os);

    if (!oversampled_image) {
        throw std::runtime_error("Error: cuNFFT_impl::deapodize: ERROR: oversampled image not provided as input.");
    }
    if (fourier_domain) {
        if (!deapodization_filterFFT)
            deapodization_filterFFT = compute_deapodization_filter(true);
        image *= *deapodization_filterFFT;
    } else {
        if (!deapodization_filter)
            deapodization_filter = compute_deapodization_filter(false);
        image *= *deapodization_filter;
    }
}

//
// Private class methods
//

template<class REAL, unsigned int D, ConvolutionType CONV>
void
Gadgetron::cuNFFT_impl<REAL, D, CONV>::check_consistency(const cuNDArray <complext<REAL>> *samples,
                                                         const cuNDArray <complext<REAL>> *image,
                                                         const cuNDArray <REAL> *weights) {

//
//  if( (components & _NFFT_conv_mode::C2NC ) && !preprocessed_C2NC ){
//    throw std::runtime_error("Error: cuNFFT_impl: Unable to compute NFFT before preprocessing.");
//  }
//
//  if( (components & _NFFT_conv_mode::NC2C ) && !(preprocessed_NC2C || (preprocessed_C2NC && CONV == ConvolutionType::ATOMIC) ) ){
//    throw std::runtime_error("Error: cuNFFT_impl: Unable to compute NFFT before preprocessing.");
//  }
//
//  if( ((components & _NFFT_conv_mode::C2NC ) || (components & _NFFT_conv_mode::NC2C )) && !(image && samples) ){
//    throw std::runtime_error("Error: cuNFFT_impl: Unable to process 0x0 input/output.");
//  }
//
//  if( ((components & _NFFT_FFT) || (components & _NFFT_DEAPODIZATION )) && !image ){
//    throw std::runtime_error("Error: cuNFFT_impl: Unable to process 0x0 input.");
//  }

    if (image->get_number_of_dimensions() < D) {
        throw std::runtime_error("Error: cuNFFT_impl: Number of image dimensions mismatch the plan.");
    }

    typename uint64d<D>::Type image_dims = from_std_vector<size_t, D>(*image->get_dimensions());
    bool oversampled_image = (image_dims == this->matrix_size_os);

    if (!((oversampled_image) ? (image_dims == this->matrix_size_os) : (image_dims == this->matrix_size))) {
        throw std::runtime_error("Error: cuNFFT_impl: Image dimensions mismatch.");
    }

    if ((samples->get_number_of_elements() == 0) ||
        (samples->get_number_of_elements() % (this->number_of_frames * this->number_of_samples))) {
        printf("\ncuNFFT::check_consistency() failed:\n#elements in the samples array: %ld.\n#samples from preprocessing: %d.\n#frames from preprocessing: %d.\n",
               samples->get_number_of_elements(), this->number_of_samples, this->number_of_frames);
        fflush(stdout);
        throw std::runtime_error(
                "Error: cuNFFT_impl: The number of samples is not a multiple of #samples/frame x #frames as requested through preprocessing");
    }

    unsigned int num_batches_in_samples_array =
            samples->get_number_of_elements() / (this->number_of_frames * this->number_of_samples);
    unsigned int num_batches_in_image_array = 1;

    for (unsigned int d = D; d < image->get_number_of_dimensions(); d++) {
        num_batches_in_image_array *= image->get_size(d);
    }
    num_batches_in_image_array /= this->number_of_frames;

    if (num_batches_in_samples_array != num_batches_in_image_array) {
        printf("\ncuNFFT::check_consistency() failed:\n#elements in the samples array: %ld.\n#samples from preprocessing: %d.\n#frames from preprocessing: %d.\nLeading to %d batches in the samples array.\nThe number of batches in the image array is %d.\n",
               samples->get_number_of_elements(), this->number_of_samples, this->number_of_frames, num_batches_in_samples_array,
               num_batches_in_image_array);
        fflush(stdout);
        throw std::runtime_error("Error: cuNFFT_impl: Number of batches mismatch between samples and image arrays");
    }


//  if( components & _NFFT_conv_mode::NC2C ){
//    if( weights ){
//      if( weights->get_number_of_elements() == 0 ||
//          !( weights->get_number_of_elements() == number_of_samples ||
//             weights->get_number_of_elements() == number_of_frames*number_of_samples) ){
//        printf("\ncuNFFT::check_consistency() failed:\n#elements in the samples array: %ld.\n#samples from preprocessing: %d.\n#frames from preprocessing: %d.\n#weights: %ld.\n",samples->get_number_of_elements(), number_of_samples, number_of_frames, weights->get_number_of_elements() ); fflush(stdout);
//        throw std::runtime_error("Error: cuNFFT_impl: The number of weights should match #samples/frame x #frames as requested through preprocessing");
//      }
//    }
//  }
}

template<class REAL, unsigned int D, ConvolutionType CONV>
void Gadgetron::cuNFFT_impl<REAL, D, CONV>::barebones() {
    // These are the fundamental booleans checked before accessing the various member pointers
    preprocessed_C2NC = preprocessed_NC2C = false;
    // and specify the this->device
    if (cudaGetDevice(&this->device) != cudaSuccess) {
        throw cuda_error("Error: cuNFFT_impl::barebones:: unable to get this->device no");
    }
}



template<class REAL, unsigned int D, ConvolutionType CONV>
void Gadgetron::cuNFFT_impl<REAL, D, CONV>::compute_beta() {
    // Compute Kaiser-Bessel beta paramter according to the formula provided in
    // Beatty et. al. IEEE TMI 2005;24(6):799-808.
    for (unsigned int d = 0; d < D; d++)
        beta[d] = (M_PI * std::sqrt(
                (this->W * this->W) / (alpha[d] * alpha[d]) * (alpha[d] - REAL(0.5)) * (alpha[d] - REAL(0.5)) -
                REAL(0.8)));
}

//
// Grid fictitious trajectory with a single sample at the origin
//

template<class REAL, unsigned int D>
__global__ void
compute_deapodization_filter_kernel(typename uintd<D>::Type matrix_size_os,
                                    typename reald<REAL, D>::Type matrix_size_os_real,
                                    REAL W, REAL half_W, REAL one_over_W,
                                    typename reald<REAL, D>::Type beta, complext <REAL> *__restrict__ image_os) {
    const unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
    const unsigned int num_elements = prod(matrix_size_os);

    if (idx < num_elements) {

        // Compute weight from Kaiser-Bessel filter
        const typename uintd<D>::Type cell_pos = idx_to_co<D>(idx, matrix_size_os);

        // Sample position ("origin")
        const vector_td<REAL, D> sample_pos = REAL(0.5) * matrix_size_os_real;

        // Calculate the distance between the cell and the sample
        vector_td<REAL, D>
        cell_pos_real = vector_td<REAL, D>(cell_pos);
        const typename reald<REAL, D>::Type delta = abs(sample_pos - cell_pos_real);

        // Compute convolution weight.
        REAL weight;
        vector_td<REAL, D>
        half_W_vec(half_W);

        if (weak_greater(delta, half_W_vec))
            weight = 0.0f;
        else {
            weight = KaiserBessel < REAL > (delta, matrix_size_os_real, one_over_W, beta);
            //if( !isfinite(weight) )
            //weight = zero;
        }

        // Output weight
        image_os[idx] = complext<REAL>(weight, 0.0f);
    }
}

//
// Function to calculate the deapodization filter
//

template<class REAL, unsigned int D, ConvolutionType CONV> boost::shared_ptr<cuNDArray < complext < REAL> > >

Gadgetron::cuNFFT_impl<REAL, D, CONV>::compute_deapodization_filter(bool FFTed) {
    std::vector<size_t> tmp_vec_os = to_std_vector(this->matrix_size_os);

    auto  filter = boost::make_shared<cuNDArray<complext<REAL>>>(tmp_vec_os);
    vector_td<REAL, D>
    matrix_size_os_real = vector_td<REAL, D>(this->matrix_size_os);

    // Find dimensions of grid/blocks.
    dim3 dimBlock(256);
    dim3 dimGrid((prod(this->matrix_size_os) + dimBlock.x - 1) / dimBlock.x);

    // Invoke kernel
    compute_deapodization_filter_kernel<REAL, D> << < dimGrid, dimBlock >> >
                                                               (vector_td<unsigned int, D>(
                                                                       this->matrix_size_os), matrix_size_os_real, this->W,
                                                                       REAL(0.5) * this->W, REAL(1) /
                                                                                            this->W, beta, filter->get_data_ptr());

    CHECK_FOR_CUDA_ERROR();

    // FFT
    if (FFTed) {
        fft(*filter, NFFT_fft_mode::FORWARDS, false);
    } else {
        fft(*filter, NFFT_fft_mode::BACKWARDS, false);
    }
    // Reciprocal
    reciprocal_inplace(filter.get());
    return filter;
}





template<class REAL, unsigned int D>
    void cuNFFT::convolverNC2C<REAL, D, ConvolutionType::SPARSE_MATRIX>::prepare(cuNFFT_impl<REAL, D, ConvolutionType::SPARSE_MATRIX> *plan,
                                         const thrust::device_vector<vector_td<REAL, D>> &trajectory) {
        matrix = std::make_unique<cuCsrMatrix<complext<REAL>>>(::make_NFFT_matrix(trajectory,plan->get_matrix_size_os(),plan->beta,plan->W));
}

template<class REAL, unsigned int D>
void cuNFFT::convolverNC2C<REAL, D, ConvolutionType::SPARSE_MATRIX>::convolve_NC2C(
        Gadgetron::cuNFFT_impl<REAL, D, ConvolutionType::SPARSE_MATRIX> *plan,
        const Gadgetron::cuNDArray<complext < REAL>> *samples, Gadgetron::cuNDArray<complext < REAL>> *image, bool accumulate) {

unsigned int number_of_frames = plan->number_of_frames;
unsigned int number_of_samples = plan->number_of_samples;
typename uint64d<D>::Type matrix_size_os = plan->matrix_size_os;
unsigned int num_batches = 1;

for (unsigned int d = D; d<image->get_number_of_dimensions(); d++){
    num_batches *= image->get_size(d);
}
num_batches /= number_of_frames;

std::vector<size_t> sample_dims = {number_of_samples*number_of_frames, num_batches};
//std::vector<size_t> image_dims = to_std_vector(matrix_size_os);
std::vector<size_t> image_dims = {prod(matrix_size_os),num_batches};
image_dims.
push_back(number_of_frames);

//size_t image_elements = std::accumulate(image_dims.begin(), image_dims.end(), 1, std::multiplies<size_t>());
//size_t sample_elements = std::accumulate(sample_dims.begin(), sample_dims.end(), 1, std::multiplies<size_t>());

//for (int batch = 0; batch<num_batches; batch++) {
//
cuNDArray <complext<REAL>> image_view(image_dims, image->get_data_ptr());
cuNDArray <complext<REAL>> samples_view(sample_dims, const_cast<complext<REAL>*>(samples->get_data_ptr()));
//
//sparseMV(complext<REAL>(1.0), complext<REAL>(1.0), transposed, samples_view, image_view,false);
//}
sparseMM(complext<REAL>(1.0),complext<REAL>(1.0),*matrix,samples_view,image_view,true);


}


template<unsigned int D>
void cuNFFT::convolverNC2C<float, D, ConvolutionType::ATOMIC>::convolve_NC2C(
        cuNFFT_impl<float, D, ConvolutionType::ATOMIC> *plan,
        const cuNDArray <complext<float>> *samples,
        cuNDArray <complext<float>> *image,
        bool accumulate) {
    //
    // Bring in some variables from the plan

    unsigned int device = plan->device;
    unsigned int number_of_frames = plan->number_of_frames;
    unsigned int number_of_samples = plan->number_of_samples;
    typename uint64d<D>::Type matrix_size_os = plan->matrix_size_os;
    typename uint64d<D>::Type matrix_size_wrap = plan->matrix_size_wrap;
    typename reald<float, D>::Type alpha = plan->alpha;
    typename reald<float, D>::Type beta = plan->beta;
    float W = plan->W;
    auto & trajectory_positions = plan->trajectory_positions;

    //
    // Atomic operations are only supported in compute model 2.0 and up
    //

    if (cudaDeviceManager::Instance()->major_version(device) == 1) {
        throw cuda_error("Error: Atomic NC2C NFFT only supported on this->device with compute model 2.0 or higher");
    }

    // Check if warp_size is a power of two. We do some modulus tricks in the kernels that depend on this...
    if (!((cudaDeviceManager::Instance()->warp_size(device) & (cudaDeviceManager::Instance()->warp_size(device) - 1)) ==
          0)) {
        throw cuda_error("cuNFFT: unsupported hardware (warpSize is not a power of two)");
    }

    unsigned int num_batches = 1;
    for (unsigned int d = D; d < image->get_number_of_dimensions(); d++)
        num_batches *= image->get_size(d);
    num_batches /= number_of_frames;

    //
    //  Setup grid and threads
    //

    size_t threads_per_block;
    unsigned int max_coils;

    threads_per_block = NFFT_THREADS_PER_KERNEL;
    max_coils = NFFT_MAX_COILS_COMPUTE_2x;

    // We can (only) convolve domain_size_coils batches per run due to shared memory issues. 
    unsigned int domain_size_coils_desired = num_batches;
    unsigned int num_repetitions = domain_size_coils_desired / max_coils +
                                   (((domain_size_coils_desired % max_coils) == 0) ? 0 : 1);
    unsigned int domain_size_coils = (num_repetitions == 1) ? domain_size_coils_desired : max_coils;
    unsigned int domain_size_coils_tail = (num_repetitions == 1) ? domain_size_coils_desired :
                                          domain_size_coils_desired - (num_repetitions - 1) * domain_size_coils;

    // Block and Grid dimensions
    dim3 dimBlock((unsigned int) threads_per_block);
    dim3 dimGrid((number_of_samples + dimBlock.x - 1) / dimBlock.x, number_of_frames);

    // Calculate how much shared memory to use per thread
    size_t bytes_per_thread = domain_size_coils * sizeof(vector_td<float, D>);
    size_t bytes_per_thread_tail = domain_size_coils_tail * sizeof(vector_td<float, D>);

    unsigned int double_warp_size_power = 0, __tmp = cudaDeviceManager::Instance()->warp_size(device) << 1;
    while (__tmp != 1) {
        __tmp >>= 1;
        double_warp_size_power++;
    }

    vector_td<float, D>
    matrix_size_os_real = vector_td<float, D>(matrix_size_os);

    if (!accumulate) {
        clear(image);
    }

    //
    // Invoke kernel
    //

    for (unsigned int repetition = 0; repetition < num_repetitions; repetition++) {

        NFFT_H_atomic_convolve_kernel<float, D>
                << < dimGrid, dimBlock,
                ((repetition == num_repetitions - 1) ? dimBlock.x * bytes_per_thread_tail : dimBlock.x *
                                                                                            bytes_per_thread) >> >
                (alpha, beta, W, vector_td<unsigned int, D>(matrix_size_os), vector_td<unsigned int, D>(
                        matrix_size_wrap), number_of_samples,
                        (repetition == num_repetitions - 1) ? domain_size_coils_tail : domain_size_coils,
                        raw_pointer_cast(&trajectory_positions[0]),
                        samples->get_data_ptr() + repetition * number_of_samples * number_of_frames * domain_size_coils,
                        image->get_data_ptr() +
                        repetition * prod(matrix_size_os) * number_of_frames * domain_size_coils,
                        double_warp_size_power, float(0.5) * W, float(1) / (W), matrix_size_os_real);
    }

    CHECK_FOR_CUDA_ERROR();
}


template<class REAL, unsigned int D>
void cuNFFT::convolverNC2C<REAL, D, ConvolutionType::STANDARD>::convolve_NC2C(cuNFFT_impl<REAL, D> *plan,
                                                                          const cuNDArray <complext<REAL>> *samples,
                                                                          cuNDArray <complext<REAL>> *image,
                                                                          bool accumulate) {
    // Bring in some variables from the plan

    unsigned int device = plan->device;
    unsigned int number_of_frames = plan->number_of_frames;
    unsigned int number_of_samples = plan->number_of_samples;
    typename uint64d<D>::Type matrix_size_os = plan->matrix_size_os;
    typename uint64d<D>::Type matrix_size_wrap = plan->matrix_size_wrap;
    typename reald<REAL, D>::Type alpha = plan->alpha;
    typename reald<REAL, D>::Type beta = plan->beta;
    REAL W = plan->W;
    auto trajectory_positions = plan->trajectory_positions;


    // private method - no consistency check. We trust in ourselves.
    // Check if warp_size is a power of two. We do some modulus tricks in the kernels that depend on this...
    if (!((cudaDeviceManager::Instance()->warp_size(device) & (cudaDeviceManager::Instance()->warp_size(device) - 1)) ==
          0)) {
        throw cuda_error("cuNFFT: unsupported hardware (warpSize is not a power of two)");

    }
    unsigned int num_batches = 1;
    for (unsigned int d = D; d < image->get_number_of_dimensions(); d++)
        num_batches *= image->get_size(d);
    num_batches /= number_of_frames;

    //
    // Setup grid and threads
    //

    size_t threads_per_block;
    unsigned int max_coils;

    threads_per_block = NFFT_THREADS_PER_KERNEL;

    if (cudaDeviceManager::Instance()->major_version(device) == 1) {
        max_coils = NFFT_MAX_COILS_COMPUTE_1x;
    } else {
        max_coils = NFFT_MAX_COILS_COMPUTE_2x;
    }

    // We can (only) convolve domain_size_coils batches per run due to shared memory issues.
    unsigned int domain_size_coils_desired = num_batches;
    unsigned int num_repetitions = domain_size_coils_desired / max_coils +
                                   (((domain_size_coils_desired % max_coils) == 0) ? 0 : 1);
    unsigned int domain_size_coils = (num_repetitions == 1) ? domain_size_coils_desired : max_coils;
    unsigned int domain_size_coils_tail = (num_repetitions == 1) ? domain_size_coils_desired :
                                          domain_size_coils_desired - (num_repetitions - 1) * domain_size_coils;

    // Block and Grid dimensions
    dim3 dimBlock((unsigned int) threads_per_block);
    dim3 dimGrid((prod(matrix_size_os + matrix_size_wrap) + dimBlock.x - 1) / dimBlock.x, number_of_frames);

    // Calculate how much shared memory to use per thread
    size_t bytes_per_thread = domain_size_coils * sizeof(complext<REAL>);
    size_t bytes_per_thread_tail = domain_size_coils_tail * sizeof(complext<REAL>);

    unsigned int double_warp_size_power = 0, __tmp = cudaDeviceManager::Instance()->warp_size(device) << 1;
    while (__tmp != 1) {
        __tmp >>= 1;
        double_warp_size_power++;
    }

    vector_td<REAL, D>
    matrix_size_os_real = vector_td<REAL, D>(matrix_size_os);

    // Define temporary image that includes a wrapping zone

    vector<size_t> vec_dims = to_std_vector(matrix_size_os + matrix_size_wrap);
    if (number_of_frames > 1)
        vec_dims.push_back(number_of_frames);
    if (num_batches > 1)
        vec_dims.push_back(num_batches);

    auto tmp_image = cuNDArray<complext<REAL>>(vec_dims);
    //
    // Invoke kernel
    //
    cudaFuncSetCacheConfig(NFFT_H_convolve_kernel<REAL,D>, cudaFuncCachePreferShared);

    for (unsigned int repetition = 0; repetition < num_repetitions; repetition++) {

        size_t num_coils = (repetition == num_repetitions - 1) ? domain_size_coils_tail : domain_size_coils;

        auto samples_view = cuNDArray<complext<REAL>>(std::vector<size_t>{number_of_samples, number_of_frames,num_coils},
                const_cast<complext<REAL>*>(samples->get_data_ptr())+repetition * number_of_samples * number_of_frames * domain_size_coils);

        std::vector<size_t> permuation = {2,0,1};
        auto samples_permuted = permute(samples_view,permuation);

        NFFT_H_convolve_kernel<REAL, D>
                << < dimGrid, dimBlock,
                ((repetition == num_repetitions - 1) ? dimBlock.x * bytes_per_thread_tail : dimBlock.x *
                                                                                            bytes_per_thread) >> >
                (alpha, beta, W, vector_td<unsigned int, D>(matrix_size_os + matrix_size_wrap), number_of_samples,
                        num_coils,
                        raw_pointer_cast(&trajectory_positions[0]),
                        tmp_image.get_data_ptr() +
                        repetition * prod(matrix_size_os + matrix_size_wrap) * number_of_frames * domain_size_coils,
                        samples_permuted.get_data_ptr(),
                        raw_pointer_cast(&tuples_last[0]), raw_pointer_cast(&bucket_begin[0]), raw_pointer_cast(
                        &bucket_end[0]),
                        double_warp_size_power, REAL(0.5) * W, REAL(1) / (W), matrix_size_os_real);

    }

    CHECK_FOR_CUDA_ERROR();

    plan->image_wrap(&tmp_image, image, accumulate);
}


// Image wrap kernels

template<class REAL, unsigned int D>
__global__ void
image_wrap_kernel(typename uintd<D>::Type matrix_size_os, typename uintd<D>::Type matrix_size_wrap, bool accumulate,
                  const complext <REAL> *__restrict__ in, complext <REAL> *__restrict__ out) {
    unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
    const unsigned int num_elements_per_image_src = prod(matrix_size_os + matrix_size_wrap);
    const unsigned int image_offset_src = blockIdx.y * num_elements_per_image_src;

    const typename uintd<D>::Type co = idx_to_co<D>(idx, matrix_size_os);
    const typename uintd<D>::Type half_wrap = matrix_size_wrap >> 1;

    // Make "boolean" vectors denoting whether wrapping needs to be performed in a given direction (forwards/backwards)
    vector_td<bool, D>
    B_l = vector_less(co, half_wrap);
    vector_td<bool, D>
    B_r = vector_greater_equal(co, matrix_size_os - half_wrap);

    complext <REAL> result = in[co_to_idx<D>(co + half_wrap, matrix_size_os + matrix_size_wrap) + image_offset_src];

    if (sum(B_l + B_r) > 0) {

        // Fold back the wrapping zone onto the image ("periodically")
        //
        // There is 2^D-1 ways to pick combinations of dimensions in D-dimensionsal space, e.g.
        //
        //  { x, y, xy } in 2D
        //  { x, y, x, xy, xz, yz, xyz } in 3D
        //
        // Every "letter" in each combination provides two possible wraps (eiher end of the dimension)
        //
        // For every 2^D-1 combinations DO
        //   - find the number of dimensions, d, in the combination
        //   - create 2^(d) stride vectors and test for wrapping using the 'B'-vectors above.
        //   - accumulate the contributions
        //
        //   The following code represents dimensions as bits in a char.
        //

        for (unsigned char combination = 1; combination < (1 << D); combination++) {

            // Find d
            unsigned char d = 0;
            for (unsigned char i = 0; i < D; i++)
                d += ((combination & (1 << i)) > 0);

            // Create stride vector for each wrapping test
            for (unsigned char s = 0; s < (1 << d); s++) {

                // Target for stride
                typename intd<D>::Type stride;
                char wrap_requests = 0;
                char skipped_dims = 0;

                // Fill dimensions of the stride
                for (unsigned char i = 1; i < D + 1; i++) {

                    // Is the stride dimension present in the current combination?
                    if (i & combination) {

                        // A zero bit in s indicates "check for left wrap" and a one bit is interpreted as "check for right wrap"
                        // ("left/right" for the individual dimension meaning wrapping on either side of the dimension).

                        if (i & (s << (skipped_dims))) {
                            if (B_r.vec[i - 1]) { // Wrapping required
                                stride[i - 1] = -1;
                                wrap_requests++;
                            } else
                                stride[i - 1] = 0;
                        } else {
                            if (B_l.vec[i - 1]) { // Wrapping required
                                stride[i - 1] = 1;
                                wrap_requests++;
                            } else
                                stride[i - 1] = 0;
                        }
                    } else {
                        // Do not test for wrapping in dimension 'i-1' (for this combination)
                        stride[i - 1] = 0;
                        skipped_dims++;
                    }
                }

                // Now it is time to do the actual wrapping (if needed)
                if (wrap_requests == d) {
                    typename intd<D>::Type src_co_int = vector_td<int, D>(co + half_wrap);
                    typename intd<D>::Type matrix_size_os_int = vector_td<int, D>(matrix_size_os);
                    typename intd<D>::Type co_offset_int =
                            src_co_int + component_wise_mul<int, D>(stride, matrix_size_os_int);
                    typename uintd<D>::Type co_offset = vector_td<unsigned int, D>(co_offset_int);
                    result += in[co_to_idx<D>(co_offset, matrix_size_os + matrix_size_wrap) + image_offset_src];
                    break; // only one stride per combination can contribute (e.g. one edge, one corner)
                }
            }
        }
    }

    // Output
    const unsigned int image_offset_tgt = blockIdx.y * prod(matrix_size_os);
    if (accumulate) result += out[idx + image_offset_tgt];
    out[idx + image_offset_tgt] = result;
}

template<class REAL, unsigned int D, ConvolutionType CONV>
void
Gadgetron::cuNFFT_impl<REAL, D, CONV>::image_wrap(cuNDArray <complext<REAL>> *source,
                                                  cuNDArray <complext<REAL>> *target, bool accumulate) {
    unsigned int num_batches = 1;
    for (unsigned int d = D; d < source->get_number_of_dimensions(); d++)
        num_batches *= source->get_size(d);
    num_batches /= this->number_of_frames;

    // Set dimensions of grid/blocks.
    unsigned int bdim = 256;
    dim3 dimBlock(bdim);
    dim3 dimGrid(prod(this->matrix_size_os) / bdim, this->number_of_frames * num_batches);

    // Safety check
    if ((prod(this->matrix_size_os) % bdim) != 0) {
        std::stringstream ss;
        ss << "Error: cuNFFT : the number of oversampled image elements must be a multiplum of the block size: "
           << bdim;
        throw std::runtime_error(ss.str());
    }

    // Invoke kernel
    image_wrap_kernel<REAL, D> << < dimGrid, dimBlock >> >
                                             (vector_td<unsigned int, D>(
                                                     this->matrix_size_os), vector_td<unsigned int, D>(
                                                     this->matrix_size_wrap), accumulate, source->get_data_ptr(), target->get_data_ptr());

    CHECK_FOR_CUDA_ERROR();
}




template<class REAL, unsigned int D, ConvolutionType CONV>
void Gadgetron::cuNFFT::convolverC2NC<REAL,D,CONV>::convolve_C2NC(Gadgetron::cuNFFT_impl<REAL, D, CONV> *plan,
                                                     const cuNDArray<complext<REAL>> *image,
                                                     cuNDArray<complext<REAL>> *samples, bool accumulate) {
   // private method - no consistency check. We trust in ourselves.

    unsigned int num_batches = 1;
    for (unsigned int d = D; d < image->get_number_of_dimensions(); d++)
        num_batches *= image->get_size(d);
    num_batches /= plan->number_of_frames;
    auto number_of_frames = plan->number_of_frames;
    auto number_of_samples = plan->number_of_samples;


    /*
      Setup grid and threads
    */

    size_t threads_per_block;
    unsigned int max_coils;

    threads_per_block = NFFT_THREADS_PER_KERNEL;

    if (cudaDeviceManager::Instance()->major_version(plan->device) == 1) {
        max_coils = NFFT_MAX_COILS_COMPUTE_1x;
    } else {
        max_coils = NFFT_MAX_COILS_COMPUTE_2x;
    }

    // We can (only) convolve max_coils batches per run due to shared memory issues.
    unsigned int domain_size_coils_desired = num_batches;
    unsigned int num_repetitions = domain_size_coils_desired / max_coils +
                                   (((domain_size_coils_desired % max_coils) == 0) ? 0 : 1);
    unsigned int domain_size_coils = (num_repetitions == 1) ? domain_size_coils_desired : max_coils;
    unsigned int domain_size_coils_tail = (num_repetitions == 1) ? domain_size_coils_desired :
                                          domain_size_coils_desired - (num_repetitions - 1) * domain_size_coils;

    // Block and Grid dimensions
    dim3 dimBlock((unsigned int) threads_per_block);
    dim3 dimGrid((number_of_samples + dimBlock.x - 1) / dimBlock.x, number_of_frames);

    // Calculate how much shared memory to use per thread
    size_t bytes_per_thread = domain_size_coils * sizeof(complext<REAL>);
    size_t bytes_per_thread_tail = domain_size_coils_tail * sizeof(complext<REAL>);

    unsigned int double_warp_size_power = 0;
    unsigned int __tmp = cudaDeviceManager::Instance()->warp_size(plan->device) << 1;
    while (__tmp != 1) {
        __tmp >>= 1;
        double_warp_size_power++;
    }

    auto matrix_size_os_real = vector_td<REAL, D>(plan->matrix_size_os);
    /*
      Invoke kernel
    */

    for (unsigned int repetition = 0; repetition < num_repetitions; repetition++) {

        size_t num_coils = (repetition == num_repetitions - 1) ? domain_size_coils_tail : domain_size_coils;
        auto image_dims = to_std_vector(plan->matrix_size_os);
        image_dims.push_back(number_of_frames);
        image_dims.push_back(num_coils);
        size_t image_view_elements = std::accumulate(image_dims.begin(),image_dims.end(),size_t(1),std::multiplies<size_t>());
        auto image_view = cuNDArray<complext<REAL>>(image_dims,
                                                    const_cast<complext<REAL>*>(image->get_data_ptr())+repetition*image_view_elements);

        auto permutation = std::vector<size_t>(D+2);
        permutation[0] = D+1;
        std::iota(permutation.begin()+1,permutation.end(),0);

        auto image_permuted  = permute(image_view,permutation);


        NFFT_convolve_kernel<REAL, D>
                << < dimGrid, dimBlock,
                ((repetition == num_repetitions - 1) ? dimBlock.x * bytes_per_thread_tail : dimBlock.x *
                                                                                            bytes_per_thread) >> >
                (plan->alpha, plan->beta, plan->W, vector_td<unsigned int, D>(plan->matrix_size_os), vector_td<unsigned int, D>(
                        plan->matrix_size_wrap), number_of_samples,
                        num_coils,
                        raw_pointer_cast(&plan->trajectory_positions[0]),
                        image_permuted.get_data_ptr(),
                        samples->get_data_ptr() + repetition * number_of_samples * number_of_frames * domain_size_coils,
                        double_warp_size_power, REAL(0.5) * plan->W, REAL(1) /
                                                                     (plan->W), accumulate, matrix_size_os_real);

        CHECK_FOR_CUDA_ERROR();
    }
}

namespace Gadgetron {
    template<class REAL, unsigned int D>
    boost::shared_ptr<cuNFFT_plan<REAL, D>> NFFT<cuNDArray, REAL, D>::make_plan(const vector_td<size_t,D>& matrix_size, const vector_td<size_t,D>& matrix_size_os, REAL W, ConvolutionType conv) {
        switch (conv) {
            case ConvolutionType::STANDARD:
                return boost::make_shared<cuNFFT_impl<REAL, D, ConvolutionType::STANDARD> >(matrix_size,matrix_size_os,W);
            case ConvolutionType::ATOMIC:
                return boost::make_shared<cuNFFT_impl<REAL, D, ConvolutionType::ATOMIC> >(matrix_size,matrix_size_os,W);
            case ConvolutionType::SPARSE_MATRIX:
                return boost::make_shared<cuNFFT_impl<REAL, D, ConvolutionType::SPARSE_MATRIX>>(matrix_size,matrix_size_os,W);
        }
        throw std::runtime_error(
                "Invalid convolution type provided. If you're reading this, you may have broken your computer quite badly");
    }

    template<unsigned int D>

    boost::shared_ptr<cuNFFT_plan<double, D>> NFFT<cuNDArray, double, D>::make_plan(const vector_td<size_t,D>& matrix_size, const vector_td<size_t,D>& matrix_size_os, double W, ConvolutionType conv) {
        if (conv == ConvolutionType::STANDARD) {
            return boost::make_shared<cuNFFT_impl<double, D, ConvolutionType::STANDARD>>(matrix_size,matrix_size_os,W);
        }
        throw std::runtime_error("Only standard convolution type supported for doubles");
    }
}

//
// Template instantion
//

template
class EXPORTGPUNFFT Gadgetron::cuNFFT_impl<float, 1, ConvolutionType::ATOMIC>;


template
class EXPORTGPUNFFT Gadgetron::cuNFFT_impl<float, 1, ConvolutionType::SPARSE_MATRIX>;

template
class EXPORTGPUNFFT Gadgetron::cuNFFT_impl<float, 1>;

template
class EXPORTGPUNFFT Gadgetron::cuNFFT_impl<double, 1>;

template
class EXPORTGPUNFFT Gadgetron::cuNFFT_impl<float, 2, ConvolutionType::ATOMIC>;


template
class EXPORTGPUNFFT Gadgetron::cuNFFT_impl<float, 2, ConvolutionType::SPARSE_MATRIX>;

template
class EXPORTGPUNFFT Gadgetron::cuNFFT_impl<float, 2>;

template
class EXPORTGPUNFFT Gadgetron::cuNFFT_impl<double, 2>;

template
class EXPORTGPUNFFT Gadgetron::cuNFFT_impl<float, 3, ConvolutionType::ATOMIC>;

template
class EXPORTGPUNFFT Gadgetron::cuNFFT_impl<float, 3, ConvolutionType::SPARSE_MATRIX>;

template
class EXPORTGPUNFFT Gadgetron::cuNFFT_impl<float, 3>;

template
class EXPORTGPUNFFT Gadgetron::cuNFFT_impl<double, 3>;

template
class EXPORTGPUNFFT Gadgetron::cuNFFT_impl<float, 4, ConvolutionType::ATOMIC>;


template
class EXPORTGPUNFFT Gadgetron::cuNFFT_impl<float, 4, ConvolutionType::SPARSE_MATRIX>;
template
class EXPORTGPUNFFT Gadgetron::cuNFFT_impl<float, 4>;

template
class EXPORTGPUNFFT Gadgetron::cuNFFT_impl<double, 4>;


template class EXPORTGPUNFFT Gadgetron::NFFT<cuNDArray,float,1>;
template class EXPORTGPUNFFT Gadgetron::NFFT<cuNDArray,float,2>;
template class EXPORTGPUNFFT Gadgetron::NFFT<cuNDArray,float,3>;
template class EXPORTGPUNFFT Gadgetron::NFFT<cuNDArray,float,4>;



template class EXPORTGPUNFFT Gadgetron::NFFT<cuNDArray,double,1>;
template class EXPORTGPUNFFT Gadgetron::NFFT<cuNDArray,double,2>;
template class EXPORTGPUNFFT Gadgetron::NFFT<cuNDArray,double,3>;
template class EXPORTGPUNFFT Gadgetron::NFFT<cuNDArray,double,4>;
