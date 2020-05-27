/**
 * \file cuSDC.cu
 * \brief Sampling density compensation (CUDA specialization).
 */

#include "cuSDC.h"
#include "cuSDC_kernel.h"

#include "cudaDeviceManager.h"
#include "cuNDArray.h"
#include "cuNDArray_elemwise.h"
#include "cuNDArray_utils.h"
#include "vector_td.h"

#include <numeric>

// CUDA configuration  
#define SDC_MAX_COILS           16
#define SDC_THREADS_PER_KERNEL  192

// Reference to shared memory
extern __shared__ char _shared_mem[];

#include "cuSDC_conv_C2NC.cu"
#include "cuSDC_conv_NC2C.cu"


namespace Gadgetron
{
    namespace
    {
        template<class REAL, unsigned int D>
        struct traj_scale
        {
            traj_scale(const typename reald<REAL,D>::Type& m, const typename reald<REAL,D>::Type& b, const typename reald<REAL,D>::Type& w)
            {
                matrix = m;
                bias = b;
                win_len = w;
            }
            
            __host__ __device__
            typename reald<REAL,D>::Type operator()(const typename reald<REAL,D>::Type &in) const
            { 
                return component_wise_mul<REAL, D>(component_wise_mul<REAL, D>(in, win_len), matrix) + bias;
            }
            
            typename reald<REAL, D>::Type matrix;
            typename reald<REAL, D>::Type bias;
            typename reald<REAL, D>::Type win_len;
        };
    }
        
    template<class REAL, unsigned int D>
    cuSDC_impl<REAL, D>::cuSDC_impl(const vector_td<size_t, D>& matrix_size, REAL os_factor, size_t num_iterations, int device)
      : SDC_impl<cuNDArray, REAL, D>(matrix_size, os_factor, num_iterations)
    {
        // minimal initialization
        barebones();

        // check if the device is valid
        if (device < 0)
        {
            if (cudaGetDevice(&this->device_) != cudaSuccess)
            {
                throw cuda_error("Error: cuSDC_impl::setup: unable to determine this->device properties.");
            }
        }
        else
        {
            this->device_ = device;
        }

        // make sure grid size is a multiple of the warp size
        size_t warp_size = static_cast<size_t>(cudaDeviceManager::Instance()->warp_size(this->device_));
        typename uint64d<D>::Type warp_size_vec(warp_size);
        this->grid_size_ = ((this->grid_size_ + warp_size - size_t(1)) / warp_size) * warp_size;

        // initialize kernel
        this->kernel_ = cuSDC_kernel<REAL, D>(this->matrix_size_, this->grid_size_);

        // check input against certain requirements
        if (sum(matrix_size % warp_size_vec) || sum(this->grid_size_ % warp_size_vec))
        {
            throw std::runtime_error("Error: Illegal matrix size for the cuSDC plan (not a multiple of the warp size)");
        }

        // set up private variables
        this->grid_padding_ = vector_td<size_t, D>(ceil(vector_td<REAL, D>(this->kernel_.get_rfp())));
        this->grid_padding_ <<= 1;

        int device_no_old;
        if (cudaGetDevice(&device_no_old) != cudaSuccess)
        {
            throw cuda_error("Error: cuSDC_impl::setup: unable to get this->device no");
        }

        if (this->device_ != device_no_old && cudaSetDevice(this->device_) != cudaSuccess)
        {
            throw cuda_error("Error: cuSDC_impl::setup: unable to set this->device");
        }

        if (this->device_ != device_no_old && cudaSetDevice(device_no_old) != cudaSuccess)
        {
            throw cuda_error("Error: cuSDC_impl::setup: unable to restore this->device");
        }
    }


    template<class REAL, unsigned int D>
    void cuSDC_impl<REAL, D>::preprocess(const cuNDArray<vector_td<REAL, D>>& traj)
    {
        // call base class
        SDC_impl<cuNDArray, REAL, D>::preprocess(traj);
        cuNDArray<REAL> trajectory_view(std::vector<size_t>{traj.get_number_of_elements() * D}, (REAL*)traj.get_data_ptr());

        // Make sure that the traj values are within range [-1/2;1/2]
        thrust::pair<thrust::device_ptr<REAL>, thrust::device_ptr<REAL> > mm_pair =
                thrust::minmax_element(trajectory_view.begin(), trajectory_view.end());

        if (*mm_pair.first < REAL(-0.5) || *mm_pair.second > REAL(0.5)) {
            std::stringstream ss;
            ss << "Error: cuNFFT::preprocess : traj [" << *mm_pair.first << "; " << *mm_pair.second
            << "] out of range [-1/2;1/2]";
            throw std::runtime_error(ss.str());
        }

        traj_ = thrust::device_vector<vector_td<REAL, D>>(traj.get_number_of_elements());

        CHECK_FOR_CUDA_ERROR();

        vector_td<REAL, D> grid_size = vector_td<REAL, D>(this->grid_size_);
        vector_td<REAL, D> half_padded_grid_size = vector_td<REAL, D>((this->grid_size_ + this->grid_padding_) >> 1);
        vector_td<REAL, D> win_length = (grid_size + this->kernel_.get_rfp() * REAL(2)) / grid_size;

        // convert input traj in [-1/2;1/2] to [0;this->grid_size_]
        thrust::transform(traj.begin(), traj.end(), traj_.begin(), traj_scale<REAL, D>(grid_size, half_padded_grid_size, win_length));

        CHECK_FOR_CUDA_ERROR();

    }

    template<class REAL, unsigned int D>
    void cuSDC_impl<REAL, D>::convolve(const cuNDArray<REAL>& in, cuNDArray<REAL>& out, SDC_conv_mode mode)
    {
        {
            const cuNDArray <REAL> *samples, *image;

            if (mode == SDC_conv_mode::C2NC)
            {
                image = &in;
                samples = &out;
            }
            else
            {
                image = &out;
                samples = &in;
            }

            check_consistency(samples, image, 0x0);
        }

        typename uint64d<D>::Type image_dims = from_std_vector<size_t, D>
                (*(((mode == SDC_conv_mode::C2NC) ? in : out).get_dimensions()));
        bool oversampled_image = (image_dims == this->grid_size_);

        if (!oversampled_image)
        {
            throw std::runtime_error("cuSDC_impl::convolve(): input matrix has unexpected size.");
        }

        switch (mode)
        {
            case SDC_conv_mode::C2NC:   convolve_C2NC(&in, &out);   break;
            case SDC_conv_mode::NC2C:   convolve_NC2C(&in, &out);   break;
        }
    }


    template<class REAL, unsigned int D>
    void cuSDC_impl<REAL, D>::update(const cuNDArray<REAL>& in, cuNDArray<REAL>& out)
    {
        thrust::transform(out.begin(), out.end(), in.begin(), out.begin(), SDC_internal::safe_divides<REAL>());
    }


    template<class REAL, unsigned int D>
    void cuSDC_impl<REAL, D>::convolve_NC2C(const cuNDArray<REAL>* in, cuNDArray<REAL>* out)
    {
        // atomic operations are only supported in compute model 2.0 and up
        if (cudaDeviceManager::Instance()->major_version(this->device_) == 1)
        {
            throw cuda_error("Error: Atomic NC2C NFFT only supported on this->device with compute model 2.0 or higher");
        }

        // check if warp_size is a power of two. We do some modulus tricks in the kernels that depend on this...
        if (!((cudaDeviceManager::Instance()->warp_size(this->device_) & (cudaDeviceManager::Instance()->warp_size(this->device_) - 1)) == 0))
        {
            throw cuda_error("cuNFFT: unsupported hardware (warpSize is not a power of two)");
        }

        unsigned int num_batches = 1;
        for (unsigned int d = D; d < out->get_number_of_dimensions(); d++)
            num_batches *= out->get_size(d);
        num_batches /= this->num_frames_;

        //  setup grid and threads
        unsigned int threads_per_block = SDC_THREADS_PER_KERNEL;
        unsigned int max_coils = SDC_MAX_COILS;

        // block and grid dimensions
        dim3 dimBlock(threads_per_block);
        dim3 dimGrid((this->num_samples_ + dimBlock.x - 1) / dimBlock.x, this->num_frames_);

        // we can (only) convolve domain_size_coils batches per run due to shared memory issues. 
        unsigned int domain_size_coils_desired = num_batches;
        unsigned int num_repetitions = domain_size_coils_desired / max_coils + (((domain_size_coils_desired % max_coils) == 0) ? 0 : 1);
        unsigned int domain_size_coils = (num_repetitions == 1) ? domain_size_coils_desired : max_coils;
        unsigned int domain_size_coils_tail = (num_repetitions == 1) ? domain_size_coils_desired : domain_size_coils_desired - (num_repetitions - 1) * domain_size_coils;

        // calculate how much shared memory to use per thread
        size_t bytes_per_thread = domain_size_coils * sizeof(vector_td<REAL, D>);
        size_t bytes_per_thread_tail = domain_size_coils_tail * sizeof(vector_td<REAL, D>);

        clear(out);

        for (unsigned int repetition = 0; repetition < num_repetitions; repetition++)
        {
            unsigned int num_samples = static_cast<unsigned int>(this->num_samples_);
            unsigned int num_coils = (repetition == num_repetitions - 1) ? domain_size_coils_tail : domain_size_coils;
            const vector_td<REAL, D>* traj_ptr = raw_pointer_cast(&this->traj_[0]);
            const REAL* in_ptr = in->get_data_ptr() + repetition * this->num_samples_ * this->num_frames_ * domain_size_coils;
            REAL* out_ptr = out->get_data_ptr() + repetition * prod(this->grid_size_) * this->num_frames_ * domain_size_coils;

            size_t sharedMemSize = ((repetition == num_repetitions - 1) ? dimBlock.x * bytes_per_thread_tail : dimBlock.x * bytes_per_thread);
            
            SDC_internal::convolve_NC2C_kernel<REAL, D>
                <<<dimGrid, dimBlock, sharedMemSize>>>(
                    this->kernel_.get_rfp(), vector_td<unsigned int, D>(this->grid_size_), vector_td<unsigned int, D>(this->grid_padding_),
                    num_samples, num_coils, traj_ptr, in_ptr, out_ptr);
        }

        CHECK_FOR_CUDA_ERROR();
    }


    template<class REAL, unsigned int D>
    void cuSDC_impl<REAL, D>::convolve_C2NC(const cuNDArray<REAL>* in, cuNDArray<REAL>* out)
    {
        unsigned int num_batches = 1;
        for (unsigned int d = D; d < in->get_number_of_dimensions(); d++)
            num_batches *= in->get_size(d);
        num_batches /= this->num_frames_;

        // setup grid and threads
        unsigned int threads_per_block = SDC_THREADS_PER_KERNEL;
        unsigned int max_coils = SDC_MAX_COILS;

        // block and grid dimensions
        dim3 dimBlock(threads_per_block);
        dim3 dimGrid((this->num_samples_ + dimBlock.x - 1) / dimBlock.x, this->num_frames_);

        // we can (only) convolve max_coils batches per run due to shared memory issues.
        unsigned int domain_size_coils_desired = num_batches;
        unsigned int num_repetitions = domain_size_coils_desired / max_coils + (((domain_size_coils_desired % max_coils) == 0) ? 0 : 1);
        unsigned int domain_size_coils = (num_repetitions == 1) ? domain_size_coils_desired : max_coils;
        unsigned int domain_size_coils_tail = (num_repetitions == 1) ? domain_size_coils_desired : domain_size_coils_desired - (num_repetitions - 1) * domain_size_coils;

        // calculate how much shared memory to use per thread
        size_t bytes_per_thread = domain_size_coils * sizeof(REAL);
        size_t bytes_per_thread_tail = domain_size_coils_tail * sizeof(REAL);

        unsigned int warp_size_power = 0;
        unsigned int __tmp = cudaDeviceManager::Instance()->warp_size(this->device_);
        while (__tmp != 1)
        {
            __tmp >>= 1;
            warp_size_power++;
        }
        
        // invoke kernel
        for (unsigned int repetition = 0; repetition < num_repetitions; repetition++)
        {
            unsigned int num_coils = (repetition == num_repetitions - 1) ? domain_size_coils_tail : domain_size_coils;
            auto image_dims = to_std_vector(this->grid_size_);
            image_dims.push_back(this->num_frames_);
            image_dims.push_back(num_coils);
            size_t image_view_elements = std::accumulate(image_dims.begin(), image_dims.end(), size_t(1), std::multiplies<size_t>());
            auto image_view = cuNDArray<REAL>(image_dims, const_cast<REAL*>(in->get_data_ptr()) + repetition * image_view_elements);

            auto permutation = std::vector<size_t>(D+2);
            permutation[0] = D+1;
            std::iota(permutation.begin() + 1, permutation.end(), 0);

            auto image_permuted = permute(image_view, permutation);

            unsigned int num_samples = static_cast<unsigned int>(this->num_samples_);
            const vector_td<REAL, D>* traj_ptr = raw_pointer_cast(&this->traj_[0]);
            const REAL* in_ptr = image_permuted.get_data_ptr();
            REAL* out_ptr = out->get_data_ptr() + repetition * this->num_samples_ * this->num_frames_ * domain_size_coils;

            size_t sharedMemSize = ((repetition == num_repetitions - 1) ? dimBlock.x * bytes_per_thread_tail : dimBlock.x * bytes_per_thread);

            SDC_internal::convolve_C2NC_kernel<REAL, D>
                <<<dimGrid, dimBlock, sharedMemSize>>>(
                    this->kernel_.get_rfp(), vector_td<unsigned int, D>(this->grid_size_), vector_td<unsigned int, D>(this->grid_padding_),
                    num_samples, num_coils, traj_ptr, in_ptr, out_ptr, warp_size_power);

            CHECK_FOR_CUDA_ERROR();
        }
    }

    template<class REAL, unsigned int D>
    void cuSDC_impl<REAL, D>::check_consistency(const cuNDArray<REAL> *samples, const cuNDArray<REAL> *image, const cuNDArray<REAL> *weights)
    {

        if (image->get_number_of_dimensions() < D) {
            throw std::runtime_error("Error: cuSDC_impl: Number of image dimensions mismatch the plan.");
        }

        typename uint64d<D>::Type image_dims = from_std_vector<size_t, D>(*image->get_dimensions());
        bool oversampled_image = (image_dims == this->grid_size_);

        if (!((oversampled_image) ? (image_dims == this->grid_size_) : (image_dims == this->matrix_size_)))
        {
            throw std::runtime_error("Error: cuSDC_impl: Image dimensions mismatch.");
        }

        if ((samples->get_number_of_elements() == 0) || (samples->get_number_of_elements() % (this->num_frames_ * this->num_samples_)))
        {
            printf("\ncuNFFT::check_consistency() failed:\n#elements in the samples array: %d.\n#samples from preprocessing: %d.\n#frames from preprocessing: %d.\n",
                (int)samples->get_number_of_elements(), (int)this->num_samples_, (int)this->num_frames_);
            fflush(stdout);
            throw std::runtime_error("Error: cuSDC_impl: The number of samples is not a multiple of #samples/frame x #frames as requested through preprocessing");
        }

        unsigned int num_batches_in_samples_array = samples->get_number_of_elements() / (this->num_frames_ * this->num_samples_);
        unsigned int num_batches_in_image_array = 1;

        for (unsigned int d = D; d < image->get_number_of_dimensions(); d++)
        {
            num_batches_in_image_array *= image->get_size(d);
        }

        num_batches_in_image_array /= this->num_frames_;

        if (num_batches_in_samples_array != num_batches_in_image_array)
        {
            printf("\ncuNFFT::check_consistency() failed:\n#elements in the samples array: %d.\n#samples from preprocessing: %d.\n#frames from preprocessing: %d.\nLeading to %d batches in the samples array.\nThe number of batches in the image array is %d.\n",
                (int)samples->get_number_of_elements(), (int)this->num_samples_, (int)this->num_frames_, (int)num_batches_in_samples_array, (int)num_batches_in_image_array);
            fflush(stdout);
            throw std::runtime_error("Error: cuSDC_impl: Number of batches mismatch between samples and image arrays");
        }
    }


    template<class REAL, unsigned int D>
    void cuSDC_impl<REAL, D>::barebones()
    {
        // and specify the this->device
        if (cudaGetDevice(&this->device_) != cudaSuccess)
        {
            throw cuda_error("Error: cuSDC_impl::barebones:: unable to get this->device no");
        }
    }

    template<class REAL, unsigned int D>
    boost::shared_ptr<cuNDArray<REAL>> estimate_dcw(
        const cuNDArray<vector_td<REAL, D>>& traj,
        const vector_td<size_t, D>& matrix_size,
        REAL os_factor,
        size_t num_iterations)
    { 
        cuSDC_impl<REAL, D> impl(matrix_size, os_factor, num_iterations);
        return impl.compute(traj);
    }

    template<class REAL, unsigned int D>
    boost::shared_ptr<cuNDArray<REAL>> estimate_dcw(
        const cuNDArray<vector_td<REAL, D>>& traj,
        const cuNDArray<REAL>& initial_dcw,
        const vector_td<size_t, D>& matrix_size,
        REAL os_factor,
        size_t num_iterations)
    { 
        cuSDC_impl<REAL, D> impl(matrix_size, os_factor, num_iterations);
        return impl.compute(traj, initial_dcw);
    }

}   // namespace Gadgetron


template EXPORTSDC boost::shared_ptr<Gadgetron::cuNDArray<float>> Gadgetron::estimate_dcw<float, 2>(
    const Gadgetron::cuNDArray<Gadgetron::vector_td<float, 2>>& traj,
    const Gadgetron::vector_td<size_t, 2>& matrix_size,
    float os_factor,
    size_t num_iterations);

template EXPORTSDC boost::shared_ptr<Gadgetron::cuNDArray<float>> Gadgetron::estimate_dcw<float, 3>(
    const Gadgetron::cuNDArray<Gadgetron::vector_td<float, 3>>& traj,
    const Gadgetron::vector_td<size_t, 3>& matrix_size,
    float os_factor,
    size_t num_iterations);

template EXPORTSDC boost::shared_ptr<Gadgetron::cuNDArray<float>> Gadgetron::estimate_dcw<float, 2>(
    const Gadgetron::cuNDArray<Gadgetron::vector_td<float, 2>>& traj,
    const Gadgetron::cuNDArray<float>& initial_dcw,
    const Gadgetron::vector_td<size_t, 2>& matrix_size,
    float os_factor,
    size_t num_iterations);

template EXPORTSDC boost::shared_ptr<Gadgetron::cuNDArray<float>> Gadgetron::estimate_dcw<float, 3>(
    const Gadgetron::cuNDArray<Gadgetron::vector_td<float, 3>>& traj,
    const Gadgetron::cuNDArray<float>& initial_dcw,
    const Gadgetron::vector_td<size_t, 3>& matrix_size,
    float os_factor,
    size_t num_iterations);
