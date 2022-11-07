
#include "cuGriddingConvolution.h"

#include <thrust/binary_search.h>

#include "cuNDArray_elemwise.h"
#include "cuNDArray_utils.h"
#include "cudaDeviceManager.h"

#include "ConvolverC2NC.cuh"
#include "ConvolverNC2C_atomic.cuh"
#include "ConvolverNC2C_sparse.cuh"
#include "ConvolverNC2C_standard.cuh"

#define CUDA_CONV_MAX_COILS             (16)
#define CUDA_CONV_THREADS_PER_KERNEL    (192)

namespace Gadgetron
{
    template<class T, unsigned int D, template<class, unsigned int> class K>
    cuGriddingConvolution<T, D, K>::cuGriddingConvolution(
        const vector_td<size_t, D>& matrix_size,
        const vector_td<size_t, D>& matrix_size_os,
        const K<REAL, D>& kernel,
        ConvolutionType conv_type)
      : GriddingConvolutionBase<cuNDArray, T, D, K>(
            matrix_size, matrix_size_os, kernel)
    {
        this->initialize(conv_type);
    }


    template<class T, unsigned int D, template<class, unsigned int> class K>
    cuGriddingConvolution<T, D, K>::cuGriddingConvolution(
        const vector_td<size_t, D>& matrix_size,
        const REAL os_factor,
        const K<REAL, D>& kernel,
        ConvolutionType conv_type)
      : GriddingConvolutionBase<cuNDArray, T, D, K>(
            matrix_size, os_factor, kernel)
    {
        this->initialize(conv_type);
    }


    template<class T, unsigned int D, template<class, unsigned int> class K>
    cuGriddingConvolution<T, D, K>::~cuGriddingConvolution()
    {
        // Release kernel on device.
        cudaFree(this->d_kernel_);
    }


    template<class T, unsigned int D, template<class, unsigned int> class K>
    void cuGriddingConvolution<T, D, K>::initialize(ConvolutionType conv_type)
    {
        // Get device number.
        if (cudaGetDevice(&this->device_) != cudaSuccess)
        {
            throw cuda_error("Could not get device number.");
        }

        // The convolution does not work properly for very small convolution
        // kernel width (experimentally observed limit).
        if (this->kernel_.get_width() < REAL(1.8))
            throw std::runtime_error("Kernel width must be equal to or larger "
                                     "than 1.8.");

        // Matrix size must be a multiple of warp size.
        vector_td<size_t, D> warp_size(
            (size_t)cudaDeviceManager::Instance()->warp_size(this->device_));
        if (sum(this->matrix_size_ % warp_size) ||
            sum(this->matrix_size_os_ % warp_size))
        {
            std::stringstream ss;
            ss << "Matrix size must be a multiple of " << warp_size[0] << ".";
            throw std::runtime_error(ss.str());
        }
        
        // Compute matrix padding.
        vector_td<REAL, D> radius(this->kernel_.get_radius());
        this->matrix_padding_ = vector_td<size_t, D>(ceil(radius));
        this->matrix_padding_ <<= 1;

        // Compute warp size power.
        unsigned int warp_size_power = 0;
        unsigned int tmp = cudaDeviceManager::Instance()->warp_size(
            this->device_);
        while (tmp != 1)
        {
            tmp >>= 1;
            warp_size_power++;
        }
        this->warp_size_power_ = warp_size_power;

        // Copy kernel to device.
        cudaMalloc((void**)&this->d_kernel_, sizeof(K<REAL, D>));
        cudaMemcpy((void*)this->d_kernel_, (const void*)&this->kernel_,
                   sizeof(K<REAL, D>), cudaMemcpyHostToDevice);

        // Set up convolvers.
        switch (conv_type)
        {
            case ConvolutionType::STANDARD:
            {
                this->conv_C2NC_ = std::make_unique<ConvolverC2NC<T, D, K, ConvolutionType::STANDARD>>(*this);
                this->conv_NC2C_ = std::make_unique<ConvolverNC2C<T, D, K, ConvolutionType::STANDARD>>(*this);
                break;
            }
            case ConvolutionType::ATOMIC:
            {
                this->conv_C2NC_ = std::make_unique<ConvolverC2NC<T, D, K, ConvolutionType::ATOMIC>>(*this);
                this->conv_NC2C_ = std::make_unique<ConvolverNC2C<T, D, K, ConvolutionType::ATOMIC>>(*this);
                break;
            }
            case ConvolutionType::SPARSE_MATRIX:
            {
                this->conv_C2NC_ = std::make_unique<ConvolverC2NC<T, D, K, ConvolutionType::SPARSE_MATRIX>>(*this);
                this->conv_NC2C_ = std::make_unique<ConvolverNC2C<T, D, K, ConvolutionType::SPARSE_MATRIX>>(*this);
                break;
            }
        }
    }


    template<class T, unsigned int D, template<class, unsigned int> class K>
    void cuGriddingConvolution<T, D, K>::preprocess(
        const cuNDArray<vector_td<REAL, D>>& trajectory,
        GriddingConvolutionPrepMode prep_mode)
    {
        // Call base class method.
        GriddingConvolutionBase<cuNDArray, T, D, K>::preprocess(
            trajectory, prep_mode);
    
        // Make sure that the trajectory values are within range [-1/2, 1/2].
        cuNDArray<REAL> traj_view(
            std::vector<size_t>{trajectory.get_number_of_elements() * D},
            (REAL*) trajectory.get_data_ptr());

        thrust::pair<thrust::device_ptr<REAL>, thrust::device_ptr<REAL>> mm_pair =
            thrust::minmax_element(traj_view.begin(), traj_view.end());
        
        if (*mm_pair.first < REAL(-0.5) || *mm_pair.second > REAL(0.5))
        {
            std::stringstream ss;
            ss << "Error: cuGriddingConvolution::preprocess: trajectory [" <<
                *mm_pair.first << ", " << *mm_pair.second <<
                "] out of range [-1/2, 1/2].";
            throw std::runtime_error(ss.str());
        }
        
        // Allocate trajectory.
        this->trajectory_ = thrust::device_vector<vector_td<REAL, D>>(
            trajectory.get_number_of_elements());

        CHECK_FOR_CUDA_ERROR();
        
        // Cast matrix size values to floating-point type.
        vector_td<REAL, D> matrix_size_os_fp =
            vector_td<REAL, D>(this->matrix_size_os_);
        vector_td<REAL, D> matrix_size_os_padded_fp =
            vector_td<REAL, D>((this->matrix_size_os_ + this->matrix_padding_) >> 1);
    
        // Convert input trajectory from range [-1/2, 1/2] to
        // [0, this->matrix_size_os_], and copy to class member.
        thrust::transform(trajectory.begin(),
                          trajectory.end(),
                          this->trajectory_.begin(),
                          trajectory_scale<REAL, D>(matrix_size_os_fp,
                                                    matrix_size_os_padded_fp));
        
        // Prepare convolution.
        if (prep_mode == GriddingConvolutionPrepMode::C2NC ||
            prep_mode == GriddingConvolutionPrepMode::ALL)
        {
            this->conv_C2NC_->prepare(this->trajectory_);
            CHECK_FOR_CUDA_ERROR(); 
        }

        if (prep_mode == GriddingConvolutionPrepMode::NC2C ||
            prep_mode == GriddingConvolutionPrepMode::ALL)
        {
            this->conv_NC2C_->prepare(this->trajectory_);
            CHECK_FOR_CUDA_ERROR(); 
        }
    }


    template<class T, unsigned int D, template<class, unsigned int> class K>
    void cuGriddingConvolution<T, D, K>::compute_C2NC(
        const cuNDArray<T>& image,
        cuNDArray<T>& samples,
        bool accumulate)
    {
        this->check_inputs(samples, image);
        this->conv_C2NC_->compute(image, samples, accumulate);
    }


    template<class T, unsigned int D, template<class, unsigned int> class K>
    void cuGriddingConvolution<T, D, K>::compute_NC2C(
        const cuNDArray<T>& samples,
        cuNDArray<T>& image,
        bool accumulate)
    {
        this->check_inputs(samples, image);
        this->conv_NC2C_->compute(samples, image, accumulate);
    }


    template<class T, unsigned int D, template<class, unsigned int> class K>
    void cuGriddingConvolution<T, D, K>::check_inputs(
        const cuNDArray<T>& samples,
        const cuNDArray<T>& image)
    {
        if (image.get_number_of_dimensions() < D)
        {
            throw std::runtime_error("Number of image dimensions does not "
                                     "match the plan.");
        }

        vector_td<size_t, D> image_dims =
            from_std_vector<size_t, D>(*image.get_dimensions());
        if (image_dims != this->matrix_size_os_)
            throw std::runtime_error("Image dimensions mismatch.");

        if ((samples.get_number_of_elements() == 0) ||
            (samples.get_number_of_elements() % (this->num_frames_ * this->num_samples_)))
        {
            printf("\nConsistency check failed:\n"
                "#elements in the samples array: %ld.\n"
                "#samples from preprocessing: %zu.\n"
                "#frames from preprocessing: %zu.\n",
                samples.get_number_of_elements(), this->num_samples_, this->num_frames_);
            fflush(stdout);
            throw std::runtime_error(
                "The number of samples is not a multiple of #samples/frame x "
                "#frames as requested through preprocessing.");
        }

        unsigned int num_batches_in_samples_array =
            samples.get_number_of_elements() / (this->num_frames_ * this->num_samples_);
        unsigned int num_batches_in_image_array = 1;

        for (unsigned int d = D; d < image.get_number_of_dimensions(); d++) {
            num_batches_in_image_array *= image.get_size(d);
        }
        num_batches_in_image_array /= this->num_frames_;

        if (num_batches_in_samples_array != num_batches_in_image_array) {
            printf("\nConsistency check failed:\n"
                "#elements in the samples array: %ld.\n"
                "#samples from preprocessing: %zu.\n"
                "#frames from preprocessing: %zu.\n"
                "Leading to %d batches in the samples array.\n"
                "The number of batches in the image array is %u.\n",
                samples.get_number_of_elements(),
                this->num_samples_,
                this->num_frames_,
                num_batches_in_samples_array,
                num_batches_in_image_array);
            fflush(stdout);
            throw std::runtime_error(
                "Number of batches mismatch between samples and image arrays.");
        }   
    }


    template<class T, unsigned int D, template<class, unsigned int> class K, ConvolutionType C>
    void ConvolverC2NC<T, D, K, C>::prepare(
        const thrust::device_vector<vector_td<REAL, D>>& trajectory)
    {
        // No need to do anything here.
        // Defined for completeness.
    }


    template<class T, unsigned int D, template<class, unsigned int> class K, ConvolutionType C>
    void ConvolverC2NC<T, D, K, C>::compute(
        const cuNDArray<T>& image,
        cuNDArray<T>& samples,
        bool accumulate)
    {
        // Compute number of batches.
        unsigned int num_batches = 1;
        for (unsigned int d = D; d < image.get_number_of_dimensions(); d++)
            num_batches *= image.get_size(d);
        num_batches /= this->plan_.num_frames_;

        // Set up grid and threads. We can convolve only max_coils batches per
        // run due to shared memory issues.
        unsigned int threads_per_block = CUDA_CONV_THREADS_PER_KERNEL;
        unsigned int max_coils = CUDA_CONV_MAX_COILS;
        unsigned int domain_size_coils_desired = num_batches;
        unsigned int num_repetitions = domain_size_coils_desired / max_coils +
            (((domain_size_coils_desired % max_coils) == 0) ? 0 : 1);
        unsigned int domain_size_coils =
            (num_repetitions == 1) ? domain_size_coils_desired : max_coils;
        unsigned int domain_size_coils_tail =
            (num_repetitions == 1) ? domain_size_coils_desired :
            domain_size_coils_desired - (num_repetitions - 1) * domain_size_coils;

        // Block and grid dimensions.
        dim3 dimBlock(threads_per_block);
        dim3 dimGrid((this->plan_.num_samples_ + dimBlock.x - 1) / dimBlock.x,
                     this->plan_.num_frames_);

        // Calculate how much shared memory to use per thread.
        size_t bytes_per_thread = domain_size_coils * sizeof(T);
        size_t bytes_per_thread_tail = domain_size_coils_tail * sizeof(T);
        
        // Image view dimensions.
        auto view_dims = to_std_vector(this->plan_.matrix_size_os_);
        view_dims.push_back(this->plan_.num_frames_);
        view_dims.push_back(0); // Placeholder for num_coils.

        for (unsigned int repetition = 0; repetition < num_repetitions; repetition++)
        {
            // Number of coils in this repetition.
            size_t num_coils = (repetition == num_repetitions - 1) ?
                domain_size_coils_tail : domain_size_coils;
            view_dims.back() = num_coils; 
            
            // Image view for this repetition.
            size_t image_view_elements = std::accumulate(
                view_dims.begin(), view_dims.end()-1, size_t(1), std::multiplies<size_t>());
            auto image_view = cuNDArray<T>(view_dims, const_cast<T*>(
                image.get_data_ptr()) + repetition * image_view_elements * domain_size_coils);
            auto permutation = std::vector<size_t>(D+2);
            permutation[0] = D+1;
            std::iota(permutation.begin() + 1, permutation.end(), 0);
            auto image_permuted = permute(image_view, permutation);
            
            // Size of shared memory.
            size_t sharedMemSize = (repetition == num_repetitions - 1) ?
                dimBlock.x * bytes_per_thread_tail :
                dimBlock.x * bytes_per_thread;

            // Launch CUDA kernel.
            NFFT_convolve_kernel<T, D, K>
                <<<dimGrid, dimBlock, sharedMemSize>>>(
                vector_td<unsigned int, D>(this->plan_.matrix_size_os_),
                vector_td<unsigned int, D>(this->plan_.matrix_padding_),
                this->plan_.num_samples_,
                num_coils,
                raw_pointer_cast(&this->plan_.trajectory_[0]),
                image_permuted.get_data_ptr(),
                samples.get_data_ptr() + repetition * this->plan_.num_samples_ *
                    this->plan_.num_frames_ * domain_size_coils,
                this->plan_.warp_size_power_,
                accumulate,
                this->plan_.d_kernel_);

            CHECK_FOR_CUDA_ERROR();
        }
    }


    template<class T, unsigned int D, template<class, unsigned int> class K>
    void ConvolverNC2C<T, D, K, ConvolutionType::STANDARD>::prepare(
        const thrust::device_vector<vector_td<REAL, D>>& trajectory)
    {
        // Allocate storage for and compute temporary prefix-sum variable
        // (#cells influenced per sample).
        thrust::device_vector<unsigned int> c_p_s(trajectory.size());
        thrust::device_vector<unsigned int> c_p_s_ps(trajectory.size());
        CHECK_FOR_CUDA_ERROR();

        REAL radius = this->plan_.kernel_.get_radius();
        transform(trajectory.begin(), trajectory.end(),
                  c_p_s.begin(), compute_num_cells_per_sample<REAL, D>(radius));
        inclusive_scan(c_p_s.begin(), c_p_s.end(), c_p_s_ps.begin(),
                       thrust::plus<unsigned int>()); // Prefix sum.

        // Build the vector of (grid_idx, sample_idx) tuples. Actually kept in
        // two separate vectors.
        unsigned int num_pairs = c_p_s_ps.back();
        c_p_s.clear();

        tuples_first = thrust::device_vector<unsigned int>(num_pairs);
        tuples_last = thrust::device_vector<unsigned int>(num_pairs);

        CHECK_FOR_CUDA_ERROR();

        // Fill tuple vector.
        write_pairs<REAL, D>(vector_td<unsigned int, D>(this->plan_.matrix_size_os_),
                             vector_td<unsigned int, D>(this->plan_.matrix_padding_),
                             this->plan_.num_samples_, this->plan_.num_frames_,
                             this->plan_.kernel_.get_width(),
                             raw_pointer_cast(&trajectory[0]),
                             raw_pointer_cast(&c_p_s_ps[0]),
                             raw_pointer_cast(&tuples_first[0]),
                             raw_pointer_cast(&tuples_last[0]));
        c_p_s_ps.clear();

        // Sort by grid indices.
        sort_by_key(tuples_first.begin(), tuples_first.end(), tuples_last.begin());

        // Each bucket_begin[i] indexes the first element of bucket i's list of points.
        // Each bucket_end[i] indexes one past the last element of bucket i's list of points.
        bucket_begin = thrust::device_vector<unsigned int>(
            this->plan_.num_frames_ * prod(this->plan_.matrix_size_os_ + this->plan_.matrix_padding_));
        bucket_end = thrust::device_vector<unsigned int>(
            this->plan_.num_frames_ * prod(this->plan_.matrix_size_os_ + this->plan_.matrix_padding_));

        CHECK_FOR_CUDA_ERROR();

        // Find the beginning of each bucket's list of points.
        thrust::counting_iterator<unsigned int> search_begin(0);
        thrust::lower_bound(tuples_first.begin(), tuples_first.end(),
                            search_begin, search_begin + this->plan_.num_frames_ *
                            prod(this->plan_.matrix_size_os_ + this->plan_.matrix_padding_),
                            bucket_begin.begin());

        // Find the end of each bucket's list of points.
        thrust::upper_bound(tuples_first.begin(), tuples_first.end(),
                            search_begin, search_begin + this->plan_.num_frames_ *
                            prod(this->plan_.matrix_size_os_ + this->plan_.matrix_padding_),
                            bucket_end.begin());
    }


    template<class T, unsigned int D, template<class, unsigned int> class K>
    void ConvolverNC2C<T, D, K, ConvolutionType::STANDARD>::compute(
        const cuNDArray<T>& samples,
        cuNDArray<T>& image,
        bool accumulate)
    {
        // Check if warp_size is a power of two. We do some modulus tricks in
        // the kernels that depend on this.
        if (!((cudaDeviceManager::Instance()->warp_size(this->plan_.device_) &
            (cudaDeviceManager::Instance()->warp_size(this->plan_.device_) - 1)) == 0))
        {
            throw cuda_error("cuGriddingConvolution: Unsupported hardware "
                             "(warp size is not a power of 2).");
        }

        // Compute number of batches.
        unsigned int num_batches = 1;
        for (unsigned int d = D; d < image.get_number_of_dimensions(); d++)
            num_batches *= image.get_size(d);
        num_batches /= this->plan_.num_frames_;

        // Set up grid and threads. We can convolve only max_coils batches per
        // run due to shared memory issues.
        unsigned int threads_per_block = CUDA_CONV_THREADS_PER_KERNEL;
        unsigned int max_coils = CUDA_CONV_MAX_COILS;
        unsigned int domain_size_coils_desired = num_batches;
        unsigned int num_repetitions = domain_size_coils_desired / max_coils +
            (((domain_size_coils_desired % max_coils) == 0) ? 0 : 1);
        unsigned int domain_size_coils = (num_repetitions == 1) ?
            domain_size_coils_desired : max_coils;
        unsigned int domain_size_coils_tail = (num_repetitions == 1) ?
            domain_size_coils_desired : 
            domain_size_coils_desired - (num_repetitions - 1) * domain_size_coils;

        // Block and grid dimensions.
        dim3 dimBlock(threads_per_block);
        dim3 dimGrid((prod(this->plan_.matrix_size_os_ + this->plan_.matrix_padding_) +
                      dimBlock.x - 1) / dimBlock.x, this->plan_.num_frames_);

        // Calculate how much shared memory to use per thread.
        size_t bytes_per_thread = domain_size_coils * sizeof(T);
        size_t bytes_per_thread_tail = domain_size_coils_tail * sizeof(T);

        // Define temporary image that includes padding.
        std::vector<size_t> padded_image_dims = to_std_vector(
            this->plan_.matrix_size_os_ + this->plan_.matrix_padding_);
        if (this->plan_.num_frames_ > 1)
            padded_image_dims.push_back(this->plan_.num_frames_);
        if (num_batches > 1)
            padded_image_dims.push_back(num_batches);
        cuNDArray<T> padded_image(padded_image_dims);
        
        // Prioritise shared memory over L1 cache.
        cudaFuncSetCacheConfig(NFFT_H_convolve_kernel<T, D, K>,
                               cudaFuncCachePreferShared);
        
        // Samples view dimensions.
        std::vector<size_t> view_dims = {
            this->plan_.num_samples_,
            this->plan_.num_frames_,
            0}; // Placeholder for num_coils.

        for (unsigned int repetition = 0; repetition < num_repetitions; repetition++)
        {
            // Number of coils in this repetition.
            size_t num_coils = (repetition == num_repetitions - 1) ?
                domain_size_coils_tail : domain_size_coils;
            
            // Samples view for current repetition.
            view_dims.back() = num_coils;
            cuNDArray<T> samples_view(view_dims,
                const_cast<T*>(samples.get_data_ptr()) +
                repetition * this->plan_.num_samples_ *
                this->plan_.num_frames_ * domain_size_coils);
            std::vector<size_t> permute_order = {2, 0, 1};
            auto samples_permuted = permute(samples_view, permute_order);
            
            // Size of shared memory.
            size_t sharedMemSize = (repetition == num_repetitions - 1) ?
                dimBlock.x * bytes_per_thread_tail :
                dimBlock.x * bytes_per_thread;

            // Launch CUDA kernel.
            NFFT_H_convolve_kernel<T, D, K>
                <<<dimGrid, dimBlock, sharedMemSize>>>(
                vector_td<unsigned int, D>(this->plan_.matrix_size_os_ +
                                           this->plan_.matrix_padding_),
                this->plan_.num_samples_,
                num_coils,
                raw_pointer_cast(&this->plan_.trajectory_[0]),
                padded_image.get_data_ptr() + repetition *
                    prod(this->plan_.matrix_size_os_ +
                         this->plan_.matrix_padding_) *
                    this->plan_.num_frames_ *
                    domain_size_coils,
                samples_permuted.get_data_ptr(),
                raw_pointer_cast(&tuples_last[0]),
                raw_pointer_cast(&bucket_begin[0]),
                raw_pointer_cast(&bucket_end[0]),
                this->plan_.warp_size_power_,
                this->plan_.d_kernel_);
        }

        CHECK_FOR_CUDA_ERROR();

        this->wrap_image(padded_image, image, accumulate);
    }


    template<class T, unsigned int D, template<class, unsigned int> class K>
    void ConvolverNC2C<T, D, K, ConvolutionType::STANDARD>::wrap_image(
        cuNDArray<T>& source,
        cuNDArray<T>& target,
        bool accumulate)
    {
        // Compute number of batches.
        unsigned int num_batches = 1;
        for (unsigned int d = D; d < source.get_number_of_dimensions(); d++)
            num_batches *= source.get_size(d);
        num_batches /= this->plan_.num_frames_;

        // Set dimensions of grid/blocks.
        unsigned int bdim = 256;
        dim3 dimBlock(bdim);
        dim3 dimGrid(prod(this->plan_.matrix_size_os_) / bdim,
                          this->plan_.num_frames_ * num_batches);

        // Safety check.
        if ((prod(this->plan_.matrix_size_os_) % bdim) != 0)
        {
            std::stringstream ss;
            ss << "Error: cuNFFT : the number of oversampled image elements must be a multiplum of the block size: "
            << bdim;
            throw std::runtime_error(ss.str());
        }

        // Invoke kernel.
        wrap_image_kernel<T, D>
            <<<dimGrid, dimBlock>>>(
            source.get_data_ptr(),
            target.get_data_ptr(),
            vector_td<unsigned int, D>(this->plan_.matrix_size_os_),
            vector_td<unsigned int, D>(this->plan_.matrix_padding_),
            accumulate);

        CHECK_FOR_CUDA_ERROR();
    }


    template<class T, unsigned int D, template<class, unsigned int> class K>
    void ConvolverNC2C<T, D, K, ConvolutionType::ATOMIC>::prepare(
        const thrust::device_vector<vector_td<REAL, D>>& trajectory)
    {
        // No need to do anything here.
        // Defined for completeness.
    }


    template<class T, unsigned int D, template<class, unsigned int> class K>
    void ConvolverNC2C<T, D, K, ConvolutionType::ATOMIC>::compute(
        const cuNDArray<T>& samples,
        cuNDArray<T>& image,
        bool accumulate)
    {
        // We need a device with compute capability >= 2.0 to use atomics.
        if (cudaDeviceManager::Instance()->major_version(this->plan_.device_) == 1)
        {
            throw cuda_error("cuGriddingConvolution: Atomic gridding "
                             "convolution supported only on devices with "
                             "compute capability 2.0 or higher.");
        }

        // We require the warp size to be a power of 2.
        if (!((cudaDeviceManager::Instance()->warp_size(this->plan_.device_) &
            (cudaDeviceManager::Instance()->warp_size(this->plan_.device_) - 1)) == 0))
        {
            throw cuda_error("cuGriddingConvolution: Unsupported hardware "
                             "(warp size is not a power of 2).");
        }
        
        // Compute number of batches.
        unsigned int num_batches = 1;
        for (unsigned int d = D; d < image.get_number_of_dimensions(); d++)
            num_batches *= image.get_size(d);
        num_batches /= this->plan_.num_frames_;

        // Set up grid and threads. We can convolve only max_coils batches per
        // run due to shared memory issues.
        unsigned int threads_per_block = CUDA_CONV_THREADS_PER_KERNEL;
        unsigned int max_coils = CUDA_CONV_MAX_COILS;
        unsigned int domain_size_coils_desired = num_batches;
        unsigned int num_repetitions = domain_size_coils_desired / max_coils +
            (((domain_size_coils_desired % max_coils) == 0) ? 0 : 1);
        unsigned int domain_size_coils =
            (num_repetitions == 1) ? domain_size_coils_desired : max_coils;
        unsigned int domain_size_coils_tail =
            (num_repetitions == 1) ? domain_size_coils_desired :
            domain_size_coils_desired - (num_repetitions - 1) * domain_size_coils;

        // Block and grid dimensions.
        dim3 dimBlock(threads_per_block);
        dim3 dimGrid((this->plan_.num_samples_ + dimBlock.x - 1) / dimBlock.x,
                     this->plan_.num_frames_);

        // Calculate how much shared memory to use per thread.
        size_t bytes_per_thread =
            domain_size_coils * sizeof(vector_td<REAL, D>);
        size_t bytes_per_thread_tail =
            domain_size_coils_tail * sizeof(vector_td<REAL, D>);

        // Clear image if not accumulating.
        if (!accumulate)
            clear(&image);

        for (unsigned int repetition = 0; repetition < num_repetitions; repetition++)
        {
            // Number of coils in this repetition.
            size_t num_coils = (repetition == num_repetitions - 1) ?
                domain_size_coils_tail : domain_size_coils;

            // Size of shared memory.
            size_t sharedMemSize = (repetition == num_repetitions - 1) ?
                dimBlock.x * bytes_per_thread_tail :
                dimBlock.x * bytes_per_thread;

            // Launch CUDA kernel.
            NFFT_H_atomic_convolve_kernel<T, D, K>
                <<<dimGrid, dimBlock, sharedMemSize>>>(
                vector_td<unsigned int, D>(this->plan_.matrix_size_os_),
                vector_td<unsigned int, D>(this->plan_.matrix_padding_),
                this->plan_.num_samples_,
                num_coils,
                raw_pointer_cast(&this->plan_.trajectory_[0]),
                samples.get_data_ptr() + repetition * this->plan_.num_samples_ *
                    this->plan_.num_frames_ * domain_size_coils,
                image.get_data_ptr() + repetition * prod(this->plan_.matrix_size_os_) *
                    this->plan_.num_frames_ * domain_size_coils,
                this->plan_.d_kernel_);
        }

        CHECK_FOR_CUDA_ERROR();
    }


    template<class T, unsigned int D, template<class, unsigned int> class K>
    void ConvolverNC2C<T, D, K, ConvolutionType::SPARSE_MATRIX>::prepare(
        const thrust::device_vector<vector_td<REAL, D>>& trajectory)
    {
        this->conv_matrix_ = std::make_unique<cuCsrMatrix<T>>(
            make_conv_matrix<T, D, K>(
                trajectory, this->plan_.matrix_size_os_, this->plan_.d_kernel_));
    }


    template<class T, unsigned int D, template<class, unsigned int> class K>
    void ConvolverNC2C<T, D, K, ConvolutionType::SPARSE_MATRIX>::compute(
        const cuNDArray<T>& samples,
        cuNDArray<T>& image,
        bool accumulate)
    { 
        unsigned int num_batches = 1;
        for (unsigned int d = D; d < image.get_number_of_dimensions(); d++)
            num_batches *= image.get_size(d);
        num_batches /= this->plan_.num_frames_;

        std::vector<size_t> sample_dims =
            { this->plan_.num_samples_ * this->plan_.num_frames_, num_batches};
        std::vector<size_t> image_dims =
            { prod(this->plan_.matrix_size_os_), num_batches };
        image_dims.push_back(this->plan_.num_frames_);

        cuNDArray<T> image_view(
            image_dims, image.get_data_ptr());
        cuNDArray<T> samples_view(
            sample_dims, const_cast<T*>(samples.get_data_ptr()));

        sparseMM(T(1.0), T(1.0), *this->conv_matrix_, samples_view, image_view, true);
    }

} // namespace Gadgetron

template class Gadgetron::cuGriddingConvolution<float, 1, Gadgetron::KaiserKernel>;
template class Gadgetron::cuGriddingConvolution<float, 2, Gadgetron::KaiserKernel>;
template class Gadgetron::cuGriddingConvolution<float, 3, Gadgetron::KaiserKernel>;
template class Gadgetron::cuGriddingConvolution<float, 4, Gadgetron::KaiserKernel>;

template class Gadgetron::cuGriddingConvolution<double, 1, Gadgetron::KaiserKernel>;
template class Gadgetron::cuGriddingConvolution<double, 2, Gadgetron::KaiserKernel>;
template class Gadgetron::cuGriddingConvolution<double, 3, Gadgetron::KaiserKernel>;
template class Gadgetron::cuGriddingConvolution<double, 4, Gadgetron::KaiserKernel>;

template class Gadgetron::cuGriddingConvolution<Gadgetron::complext<float>, 1, Gadgetron::KaiserKernel>;
template class Gadgetron::cuGriddingConvolution<Gadgetron::complext<float>, 2, Gadgetron::KaiserKernel>;
template class Gadgetron::cuGriddingConvolution<Gadgetron::complext<float>, 3, Gadgetron::KaiserKernel>;
template class Gadgetron::cuGriddingConvolution<Gadgetron::complext<float>, 4, Gadgetron::KaiserKernel>;

template class Gadgetron::cuGriddingConvolution<Gadgetron::complext<double>, 1, Gadgetron::KaiserKernel>;
template class Gadgetron::cuGriddingConvolution<Gadgetron::complext<double>, 2, Gadgetron::KaiserKernel>;
template class Gadgetron::cuGriddingConvolution<Gadgetron::complext<double>, 3, Gadgetron::KaiserKernel>;
template class Gadgetron::cuGriddingConvolution<Gadgetron::complext<double>, 4, Gadgetron::KaiserKernel>;

template class Gadgetron::cuGriddingConvolution<float, 1, Gadgetron::JincKernel>;
template class Gadgetron::cuGriddingConvolution<float, 2, Gadgetron::JincKernel>;
template class Gadgetron::cuGriddingConvolution<float, 3, Gadgetron::JincKernel>;
template class Gadgetron::cuGriddingConvolution<float, 4, Gadgetron::JincKernel>;

template class Gadgetron::cuGriddingConvolution<double, 1, Gadgetron::JincKernel>;
template class Gadgetron::cuGriddingConvolution<double, 2, Gadgetron::JincKernel>;
template class Gadgetron::cuGriddingConvolution<double, 3, Gadgetron::JincKernel>;
template class Gadgetron::cuGriddingConvolution<double, 4, Gadgetron::JincKernel>;

template class Gadgetron::cuGriddingConvolution<Gadgetron::complext<float>, 1, Gadgetron::JincKernel>;
template class Gadgetron::cuGriddingConvolution<Gadgetron::complext<float>, 2, Gadgetron::JincKernel>;
template class Gadgetron::cuGriddingConvolution<Gadgetron::complext<float>, 3, Gadgetron::JincKernel>;
template class Gadgetron::cuGriddingConvolution<Gadgetron::complext<float>, 4, Gadgetron::JincKernel>;

template class Gadgetron::cuGriddingConvolution<Gadgetron::complext<double>, 1, Gadgetron::JincKernel>;
template class Gadgetron::cuGriddingConvolution<Gadgetron::complext<double>, 2, Gadgetron::JincKernel>;
template class Gadgetron::cuGriddingConvolution<Gadgetron::complext<double>, 3, Gadgetron::JincKernel>;
template class Gadgetron::cuGriddingConvolution<Gadgetron::complext<double>, 4, Gadgetron::JincKernel>;

template class Gadgetron::GriddingConvolution<Gadgetron::cuNDArray, float, 1, Gadgetron::KaiserKernel>;
template class Gadgetron::GriddingConvolution<Gadgetron::cuNDArray, float, 2, Gadgetron::KaiserKernel>;
template class Gadgetron::GriddingConvolution<Gadgetron::cuNDArray, float, 3, Gadgetron::KaiserKernel>;
template class Gadgetron::GriddingConvolution<Gadgetron::cuNDArray, float, 4, Gadgetron::KaiserKernel>;

template class Gadgetron::GriddingConvolution<Gadgetron::cuNDArray, double, 1, Gadgetron::KaiserKernel>;
template class Gadgetron::GriddingConvolution<Gadgetron::cuNDArray, double, 2, Gadgetron::KaiserKernel>;
template class Gadgetron::GriddingConvolution<Gadgetron::cuNDArray, double, 3, Gadgetron::KaiserKernel>;
template class Gadgetron::GriddingConvolution<Gadgetron::cuNDArray, double, 4, Gadgetron::KaiserKernel>;

template class Gadgetron::GriddingConvolution<Gadgetron::cuNDArray, Gadgetron::complext<float>, 1, Gadgetron::KaiserKernel>;
template class Gadgetron::GriddingConvolution<Gadgetron::cuNDArray, Gadgetron::complext<float>, 2, Gadgetron::KaiserKernel>;
template class Gadgetron::GriddingConvolution<Gadgetron::cuNDArray, Gadgetron::complext<float>, 3, Gadgetron::KaiserKernel>;
template class Gadgetron::GriddingConvolution<Gadgetron::cuNDArray, Gadgetron::complext<float>, 4, Gadgetron::KaiserKernel>;

template class Gadgetron::GriddingConvolution<Gadgetron::cuNDArray, Gadgetron::complext<double>, 1, Gadgetron::KaiserKernel>;
template class Gadgetron::GriddingConvolution<Gadgetron::cuNDArray, Gadgetron::complext<double>, 2, Gadgetron::KaiserKernel>;
template class Gadgetron::GriddingConvolution<Gadgetron::cuNDArray, Gadgetron::complext<double>, 3, Gadgetron::KaiserKernel>;
template class Gadgetron::GriddingConvolution<Gadgetron::cuNDArray, Gadgetron::complext<double>, 4, Gadgetron::KaiserKernel>;

template class Gadgetron::GriddingConvolution<Gadgetron::cuNDArray, float, 1, Gadgetron::JincKernel>;
template class Gadgetron::GriddingConvolution<Gadgetron::cuNDArray, float, 2, Gadgetron::JincKernel>;
template class Gadgetron::GriddingConvolution<Gadgetron::cuNDArray, float, 3, Gadgetron::JincKernel>;
template class Gadgetron::GriddingConvolution<Gadgetron::cuNDArray, float, 4, Gadgetron::JincKernel>;

template class Gadgetron::GriddingConvolution<Gadgetron::cuNDArray, double, 1, Gadgetron::JincKernel>;
template class Gadgetron::GriddingConvolution<Gadgetron::cuNDArray, double, 2, Gadgetron::JincKernel>;
template class Gadgetron::GriddingConvolution<Gadgetron::cuNDArray, double, 3, Gadgetron::JincKernel>;
template class Gadgetron::GriddingConvolution<Gadgetron::cuNDArray, double, 4, Gadgetron::JincKernel>;

template class Gadgetron::GriddingConvolution<Gadgetron::cuNDArray, Gadgetron::complext<float>, 1, Gadgetron::JincKernel>;
template class Gadgetron::GriddingConvolution<Gadgetron::cuNDArray, Gadgetron::complext<float>, 2, Gadgetron::JincKernel>;
template class Gadgetron::GriddingConvolution<Gadgetron::cuNDArray, Gadgetron::complext<float>, 3, Gadgetron::JincKernel>;
template class Gadgetron::GriddingConvolution<Gadgetron::cuNDArray, Gadgetron::complext<float>, 4, Gadgetron::JincKernel>;

template class Gadgetron::GriddingConvolution<Gadgetron::cuNDArray, Gadgetron::complext<double>, 1, Gadgetron::JincKernel>;
template class Gadgetron::GriddingConvolution<Gadgetron::cuNDArray, Gadgetron::complext<double>, 2, Gadgetron::JincKernel>;
template class Gadgetron::GriddingConvolution<Gadgetron::cuNDArray, Gadgetron::complext<double>, 3, Gadgetron::JincKernel>;
template class Gadgetron::GriddingConvolution<Gadgetron::cuNDArray, Gadgetron::complext<double>, 4, Gadgetron::JincKernel>;
