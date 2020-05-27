namespace Gadgetron
{
    namespace
    {        
        template<class REAL> __inline__ __device__
        void output(unsigned int num_samples, unsigned int num_batches, REAL* __restrict__ samples,
            unsigned int double_warp_size_power, unsigned int globalThreadId, unsigned int sharedMemFirstSampleIdx)
        {
            REAL *shared_mem = (REAL*)_shared_mem;
            
            for (unsigned int batch = 0; batch < num_batches; batch++)
            {
                REAL sample_value;
                sample_value = shared_mem[sharedMemFirstSampleIdx+(batch<<double_warp_size_power)];

                unsigned int out_idx = (batch*gridDim.y+blockIdx.y)*num_samples + globalThreadId;

                samples[out_idx] = sample_value;
            }
        }

        template<unsigned int D> __inline__ __device__
        static void resolve_wrap(vector_td<int, D> &grid_position, vector_td<unsigned int, D> &grid_size)
        {
            vector_td<int,D> zero(0);
            grid_position += vector_less(grid_position, zero) * grid_size;
            grid_position -= vector_greater_equal(grid_position, grid_size) * grid_size;
        }

        template<class REAL, unsigned int D> __inline__ __device__
        void iterate_body(REAL kernel_radius,
            vector_td<unsigned int, D> grid_size, unsigned int num_batches, const REAL* __restrict__ image,
            unsigned int double_warp_size_power, unsigned int sharedMemFirstSampleIdx,
            vector_td<REAL,D> sample_position, vector_td<int,D> grid_position)
        {

            using namespace thrust::cuda_cub;
            
            // calculate the distance between current sample and the grid cell
            vector_td<REAL,D> grid_position_real = vector_td<REAL,D>(grid_position);
            const vector_td<REAL,D> delta = abs(sample_position - grid_position_real);
            const vector_td<REAL,D> radius(kernel_radius);
            
            // if cell too distant from sample then move on to the next cell
            if (weak_greater(delta, radius))
                return;

            // compute convolution weight
            const REAL weight = cuSDC_kernel<REAL, D>::compute(delta, kernel_radius);

            // safety measure. We have occationally observed a NaN from the KaiserBessel computation
            if (!isfinite(weight))
                return;

            // resolve wrapping of grid position
            resolve_wrap<D>(grid_position, grid_size);

            REAL *shared_mem = (REAL*) _shared_mem;

            unsigned int image_idx = ((blockIdx.y) * prod(grid_size) + co_to_idx<D>( vector_td<unsigned int, D>(grid_position), grid_size ))*num_batches;

            for( unsigned int batch=0; batch < num_batches; batch++)
            {  
                // read the grid cell value from global memory
                const REAL grid_value = cub::ThreadLoad<cub::LOAD_LDG>(image + image_idx + batch);
                
                // add 'weight*grid_value' to the samples in shared memory
                shared_mem[sharedMemFirstSampleIdx + (batch << double_warp_size_power)] += (weight * grid_value);
                // shared_mem[sharedMemFirstSampleIdx + (batch << double_warp_size_power) + warpSize] += (weight * grid_value._imag);
            }
        }

        // this method is deliberately overloaded in 'UINTd' (rather than templatized) to improve performance of the loop iteration
        template<class REAL> __inline__ __device__
        void iterate(REAL kernel_radius, vector_td<unsigned int, 1> grid_size,
                unsigned int num_batches, const REAL* __restrict__ image,
                unsigned int double_warp_size_power, unsigned int sharedMemFirstSampleIdx,
                vector_td<REAL, 1> sample_position, vector_td<int,1> lower_limit, vector_td<int,1> upper_limit )
        {
            // iterate through all grid cells influencing the corresponding sample
            for (int x = lower_limit.vec[0]; x <= upper_limit.vec[0]; x++)
            {
                const intd<1>::Type grid_position(x);
                
                iterate_body<REAL,1>(kernel_radius, grid_size, num_batches, image, double_warp_size_power, 
                        sharedMemFirstSampleIdx, sample_position, grid_position);
            }
        }

        // this method is deliberately overloaded in 'UINTd' (rather than templatized) to improve performance of the loop iteration
        template<class REAL> __inline__ __device__
        void iterate(REAL kernel_radius, vector_td<unsigned int, 2> grid_size,
            unsigned int num_batches, const REAL* __restrict__ image,
                unsigned int double_warp_size_power, unsigned int sharedMemFirstSampleIdx,
                vector_td<REAL, 2> sample_position, vector_td<int, 2> lower_limit, vector_td<int, 2> upper_limit)
        {
            // iterate through all grid cells influencing the corresponding sample
            for (int y = lower_limit.vec[1]; y <= upper_limit.vec[1]; y++)
            {
                for (int x = lower_limit.vec[0]; x <= upper_limit.vec[0]; x++)
                {
                    const intd<2>::Type grid_position(x, y);
                    
                    iterate_body<REAL,2>(kernel_radius, grid_size, num_batches, image, double_warp_size_power, 
                        sharedMemFirstSampleIdx, sample_position, grid_position);
                }
            }
        }

        // this method is deliberately overloaded in 'd' (rather than templatized) to improve performance of the loop iteration
        template<class REAL> __inline__ __device__
        void iterate(REAL kernel_radius, vector_td<unsigned int, 3> grid_size,
            unsigned int num_batches, const REAL* __restrict__ image,
                unsigned int double_warp_size_power, unsigned int sharedMemFirstSampleIdx,
                vector_td<REAL, 3> sample_position, vector_td<int,3> lower_limit, vector_td<int, 3> upper_limit)
        {
            // iterate through all grid cells influencing the corresponding sample
            for (int z = lower_limit.vec[2]; z <= upper_limit.vec[2]; z++)
            {
                for (int y = lower_limit.vec[1]; y <= upper_limit.vec[1]; y++)
                {
                    for (int x = lower_limit.vec[0]; x <= upper_limit.vec[0]; x++)
                    {
                        const intd<3>::Type grid_position(x, y, z);
                        
                        iterate_body<REAL,3>(kernel_radius, grid_size, num_batches, image, double_warp_size_power, 
                                sharedMemFirstSampleIdx, sample_position, grid_position);
                    }
                }
            }
        }

        // this method is deliberately overloaded in 'd' (rather than templatized) to improve performance of the loop iteration
        template<class REAL> __inline__ __device__
        void iterate(REAL kernel_radius, vector_td<unsigned int, 4> grid_size,
            unsigned int num_batches, const REAL* __restrict__ image,
                unsigned int double_warp_size_power, unsigned int sharedMemFirstSampleIdx,
                vector_td<REAL, 4> sample_position, vector_td<int,4> lower_limit, vector_td<int, 4> upper_limit)
        {
            // iterate through all grid cells influencing the corresponding sample
            for (int w = lower_limit.vec[3]; w<=upper_limit.vec[3]; w++)
            {
                for (int z = lower_limit.vec[2]; z<=upper_limit.vec[2]; z++)
                {
                    for (int y = lower_limit.vec[1]; y<=upper_limit.vec[1]; y++)
                    {
                        for (int x = lower_limit.vec[0]; x<=upper_limit.vec[0]; x++)
                        {
                            const intd<4>::Type grid_position(x,y,z,w);
                            
                            iterate_body<REAL,4>(kernel_radius, grid_size, num_batches, image, double_warp_size_power, 
                                    sharedMemFirstSampleIdx, sample_position, grid_position);
                        }
                    }
                }
            }
        }

        template<class REAL, unsigned int D> __inline__ __device__
        void convolve(REAL kernel_radius, vector_td<unsigned int, D> grid_size, vector_td<unsigned int, D> grid_padding, 
                unsigned int num_samples, unsigned int num_batches, const vector_td<REAL,D> * __restrict__ traj, const REAL* __restrict__ image,
                unsigned int double_warp_size_power,
                unsigned int globalThreadId, unsigned int sharedMemFirstSampleIdx )
        {
            // sample position to convolve onto
            // computed in preprocessing, which included a wrap zone. Remove this wrapping.
            const vector_td<REAL,D> half_padding = vector_td<REAL,D>(grid_padding >> 1);
            const vector_td<REAL,D> sample_position = traj[globalThreadId+blockIdx.y*num_samples] - half_padding;
            
            // half the kernel width
            const vector_td<REAL,D> radius(kernel_radius);
            
            // limits of the subgrid to consider
            const vector_td<int,D> lower_limit = vector_td<int,D>(ceil(sample_position - radius));
            const vector_td<int,D> upper_limit = vector_td<int,D>(floor(sample_position + radius));

            // accumulate contributions from the grid
            iterate<REAL>(kernel_radius, grid_size, num_batches, image, double_warp_size_power, 
                    sharedMemFirstSampleIdx, sample_position, lower_limit, upper_limit);
        }
    }   // namespace


    namespace SDC_internal
    {
        // kernel main
        template<class REAL, unsigned int D> __global__
        void convolve_C2NC_kernel(
            REAL kernel_radius, 
            vector_td<unsigned int, D> grid_size, vector_td<unsigned int, D> grid_padding,
            unsigned int num_samples, unsigned int num_batches, 
            const vector_td<REAL, D>* __restrict__ traj,
            const REAL* __restrict__ image,  REAL* __restrict__ samples,
            unsigned int double_warp_size_power)
        {
            // global thread number	
            const unsigned int globalThreadId = blockIdx.x * blockDim.x + threadIdx.x;

            // check if we are within bounds
            if (globalThreadId >= num_samples)
                return;
            
            // number of reals to compute/output per thread
            const unsigned int num_reals = num_batches;
            
            // all shared memory reals corresponding to domain 'threadIdx.x' are located in bank threadIdx.x%warp_size to limit bank conflicts
            const unsigned int scatterSharedMemStart = (threadIdx.x / warpSize) * warpSize;
            const unsigned int scatterSharedMemStartOffset = threadIdx.x&(warpSize-1); // a faster way of saying (threadIdx.x%warpSize) 
            const unsigned int sharedMemFirstSampleIdx = scatterSharedMemStart * num_reals + scatterSharedMemStartOffset;

            REAL *shared_mem = (REAL*) _shared_mem;
            const REAL zero = REAL(0);

            // initialize shared memory
            for (unsigned int i = 0; i < num_reals; i++)
                shared_mem[sharedMemFirstSampleIdx+warpSize*i] = zero;
            
            // compute NFFT using arbitrary sample trajectories
            convolve<REAL,D>(kernel_radius, grid_size, grid_padding, num_samples, num_batches, 
                traj, image, double_warp_size_power, 
                globalThreadId, sharedMemFirstSampleIdx);
            
            // output k-space image to global memory
            output<REAL>(num_samples, num_batches, samples, double_warp_size_power, globalThreadId, sharedMemFirstSampleIdx);
        }
    }   // namespace SDC_internal
}   // namespace Gadgetron
