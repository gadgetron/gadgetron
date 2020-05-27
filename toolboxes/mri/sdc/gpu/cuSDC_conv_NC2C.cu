/**
 * \file cuSDC_conv_NC2C.cu
 * \brief Convolution for sampling density compensation (CUDA specialization).
 */

namespace Gadgetron
{
    namespace
    {
        template<class REAL, unsigned int D> __inline__ __device__
        static void iterate_body(REAL kernel_radius, vector_td<unsigned int, D> grid_size, 
            unsigned int num_batches, const REAL* __restrict__ samples, REAL* __restrict__ image,
            unsigned int frame, unsigned int num_frames,
            unsigned int num_samples_per_batch, unsigned int sample_idx_in_batch, 
            vector_td<REAL,D> sample_position, vector_td<int,D> grid_position )
        {
            // calculate the distance between current sample and the grid cell
            vector_td<REAL, D> grid_position_real = vector_td<REAL,D>(grid_position);
            const vector_td<REAL, D> delta = abs(sample_position - grid_position_real);
            const vector_td<REAL, D> radius(kernel_radius);
            
            // if cell too distant from sample then move on to the next cell
            if (weak_greater(delta, radius))
                return;

            // compute convolution weight.
            const REAL weight = cuSDC_kernel<REAL, D>::compute(delta, kernel_radius);
            
            // safety measure. We have occationally observed a NaN from the KaiserBessel computation
            if (!isfinite(weight))
                return;

            // resolve wrapping of grid position
            resolve_wrap<D>(grid_position, grid_size);

            for (unsigned int batch = 0; batch < num_batches; batch++)
            {
                // read the grid sample value from global memory
                REAL sample_value = samples[sample_idx_in_batch+batch*num_samples_per_batch];
                
                // determine the grid cell idx
                unsigned int grid_idx = (batch*num_frames+frame)*prod(grid_size) + co_to_idx<D>( vector_td<unsigned int, D>(grid_position), grid_size );

                // atomic update of real and imaginary component
                atomicAdd(&(((REAL*)image)[grid_idx]), weight * sample_value);
            }
        }

        // this method is deliberately overloaded in 'UINTd' (rather than templetized) to improve performance of the loop iteration
        template<class REAL> __inline__ __device__
        void iterate(REAL kernel_radius, vector_td<unsigned int,1> grid_size, 
                unsigned int num_batches, const REAL* __restrict__ samples, REAL* __restrict__ image,
                unsigned int frame, unsigned int num_frames, 
                unsigned int num_samples_per_batch, unsigned int sample_idx_in_batch, 
                vector_td<REAL,1> sample_position, vector_td<int,1> lower_limit, vector_td<int, 1> upper_limit)
        {
            // iterate through all grid cells influencing the corresponding sample
            for (int x = lower_limit.vec[0]; x <= upper_limit.vec[0]; x++)
            {
                const intd<1>::Type grid_position(x);
                
                iterate_body<REAL,1>(kernel_radius, grid_size, num_batches, samples, image,
                        frame, num_frames,
                        num_samples_per_batch, sample_idx_in_batch, sample_position, grid_position);
            }
        }

        // this method is deliberately overloaded in 'd' (rather than templetized) to improve performance of the loop iteration
        template<class REAL> __inline__ __device__
        void iterate(REAL kernel_radius, vector_td<unsigned int,2> grid_size, 
                unsigned int num_batches, const REAL* __restrict__ samples, REAL* __restrict__ image,
                unsigned int frame, unsigned int num_frames, 
                unsigned int num_samples_per_batch, unsigned int sample_idx_in_batch, 
                vector_td<REAL,2> sample_position, vector_td<int,2> lower_limit, vector_td<int,2> upper_limit)
        {
            // iterate through all grid cells influencing the corresponding sample
            for (int y = lower_limit.vec[1]; y <= upper_limit.vec[1]; y++)
            {
                for (int x = lower_limit.vec[0]; x <= upper_limit.vec[0]; x++)
                {
                    const intd<2>::Type grid_position(x, y);
                    
                    iterate_body<REAL,2>(kernel_radius, grid_size, num_batches, samples, image,
                        frame, num_frames,
                        num_samples_per_batch, sample_idx_in_batch, sample_position, grid_position);
                }
            }
        }

        // this method is deliberately overloaded in 'd' (rather than templetized) to improve performance of the loop iteration
        template<class REAL> __inline__ __device__
        void iterate(REAL kernel_radius, vector_td<unsigned int,3> grid_size, 
                unsigned int num_batches, const REAL* __restrict__ samples, REAL* __restrict__ image,
                unsigned int frame, unsigned int num_frames, 	      
                unsigned int num_samples_per_batch, unsigned int sample_idx_in_batch, 
                vector_td<REAL, 3> sample_position, vector_td<int, 3> lower_limit, vector_td<int, 3> upper_limit)
        {
            // iterate through all grid cells influencing the corresponding sample
            for (int z = lower_limit.vec[2]; z <= upper_limit.vec[2]; z++)
            {
                for (int y = lower_limit.vec[1]; y <= upper_limit.vec[1]; y++)
                {
                    for (int x = lower_limit.vec[0]; x <= upper_limit.vec[0]; x++)
                    {
                        const intd<3>::Type grid_position(x, y, z);
                        
                        iterate_body<REAL, 3>(kernel_radius, grid_size, num_batches, samples, image, 
                                frame, num_frames,
                                num_samples_per_batch, sample_idx_in_batch, sample_position, grid_position);
                    }
                }
            }
        }

        // this method is deliberately overloaded in 'd' (rather than templetized) to improve performance of the loop iteration
        template<class REAL> __inline__ __device__
        void iterate(REAL kernel_radius, vector_td<unsigned int,4> grid_size, 
                unsigned int num_batches, const REAL* __restrict__ samples, REAL* __restrict image,
                unsigned int frame, unsigned int num_frames, 
                unsigned int num_samples_per_batch, unsigned int sample_idx_in_batch, 
                vector_td<REAL,4> sample_position, vector_td<int,4> lower_limit, vector_td<int,4> upper_limit)
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
                            
                            iterate_body<REAL, 4>(kernel_radius, grid_size, num_batches, samples, image,
                                    frame, num_frames,
                                    num_samples_per_batch, sample_idx_in_batch, sample_position, grid_position);
                        }
                    }
                }
            }
        }
    }   // namespace


    namespace SDC_internal
    {
        // kernel main
        template<class REAL, unsigned int D> __global__
        void convolve_NC2C_kernel(
            REAL kernel_radius, vector_td<unsigned int, D> grid_size, vector_td<unsigned int, D> grid_padding,
            unsigned int num_samples, unsigned int num_batches, 
            const vector_td<REAL,D>* __restrict__ traj,
            const REAL* __restrict__ samples, REAL* __restrict__ image)
        {
        
            // A runtime check will prevent this kernel from being run for compute models 1.x.            
            #if(__CUDA_ARCH__>=200)
            
            const unsigned int sample_idx_in_frame = blockIdx.x * blockDim.x + threadIdx.x;
        
            // check if we are within bounds
            if (sample_idx_in_frame >= num_samples)
                return;
                
            const unsigned int frame = blockIdx.y;
            const unsigned int num_frames = gridDim.y;
            const unsigned int num_samples_per_batch = num_samples * num_frames;
            const unsigned int sample_idx_in_batch = sample_idx_in_frame + frame * num_samples;
            
            // sample position computed in preprocessing includes a wrap zone. Remove this wrapping.
            const vector_td<REAL, D> half_padding = vector_td<REAL, D>(grid_padding >> 1);
            const vector_td<REAL, D> sample_position = traj[sample_idx_in_batch] - half_padding;
            
            // half the kernel width
            const vector_td<REAL, D> radius = vector_td<REAL, D>(kernel_radius);
            
            // limits of the subgrid to consider
            const vector_td<int, D> lower_limit = vector_td<int, D>(ceil(sample_position - radius));
            const vector_td<int, D> upper_limit = vector_td<int, D>(floor(sample_position + radius));
                    
            // output to the grid
            iterate<REAL>(kernel_radius, grid_size, num_batches, samples, image, 
                    frame, num_frames, num_samples_per_batch, sample_idx_in_batch, 
                    sample_position, lower_limit, upper_limit );

            #endif
        }
    }   // namespace SDC_internal
}   // namespace Gadgetron