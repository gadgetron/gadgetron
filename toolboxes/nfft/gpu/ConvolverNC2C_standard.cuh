#pragma once
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

#include "ConvolutionKernel.h"


//
// NFFT_H preprocessing kernels
//

// convert input trajectory in [-1/2;1/2] to [0;matrix_size_os+matrix_size_wrap]

template<class REAL, unsigned int D> struct trajectory_scale
{
    typename reald<REAL,D>::Type matrix, bias;

    trajectory_scale( const typename reald<REAL,D>::Type &m, const typename reald<REAL,D>::Type &b )
    {
        matrix = m;
        bias = b;
    }

    __host__ __device__
    typename reald<REAL,D>::Type operator()(const typename reald<REAL,D>::Type &in) const {
        return component_wise_mul<REAL,D>(in,matrix)+bias;
    }
};

template<class REAL, unsigned int D>
struct compute_num_cells_per_sample
{
    __host__ __device__
    compute_num_cells_per_sample(REAL _half_W) : half_W(_half_W) {}

    __host__ __device__
    unsigned int operator()(typename reald<REAL,D>::Type p) const
    {
        unsigned int num_cells = 1;
        for( unsigned int dim=0; dim<D; dim++ ){
            unsigned int upper_limit = (unsigned int)floor((((float*)&p)[dim])+half_W);
            unsigned int lower_limit = (unsigned int)ceil((((float*)&p)[dim])-half_W);
            num_cells *= (upper_limit-lower_limit+1);
        }
        return num_cells;
    }

    REAL half_W;
};

template<class REAL> __inline__ __device__ void
output_pairs( unsigned int sample_idx, unsigned int frame,
	      typename reald<REAL,1>::Type &p, typename uintd<1>::Type &matrix_size_os, typename uintd<1>::Type &matrix_size_wrap,
	      REAL half_W, const unsigned int * __restrict__ write_offsets, unsigned int * __restrict__ tuples_first, unsigned int * __restrict__ tuples_last )
{
    unsigned int lower_limit_x = (unsigned int)ceil(p.vec[0]-half_W);
    unsigned int upper_limit_x = (unsigned int)floor(p.vec[0]+half_W);

    unsigned int pair_idx = 0;
    unsigned int write_offset = (sample_idx==0) ? 0 : write_offsets[sample_idx-1];
    unsigned int frame_offset = frame*prod(matrix_size_os+matrix_size_wrap);
    for( unsigned int x=lower_limit_x; x<=upper_limit_x; x++ ){
        typename uintd<1>::Type co; co.vec[0] = x;
        tuples_first[write_offset+pair_idx] = co_to_idx(co, matrix_size_os+matrix_size_wrap)+frame_offset;
        tuples_last[write_offset+pair_idx] = sample_idx;
        pair_idx++;
    }
}

template<class REAL> __inline__ __device__ void
output_pairs( unsigned int sample_idx, unsigned int frame,
	      typename reald<REAL,2>::Type &p, typename uintd<2>::Type &matrix_size_os, typename uintd<2>::Type &matrix_size_wrap,
	      REAL half_W, const unsigned int * __restrict__ write_offsets, unsigned int * __restrict__ tuples_first, unsigned int * __restrict__ tuples_last )
{
    unsigned int lower_limit_x = (unsigned int)ceil(p.vec[0]-half_W);
    unsigned int lower_limit_y = (unsigned int)ceil(p.vec[1]-half_W);
    unsigned int upper_limit_x = (unsigned int)floor(p.vec[0]+half_W);
    unsigned int upper_limit_y = (unsigned int)floor(p.vec[1]+half_W);

    unsigned int pair_idx = 0;
    unsigned int write_offset = (sample_idx==0) ? 0 : write_offsets[sample_idx-1];
    unsigned int frame_offset = frame*prod(matrix_size_os+matrix_size_wrap);
    for( unsigned int y=lower_limit_y; y<=upper_limit_y; y++ ){
        for( unsigned int x=lower_limit_x; x<=upper_limit_x; x++ ){
            typename uintd<2>::Type co; co.vec[0] = x; co.vec[1] = y;
            tuples_first[write_offset+pair_idx] = co_to_idx(co, matrix_size_os+matrix_size_wrap)+frame_offset;
            tuples_last[write_offset+pair_idx] = sample_idx;
            pair_idx++;
        }
    }
}

template <class REAL> __inline__ __device__ void
output_pairs( unsigned int sample_idx, unsigned int frame,
	      typename reald<REAL,3>::Type &p, typename uintd<3>::Type &matrix_size_os, typename uintd<3>::Type &matrix_size_wrap,
	      REAL half_W, const unsigned int * __restrict__ write_offsets, unsigned int * __restrict__ tuples_first, unsigned int * __restrict__ tuples_last )
{
    unsigned int lower_limit_x = (unsigned int)ceil(p.vec[0]-half_W);
    unsigned int lower_limit_y = (unsigned int)ceil(p.vec[1]-half_W);
    unsigned int lower_limit_z = (unsigned int)ceil(p.vec[2]-half_W);
    unsigned int upper_limit_x = (unsigned int)floor(p.vec[0]+half_W);
    unsigned int upper_limit_y = (unsigned int)floor(p.vec[1]+half_W);
    unsigned int upper_limit_z = (unsigned int)floor(p.vec[2]+half_W);

    unsigned int pair_idx = 0;
    unsigned int write_offset = (sample_idx==0) ? 0 : write_offsets[sample_idx-1];
    unsigned int frame_offset = frame*prod(matrix_size_os+matrix_size_wrap);
    for( unsigned int z=lower_limit_z; z<=upper_limit_z; z++ ){
        for( unsigned int y=lower_limit_y; y<=upper_limit_y; y++ ){
            for( unsigned int x=lower_limit_x; x<=upper_limit_x; x++ ){
                typename uintd<3>::Type co; co.vec[0] = x; co.vec[1] = y; co.vec[2] = z;
                tuples_first[write_offset+pair_idx] = co_to_idx(co, matrix_size_os+matrix_size_wrap)+frame_offset;
                tuples_last[write_offset+pair_idx] = sample_idx;
                pair_idx++;
            }
        }
    }
}

template <class REAL> __inline__ __device__ void
output_pairs( unsigned int sample_idx, unsigned int frame,
	      typename reald<REAL,4>::Type &p, typename uintd<4>::Type &matrix_size_os, typename uintd<4>::Type &matrix_size_wrap,
	      REAL half_W, const unsigned int * __restrict__ write_offsets, unsigned int * __restrict__ tuples_first, unsigned int * __restrict__ tuples_last )
{
    unsigned int lower_limit_x = (unsigned int)ceil(p.vec[0]-half_W);
    unsigned int lower_limit_y = (unsigned int)ceil(p.vec[1]-half_W);
    unsigned int lower_limit_z = (unsigned int)ceil(p.vec[2]-half_W);
    unsigned int lower_limit_w = (unsigned int)ceil(p.vec[3]-half_W);
    unsigned int upper_limit_x = (unsigned int)floor(p.vec[0]+half_W);
    unsigned int upper_limit_y = (unsigned int)floor(p.vec[1]+half_W);
    unsigned int upper_limit_z = (unsigned int)floor(p.vec[2]+half_W);
    unsigned int upper_limit_w = (unsigned int)floor(p.vec[3]+half_W);

    unsigned int pair_idx = 0;
    unsigned int write_offset = (sample_idx==0) ? 0 : write_offsets[sample_idx-1];
    unsigned int frame_offset = frame*prod(matrix_size_os+matrix_size_wrap);
    for( unsigned int w=lower_limit_w; w<=upper_limit_w; w++ ){
        for( unsigned int z=lower_limit_z; z<=upper_limit_z; z++ ){
            for( unsigned int y=lower_limit_y; y<=upper_limit_y; y++ ){
	            for( unsigned int x=lower_limit_x; x<=upper_limit_x; x++ ){
                    typename uintd<4>::Type co; co.vec[0] = x; co.vec[1] = y; co.vec[2] = z; co.vec[3] = w;
                    tuples_first[write_offset+pair_idx] = co_to_idx(co, matrix_size_os+matrix_size_wrap)+frame_offset;
                    tuples_last[write_offset+pair_idx] = sample_idx;
                    pair_idx++;
	            }
            }
        }
    }
}

template<class REAL, unsigned int D> __global__ void
write_pairs_kernel( typename uintd<D>::Type matrix_size_os, typename uintd<D>::Type matrix_size_wrap, unsigned int num_samples_per_frame, REAL half_W,
		   const typename reald<REAL,D>::Type * __restrict__ traj_positions, unsigned int * __restrict__ write_offsets, unsigned int * __restrict__ tuples_first, unsigned int * __restrict__ tuples_last )
{
    // Get sample idx
    unsigned int sample_idx = blockIdx.x*blockDim.x + threadIdx.x;
    unsigned int frame = blockIdx.y;

    if( sample_idx<num_samples_per_frame ){

        sample_idx += frame*num_samples_per_frame;
        typename reald<REAL,D>::Type p = traj_positions[sample_idx];
        output_pairs<REAL>( sample_idx, frame, p, matrix_size_os, matrix_size_wrap, half_W, write_offsets, tuples_first, tuples_last );
    }
};

template <class REAL, unsigned int D> void
write_pairs( typename uintd<D>::Type matrix_size_os, typename uintd<D>::Type matrix_size_wrap, unsigned int num_samples_per_frame, unsigned int num_frames, REAL W,
	     const typename reald<REAL,D>::Type * __restrict__ traj_positions, unsigned int * __restrict__ write_offsets, unsigned int * __restrict__ tuples_first, unsigned int * __restrict__ tuples_last )
{
    dim3 blockDim(256);
    dim3 gridDim((int)ceil((double)num_samples_per_frame/(double)blockDim.x), num_frames);

    REAL half_W = REAL(0.5)*W;
    write_pairs_kernel<REAL,D><<< gridDim, blockDim >>>
        ( matrix_size_os, matrix_size_wrap, num_samples_per_frame, half_W, traj_positions, write_offsets, tuples_first, tuples_last );

    CHECK_FOR_CUDA_ERROR();
}


//
// There is no header file accompanying this kernel, so it makes most sense to read the code/file from the end and upwards
//

//
// Transfer result from shared memory to global memory.
//

template<class T>
__inline__ __device__
void NFFT_H_output(unsigned int number_of_batches, T* __restrict__ image,
	       unsigned int warp_size_power, unsigned int number_of_domains,
	       unsigned int globalThreadId, unsigned int sharedMemFirstCellIdx)
{
    realType_t<T> *shared_mem = (realType_t<T>*) _shared_mem;

    for (unsigned int batch = 0; batch < number_of_batches; batch++)
    {
        T cell_coefficient;
        if constexpr (is_complex_type_v<T>)
        {
            cell_coefficient._real =
                shared_mem[sharedMemFirstCellIdx + (batch << warp_size_power)];
            cell_coefficient._imag =
                shared_mem[sharedMemFirstCellIdx + (batch << warp_size_power) + warpSize];
        }
        else
        {
            cell_coefficient = shared_mem[sharedMemFirstCellIdx + (batch << warp_size_power)];
        }
        image[(batch*gridDim.y+blockIdx.y) * number_of_domains + globalThreadId] =
            cell_coefficient;
    }
}


template<class T, unsigned int D, template<class, unsigned int> class K>
__inline__ __device__
void NFFT_H_convolve(
    unsigned int number_of_samples, unsigned int number_of_batches,
    unsigned int number_of_domains,
    const vector_td<realType_t<T>,D> * __restrict__ traj_positions, const  T* __restrict__ samples, const unsigned int * __restrict__ tuples_last,
    const unsigned int * __restrict__ bucket_begin, const unsigned int * __restrict__ bucket_end,
    unsigned int warp_size_power,
    unsigned int globalThreadId, vector_td<unsigned int,D> domainPos,
    unsigned int sharedMemFirstCellIdx,
    const ConvolutionKernel<realType_t<T>, D, K>* kernel)
{
    using namespace thrust::cuda_cub;

    realType_t<T> *shared_mem = (realType_t<T>*) _shared_mem;

    // Cell position as reald
    vector_td<realType_t<T>,D> cell_pos = vector_td<realType_t<T>,D>( domainPos );

    // Convolve samples onto the domain (shared memory)
    const unsigned int frame_offset = blockIdx.y*number_of_domains;
    for (unsigned int i = bucket_begin[globalThreadId + frame_offset];
         i < bucket_end[globalThreadId + frame_offset];
         i++)
    {
        // Safety precaution. TODO
        unsigned int sampleIdx = tuples_last[i];

        // Safety precaution. TODO
        vector_td<realType_t<T>,D> sample_pos = traj_positions[sampleIdx];

        // Calculate the distance between the cell and the sample.
        vector_td<realType_t<T>, D> delta = abs(sample_pos-cell_pos);

        // Compute convolution weights.
        float weight = kernel->get(delta);

        // Safety measure.
        if (!isfinite(weight))
            continue;

        // Apply filter to input images.
        for (int batch = 0; batch < number_of_batches; batch++)
        {
            unsigned int idx = sampleIdx * number_of_batches + batch;

            // Apply filter to shared memory domain.
            T sample_val  = cub::ThreadLoad<cub::LOAD_CS>(samples + idx);

            if constexpr (is_complex_type_v<T>)
            {
                shared_mem[sharedMemFirstCellIdx + (batch << warp_size_power)] +=
                    (weight * sample_val._real);
                shared_mem[sharedMemFirstCellIdx + (batch << warp_size_power) + warpSize] +=
                    (weight * sample_val._imag);
            }
            else
            {
                shared_mem[sharedMemFirstCellIdx + (batch << warp_size_power)] +=
                    (weight * sample_val);
            }
        }
    }
}

//
// kernel main
//

template<class T, unsigned int D, template<class, unsigned int> class K>
__global__
void NFFT_H_convolve_kernel(vector_td<unsigned int,D> domain_count_grid,
    unsigned int number_of_samples, unsigned int number_of_batches,
    const vector_td<realType_t<T>,D> * __restrict__ traj_positions,
    T* __restrict__ image, const T* __restrict__ samples,
    const unsigned int * __restrict__ tuples_last,
    const unsigned int * __restrict__ bucket_begin,
    const unsigned int * __restrict__ bucket_end,
    unsigned int warp_size_power,
    const ConvolutionKernel<realType_t<T>, D, K>* kernel)
{

    // Global thread index.
    const unsigned int index = blockIdx.x*blockDim.x + threadIdx.x;

    // Number of domains.
    const unsigned int number_of_domains = prod(domain_count_grid);

    // Check if we are within bounds.
    if (index >= number_of_domains)
        return;

    // Mapped global thread index (actually we don't use a map currently).
    const unsigned int domainIdx = index;

    // Compute global domain position.
    const vector_td<unsigned int,D> domainPos = idx_to_co(domainIdx, domain_count_grid);

    // Number of cells.
    const unsigned int num_reals = is_complex_type_v<T> ?
        number_of_batches << 1 : number_of_batches;

    // For complex numbers, we need twice as many real samples per batch.
    if constexpr (is_complex_type_v<T>)
        warp_size_power += 1;

    // All shared memory floats corresponding to domain 'threadIdx.x' is located in bank threadIdx.x%warp_size to limit bank conflicts
    const unsigned int scatterSharedMemStart = (threadIdx.x/warpSize)*warpSize;
    const unsigned int scatterSharedMemStartOffset = threadIdx.x&(warpSize-1); // a faster way of saying (threadIdx.x%warpSize)
    const unsigned int sharedMemFirstCellIdx = scatterSharedMemStart*num_reals + scatterSharedMemStartOffset;

    realType_t<T> *shared_mem = (realType_t<T>*) _shared_mem;

    // Initialize shared memory.
    for (unsigned int i = 0; i < num_reals; i++)
        shared_mem[sharedMemFirstCellIdx+warpSize*i] = realType_t<T>(0);

    // Compute NFFT using arbitrary sample trajectories.
    NFFT_H_convolve<T, D>
      (number_of_samples, number_of_batches, number_of_domains,
        traj_positions, samples, tuples_last, bucket_begin, bucket_end,
        warp_size_power, index, domainPos,
        sharedMemFirstCellIdx, kernel);

    // Output k-space image to global memory.
    NFFT_H_output<T>( number_of_batches, image, warp_size_power, number_of_domains, index, sharedMemFirstCellIdx );
}


template<class T, unsigned int D>
__global__
void wrap_image_kernel(const T* __restrict__ in,
                        T* __restrict__ out,
                        vector_td<unsigned int, D> matrix_size_os,
                        vector_td<unsigned int, D> matrix_padding,
                        bool accumulate)
{
    unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
    const unsigned int num_elements_per_image_src = prod(matrix_size_os + matrix_padding);
    const unsigned int image_offset_src = blockIdx.y * num_elements_per_image_src;

    const typename uintd<D>::Type co = idx_to_co(idx, matrix_size_os);
    const typename uintd<D>::Type half_wrap = matrix_padding >> 1;

    // Make "boolean" vectors denoting whether wrapping needs to be
    // performed in a given direction (forwards/backwards).
    vector_td<bool, D>
    B_l = vector_less(co, half_wrap);
    vector_td<bool, D>
    B_r = vector_greater_equal(co, matrix_size_os - half_wrap);

    T result =
        in[co_to_idx(co + half_wrap, matrix_size_os + matrix_padding) + image_offset_src];

    if (sum(B_l + B_r) > 0)
    {
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

        for (unsigned char combination = 1; combination < (1 << D); combination++)
        {
            // Find d
            unsigned char d = 0;
            for (unsigned char i = 0; i < D; i++)
                d += ((combination & (1 << i)) > 0);

            // Create stride vector for each wrapping test.
            for (unsigned char s = 0; s < (1 << d); s++) {

                // Target for stride.
                typename intd<D>::Type stride;
                char wrap_requests = 0;
                char skipped_dims = 0;

                // Fill dimensions of the stride.
                for (unsigned char i = 1; i < D + 1; i++)
                {
                    // Is the stride dimension present in the current
                    // combination?
                    if (i & combination)
                    {
                        // A zero bit in s indicates "check for left wrap"
                        // and a one bit is interpreted as "check for right
                        // wrap" ("left/right" for the individual dimension
                        // meaning wrapping on either side of the dimension).
                        if (i & (s << (skipped_dims)))
                        {
                            if (B_r.vec[i - 1]) { // Wrapping required.
                                stride[i - 1] = -1;
                                wrap_requests++;
                            } else
                                stride[i - 1] = 0;
                        }
                        else
                        {
                            if (B_l.vec[i - 1]) { // Wrapping required.
                                stride[i - 1] = 1;
                                wrap_requests++;
                            } else
                                stride[i - 1] = 0;
                        }
                    }
                    else
                    {
                        // Do not test for wrapping in dimension 'i-1'
                        // (for this combination).
                        stride[i - 1] = 0;
                        skipped_dims++;
                    }
                }

                // Now it is time to do the actual wrapping (if needed).
                if (wrap_requests == d)
                {
                    typename intd<D>::Type src_co_int =
                        vector_td<int, D>(co + half_wrap);
                    typename intd<D>::Type matrix_size_os_int =
                        vector_td<int, D>(matrix_size_os);
                    typename intd<D>::Type co_offset_int =
                        src_co_int + component_wise_mul<int, D>(stride, matrix_size_os_int);
                    typename uintd<D>::Type co_offset =
                        vector_td<unsigned int, D>(co_offset_int);
                    result += in[co_to_idx(co_offset, matrix_size_os + matrix_padding) + image_offset_src];
                    break; // Only one stride per combination can contribute (e.g. one edge, one corner).
                }
            }
        }
    }

    // Output.
    const unsigned int image_offset_tgt = blockIdx.y * prod(matrix_size_os);
    if (accumulate)
        result += out[idx + image_offset_tgt];
    out[idx + image_offset_tgt] = result;
}
