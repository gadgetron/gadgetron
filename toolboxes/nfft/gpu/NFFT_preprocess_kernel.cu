//
// NFFT_H preprocessing kernels
//

// convert input trajectory in [-1/2;1/2] to [0;matrix_size_os+matrix_size_wrap]

template<class REAL, unsigned int D> struct trajectory_scale
{
  typename reald<REAL,D>::Type matrix, bias;
  
  trajectory_scale( const typename reald<REAL,D>::Type &m, const typename reald<REAL,D>::Type &b ){
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
