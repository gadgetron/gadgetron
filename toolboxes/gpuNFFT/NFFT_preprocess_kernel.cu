//
// NFFT_H preprocessing kernels
//

// convert input trajectory in [-1/2;1/2] to [0;matrix_size_os+matrix_size_wrap]

template<class REAL, unsigned int D> struct trajectory_scale
{
  vector_td<REAL,D> matrix, bias;
  
  trajectory_scale( const vector_td<REAL,D> &m, const vector_td<REAL,D> &b ){
    matrix = m;
    bias = b;
  }
  
  __host__ __device__
  vector_td<REAL,D> operator()(const vector_td<REAL,D> &in) const { 
    return component_wise_mul<REAL,D>(in,matrix)+bias;
  }
};

template<class REAL, unsigned int D>
struct compute_num_cells_per_sample
{
  __host__ __device__
  compute_num_cells_per_sample(REAL _half_W) : half_W(_half_W) {}
  
  __host__ __device__
  unsigned int operator()(vector_td<REAL,D> p) const
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

// TODO: can we avoid overloading here?
template<class REAL> __inline__ __device__ void
output_pairs( unsigned int sample_idx, vector_td<REAL,2> &p, vector_td<unsigned int, 2> &matrix_size_os, vector_td<unsigned int, 2> &matrix_size_wrap, REAL half_W, unsigned int *write_offsets, unsigned int *tuples_first, unsigned int *tuples_last )
{
  unsigned int lower_limit_x = (unsigned int)ceil(p.vec[0]-half_W);
  unsigned int lower_limit_y = (unsigned int)ceil(p.vec[1]-half_W);
  unsigned int upper_limit_x = (unsigned int)floor(p.vec[0]+half_W);
  unsigned int upper_limit_y = (unsigned int)floor(p.vec[1]+half_W);

  unsigned int pair_idx = 0;
  unsigned int write_offset = (sample_idx==0) ? 0 : write_offsets[sample_idx-1];
  for( unsigned int y=lower_limit_y; y<=upper_limit_y; y++ ){
    for( unsigned int x=lower_limit_x; x<=upper_limit_x; x++ ){
      vector_td<unsigned int, 2> co; co.vec[0] = x; co.vec[1] = y;
      tuples_first[write_offset+pair_idx] = co_to_idx<2>(co, matrix_size_os+matrix_size_wrap);
      tuples_last[write_offset+pair_idx] = sample_idx;
      pair_idx++;
    }
  }
}

template <class REAL> __inline__ __device__ void
output_pairs( unsigned int sample_idx, vector_td<REAL,3> &p, vector_td<unsigned int, 3> &matrix_size_os, vector_td<unsigned int, 3> &matrix_size_wrap, REAL half_W, unsigned int *write_offsets, unsigned int *tuples_first, unsigned int *tuples_last )
{
  unsigned int lower_limit_x = (unsigned int)ceil(p.vec[0]-half_W);
  unsigned int lower_limit_y = (unsigned int)ceil(p.vec[1]-half_W);
  unsigned int lower_limit_z = (unsigned int)ceil(p.vec[2]-half_W);
  unsigned int upper_limit_x = (unsigned int)floor(p.vec[0]+half_W);
  unsigned int upper_limit_y = (unsigned int)floor(p.vec[1]+half_W);
  unsigned int upper_limit_z = (unsigned int)floor(p.vec[2]+half_W);

  unsigned int pair_idx = 0;
  unsigned int write_offset = (sample_idx==0) ? 0 : write_offsets[sample_idx-1];
  for( unsigned int z=lower_limit_z; z<=upper_limit_z; z++ ){
    for( unsigned int y=lower_limit_y; y<=upper_limit_y; y++ ){
      for( unsigned int x=lower_limit_x; x<=upper_limit_x; x++ ){
	vector_td<unsigned int, 3> co; co.vec[0] = x; co.vec[1] = y; co.vec[2] = z;
	tuples_first[write_offset+pair_idx] = co_to_idx<3>(co, matrix_size_os+matrix_size_wrap);
	tuples_last[write_offset+pair_idx] = sample_idx;
	pair_idx++;
      }
    }
  }
}

template <class REAL> __inline__ __device__ void
output_pairs( unsigned int sample_idx, vector_td<REAL,4> &p, vector_td<unsigned int, 4> &matrix_size_os, vector_td<unsigned int, 4> &matrix_size_wrap, REAL half_W, unsigned int *write_offsets, unsigned int *tuples_first, unsigned int *tuples_last )
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
  for( unsigned int w=lower_limit_w; w<=upper_limit_w; w++ ){
    for( unsigned int z=lower_limit_z; z<=upper_limit_z; z++ ){
      for( unsigned int y=lower_limit_y; y<=upper_limit_y; y++ ){
	for( unsigned int x=lower_limit_x; x<=upper_limit_x; x++ ){
	  vector_td<unsigned int, 4> co; co.vec[0] = x; co.vec[1] = y; co.vec[2] = z; co.vec[3] = w;
	  tuples_first[write_offset+pair_idx] = co_to_idx<4>(co, matrix_size_os+matrix_size_wrap);
	  tuples_last[write_offset+pair_idx] = sample_idx;
	  pair_idx++;
	}
      }
    }
  }
}

template<class REAL, unsigned int D> __global__ void
write_pairs_kernel( vector_td<unsigned int, D> matrix_size_os, vector_td<unsigned int, D> matrix_size_wrap, unsigned int num_samples, REAL half_W, vector_td<REAL,D> *traj_positions, unsigned int *write_offsets, unsigned int *tuples_first, unsigned int *tuples_last )
{
  // Get sample idx
  unsigned int sample_idx = blockIdx.x*blockDim.x + threadIdx.x;

  if( sample_idx<num_samples ){

    vector_td<REAL,D> p = traj_positions[sample_idx];
    output_pairs<REAL>( sample_idx, p, matrix_size_os, matrix_size_wrap, half_W, write_offsets, tuples_first, tuples_last );
  }
};

template <class REAL, unsigned int D> void 
write_pairs( vector_td<unsigned int, D> matrix_size_os, vector_td<unsigned int, D> matrix_size_wrap, unsigned int num_samples, REAL W, vector_td<REAL,D> *traj_positions, unsigned int *write_offsets, unsigned int *tuples_first, unsigned int *tuples_last )
{  
  dim3 blockDim(512);
  dim3 gridDim((int)ceil((double)num_samples/(double)blockDim.x));

  REAL half_W = get_half<REAL>()*W;
  write_pairs_kernel<REAL,D><<< gridDim, blockDim >>>
    ( matrix_size_os, matrix_size_wrap, num_samples, half_W, traj_positions, write_offsets, tuples_first, tuples_last );

 CHECK_FOR_CUDA_ERROR();
}
