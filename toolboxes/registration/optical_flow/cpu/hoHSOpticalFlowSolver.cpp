#include "hoHSOpticalFlowSolver.h"
#include "vector_td_utilities.h"

#ifdef USE_OMP
#include <omp.h>
#endif 

namespace Gadgetron{

  // Helpers
  //
  
  template<unsigned long long D> inline bool
  is_border_pixel_for_stride( typename intd<D>::Type stride, typename uintd<D>::Type co, typename uintd<D>::Type dims )
  {
    for( unsigned long long d=0; d<D; d++ ){
      if( stride.vec[d] == -1 ){
	if( co.vec[d] == 0 ){
	  return true;
	}
      }
      else if( stride.vec[d] == 1 ){
	if( co.vec[d] == (dims.vec[d]-1) ){
	  return true;
	}
      }
    }
    return false;
  }
  
  template<unsigned long long i, unsigned long long j> struct Pow
  {
    enum { Value = i*Pow<i,j-1>::Value };
  };
  
  template <unsigned long long i> struct Pow<i,1>
  {
    enum { Value = i };
  };

  //
  // Implementation
  //

  template<class T, unsigned long long D> boost::shared_ptr< hoNDArray<T> >
  hoHSOpticalFlowSolver<T,D>::core_solver( hoNDArray<T> *_gradient_image, hoNDArray<T> *_stencil_image )
  {
    // Sanity checks
    //
  
    if( !_gradient_image ){
      throw std::runtime_error("hoHSOpticalFlowSolver::core_solver(): illegal input gradient image received.");
    }
  
    if( _gradient_image->get_number_of_dimensions() <= D ){
      throw std::runtime_error("hoHSOpticalFlowSolver::core_solver(): number of gradient image dimensions is too small.");
    }
    
    // The dimensions of the displacement field should match the gradient field
    // - when removing the temporal gradient component (replacing D+1 with D)
    //
  
    boost::shared_ptr< std::vector<unsigned long long> > disp_dims = _gradient_image->get_dimensions();
    disp_dims->pop_back(); disp_dims->push_back(D);

    boost::shared_ptr< hoNDArray<T> > displacements_ping(new hoNDArray<T>(disp_dims.get()));
    boost::shared_ptr< hoNDArray<T> > displacements_pong(new hoNDArray<T>(disp_dims.get()));
  
    clear(displacements_ping.get());
    clear(displacements_pong.get());

    // We use "shared memory" to hold the averaged displacements
    boost::shared_ptr< hoNDArray<T> > _shared_mem(new hoNDArray<T>(disp_dims.get()));
    T *shared_mem = _shared_mem->get_data_ptr();
    clear( _shared_mem.get());
   
    typename uintd<D>::Type matrix_size = from_std_vector<unsigned long long,D>( *_gradient_image->get_dimensions() );  
    unsigned long long number_of_elements = prod(matrix_size);
    unsigned long long num_batches = 1;
    
    for( unsigned long long d=D; d<_gradient_image->get_number_of_dimensions()-1; d++ ){
      num_batches *= _gradient_image->get_size(d);
    }
    
    // Get ready...
    //

    unsigned long long iteration_no = 0;
    hoNDArray<T> *ping = displacements_ping.get();
    hoNDArray<T> *pong = displacements_pong.get();

    if( this->output_mode_ >= hoOpticalFlowSolver<T,D>::OUTPUT_VERBOSE ) {
      std::cout << std::endl;
    }

    //
    // Main Jacobi loop
    //

    while(true){
    
      if( this->output_mode_ >= hoOpticalFlowSolver<T,D>::OUTPUT_VERBOSE ) {
	std::cout << "."; std::cout.flush();
      }
    
      // Continuation flag used for early Jacobi termination      
      unsigned long long continue_flag = 0;

      // Number of elements per batch
      const unsigned long long num_elements_per_batch = prod(matrix_size);
      
      // Number of elements per dim
      const unsigned long long num_elements_per_dim = num_elements_per_batch*num_batches;
      
      T *in_disp = ping->get_data_ptr();
      T *out_disp = pong->get_data_ptr();
      T *gradient_image = _gradient_image->get_data_ptr();
      T *stencil_image = (_stencil_image) ? _stencil_image->get_data_ptr() : 0x0;

      //
      // Find the average velocities (shared memory)
      //
      
      for( unsigned long long dim = 0; dim < D; dim++ ){
#ifdef USE_OMP
#pragma omp parallel for
#endif
      for( int idx = 0; idx < num_elements_per_dim; idx++ ){
	  
	  // Index to the shared memory
	  const unsigned long long shared_idx = dim*num_elements_per_dim+idx;
	  
	  // Batch idx (second slowest varying dimension)   
	  const unsigned long long batch_idx = idx/num_elements_per_batch;
	  
	  // Local index to the image (or batch in our terminology)
	  const unsigned long long idx_in_batch = idx-batch_idx*num_elements_per_batch;
    	  
	  if( stencil_image && stencil_image[idx_in_batch] > T(0) )
	    continue;

	  // Local co to the image
	  const typename uintd<D>::Type co = idx_to_co<D>( idx_in_batch, matrix_size );	  
	  const typename intd<D>::Type zeros  = to_vector_td<long long,D>(0);
	  const typename intd<D>::Type ones   = to_vector_td<long long,D>(1);
	  const typename intd<D>::Type threes = to_vector_td<long long,D>(3);
	  
	  const int num_neighbors = Pow<3,D>::Value;
	  T num_contribs = T(0);
      	  
	  for( int i=0; i<num_neighbors; i++ ){
	    
	    // Find the stride of the neighbor {-1, 0, 1}^D
	    const typename intd<D>::Type stride = idx_to_co<D>( i, threes ) - ones;
	    
	    // Verify that the neighbor is not out of bounds (and not the thread itself)
	    if( !is_border_pixel_for_stride<D>( stride, co, matrix_size ) && !(stride==zeros) ){
	  
	      // Compute average of neighbors
	      //
	      
	      const unsigned long long base_offset = dim*num_elements_per_dim + batch_idx*num_elements_per_batch;
	      const unsigned long long neighbor_idx = (unsigned long long) co_to_idx<D>( to_intd(co)+stride, to_intd(matrix_size)) + base_offset;
	  
	      shared_mem[shared_idx] += in_disp[neighbor_idx];
	      num_contribs += T(1);
	    }
	  }
      
	  // Normalize
	  shared_mem[shared_idx] /= num_contribs;       	
	}
      }
      
      //
      // Update displacement field (Jacobi iteration)
      //
      
      const T disp_thresh_sqr = this->limit_*this->limit_;
      
      for( unsigned long long dim = 0; dim < D; dim++ ){
#ifdef USE_OMP
#pragma omp parallel for
#endif
      for( int idx = 0; idx < num_elements_per_dim; idx++ ){

	  // Batch idx (second slowest varying dimension)   
	  const unsigned long long batch_idx = idx/num_elements_per_batch;
	  
	  // Local index to the image (or batch in our terminology)
	  const unsigned long long idx_in_batch = idx-batch_idx*num_elements_per_batch;
    	  
	  if( stencil_image && stencil_image[idx_in_batch] > T(0) )
	    continue;

	  // Index to the shared memory
	  const unsigned long long shared_idx = dim*num_elements_per_dim+idx;

	  T phi = T(0);
	  T norm = T(0);
	  
	  typename reald<T,D>::Type derivatives;
	  
	  // Contributions from the spatial dimensions
	  //
	  
	  for( unsigned long long d=0; d<D; d++ ){
	    derivatives.vec[d] = gradient_image[d*num_elements_per_dim+idx];
	    const unsigned long long shared_idx_d = d*num_elements_per_dim+idx;
	    phi += (shared_mem[shared_idx_d]*derivatives.vec[d]);
	    norm += (derivatives.vec[d]*derivatives.vec[d]);
	  }
	  
	  // Contributions from the temporal dimension
	  //
	  
	  phi += gradient_image[D*num_elements_per_dim+idx];
	  
	  // Normalize
	  //
	  
	  phi /= (alpha_*alpha_+norm);
	  
	  // Form result displacement
	  //
	  
	  T result = shared_mem[shared_idx]-derivatives.vec[dim]*phi;
	  
	  // Clear the "termination" flag if the displacement field has changed above the threshold
	  //
	  
	  T delta = result-in_disp[dim*num_elements_per_dim+idx];
	  if( delta*delta > disp_thresh_sqr )
	    continue_flag = 1;
	  
	  // Output result
	  //
	  
	  out_disp[dim*num_elements_per_dim+idx] = result;
	}
      }
      
      // Swap in/out buffers
      //
      
      hoNDArray<T> *tmp = ping;
      ping = pong;
      pong = tmp;
      
      // Check termination criteria
      //

      if( continue_flag == 0 ){
	if( this->output_mode_ >= hoOpticalFlowSolver<T,D>::OUTPUT_VERBOSE ) {
	  std::cout << std::endl << "Break after " << iteration_no+1 << " iterations" << std::endl;
	}
	break;
      }
    
      if( iteration_no > this->max_num_iterations_per_level_ ) 
	break;    
      
      iteration_no++;
    }
  
    if( ping == displacements_ping.get() )   
      return displacements_ping;
    else
      return displacements_pong;
  }
     
  // 
  // Template instantiation
  //
  
  template class EXPORTCPUREG hoHSOpticalFlowSolver<float,1>;
  template class EXPORTCPUREG hoHSOpticalFlowSolver<float,2>;
  template class EXPORTCPUREG hoHSOpticalFlowSolver<float,3>;
  template class EXPORTCPUREG hoHSOpticalFlowSolver<float,4>;
  
  template class EXPORTCPUREG hoHSOpticalFlowSolver<double,1>;
  template class EXPORTCPUREG hoHSOpticalFlowSolver<double,2>;
  template class EXPORTCPUREG hoHSOpticalFlowSolver<double,3>;
  template class EXPORTCPUREG hoHSOpticalFlowSolver<double,4>;
}
