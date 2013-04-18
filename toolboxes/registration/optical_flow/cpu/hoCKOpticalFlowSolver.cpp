#include "hoCKOpticalFlowSolver.h"
#include "vector_td_utilities.h"

#include <omp.h>

namespace Gadgetron{

  // Helpers
  //
  
  template<unsigned int D> inline bool
  is_border_pixel_for_stride( typename intd<D>::Type stride, typename uintd<D>::Type co, typename uintd<D>::Type dims )
  {
    for( unsigned int d=0; d<D; d++ ){
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

  template<unsigned int i, unsigned int j> struct Pow
  {
    enum { Value = i*Pow<i,j-1>::Value };
  };
  
  template <unsigned int i> struct Pow<i,1>
  {
    enum { Value = i };
  };
  
  //
  // Implementation
  //

  template<class REAL, unsigned int D> boost::shared_ptr< hoNDArray<REAL> >
  hoCKOpticalFlowSolver<REAL,D>::core_solver( hoNDArray<REAL> *_gradient_image, hoNDArray<REAL> *_stencil_image )
  {
    // Sanity checks
    //
  
    if( !_gradient_image ){
      BOOST_THROW_EXCEPTION(runtime_error("hoCKOpticalFlowSolver::core_solver(): illegal input gradient image received."));
    }

    if( _gradient_image->get_number_of_dimensions() <= D ){
      BOOST_THROW_EXCEPTION(runtime_error("hoCKOpticalFlowSolver::core_solver(): number of gradient image dimensions is too small."));
    }
  
    // The dimensions of the displacement field should match the gradient field
    //
  
    boost::shared_ptr< std::vector<unsigned int> > disp_dims = _gradient_image->get_dimensions();
    boost::shared_ptr< hoNDArray<REAL> > displacements_ping( new hoNDArray<REAL>(disp_dims.get()) );
    boost::shared_ptr< hoNDArray<REAL> > displacements_pong( new hoNDArray<REAL>(disp_dims.get()) );
    clear(displacements_ping.get());
    clear(displacements_pong.get());
    
    // We use "shared memory" to hold the averaged displacements
    boost::shared_ptr< hoNDArray<REAL> > _shared_mem(new hoNDArray<REAL>(disp_dims.get()));
    REAL *shared_mem = _shared_mem->get_data_ptr();
    clear( _shared_mem.get());

    typename uintd<D>::Type matrix_size = from_std_vector<unsigned int,D>( *disp_dims );  
    unsigned int number_of_elements = prod(matrix_size);
    unsigned int num_batches = 1;

    for( unsigned int d=D; d<_gradient_image->get_number_of_dimensions()-1; d++ ){
      num_batches *= _gradient_image->get_size(d);
    }
  
    // Get ready
    // 

    unsigned int iteration_no = 0;
    hoNDArray<REAL> *ping = displacements_ping.get();
    hoNDArray<REAL> *pong = displacements_pong.get(); 

    if( this->output_mode_ >= hoOpticalFlowSolver<REAL,D>::OUTPUT_VERBOSE ) {
      std::cout << std::endl;
    }

    //
    // Main Jacobi loop
    //
    
    while(true){
    
      if( this->output_mode_ >= hoOpticalFlowSolver<REAL,D>::OUTPUT_VERBOSE ) {
	std::cout << "."; std::cout.flush();
      }
    
      // Continuation flag used for early Jacobi termination
      unsigned int continue_flag = 0;

      // Number of elements per batch
      const unsigned int num_elements_per_batch = prod(matrix_size);
  
      // Number of elements per dim
      const unsigned int num_elements_per_dim = num_elements_per_batch*num_batches;

      REAL *in_disp = ping->get_data_ptr();
      REAL *out_disp = pong->get_data_ptr();
      REAL *gradient_image = _gradient_image->get_data_ptr();
      REAL *stencil_image = (_stencil_image) ? _stencil_image->get_data_ptr() : 0x0;

      //
      // Find the average velocities (shared memory)
      //
      
      for( unsigned int dim = 0; dim < D+1; dim++ ){
#pragma omp parallel for
	for( unsigned int idx = 0; idx < num_elements_per_dim; idx++ ){
	  	  
	  // Index to the shared memory
	  const unsigned int shared_idx = dim*num_elements_per_dim+idx;
	  
	  // Batch idx (second slowest varying dimension)   
	  const unsigned int batch_idx = idx/num_elements_per_batch;
	  
	  // Local index to the image (or batch in our terminology)
	  const unsigned int idx_in_batch = idx-batch_idx*num_elements_per_batch;
    	  
	  if( stencil_image && stencil_image[idx_in_batch] > REAL(0) )
	    continue;

	  // Local co to the image
	  const typename uintd<D>::Type co = idx_to_co<D>( idx_in_batch, matrix_size );    
	  const typename intd<D>::Type zeros  = to_vector_td<int,D>(0);
	  const typename intd<D>::Type ones   = to_vector_td<int,D>(1);
	  const typename intd<D>::Type threes = to_vector_td<int,D>(3);
	  
	  const int num_neighbors = Pow<3,D>::Value;
	  REAL num_contribs = REAL(0);
	  
	  shared_mem[shared_idx] = REAL(0);
	  
	  // Compute average of neighbors
	  //
	  
	  for( int i=0; i<num_neighbors; i++ ){
	    
	    // Find the stride of the neighbor {-1, 0, 1}^D
	    const typename intd<D>::Type stride = idx_to_co<D>( i, threes ) - ones;
	    
	    unsigned int neighbor_idx;
	    
	    const unsigned int base_offset = dim*num_elements_per_dim + batch_idx*num_elements_per_batch;
	    
	    // Verify that the neighbor is not out of bounds (and not the thread itself)
	    if( !is_border_pixel_for_stride<D>( stride, co, matrix_size ) && !(stride==zeros) ){	
	      neighbor_idx = (unsigned int) co_to_idx<D>( to_intd(co)+stride, to_intd(matrix_size)) + base_offset;
	    }
	    else{
	      neighbor_idx = idx_in_batch + base_offset;
	    }
	    
	    shared_mem[shared_idx] += in_disp[neighbor_idx];
	    num_contribs += REAL(1);
	  }
	  
	  // Normalize
	  shared_mem[shared_idx] /= num_contribs;
	}
      }
      
      //
      // Update displacement field (Jacobi iteration)
      //
      
      const REAL disp_thresh_sqr = this->limit_*this->limit_;

      for( unsigned int dim = 0; dim < D+1; dim++ ){
#pragma omp parallel for
	for( unsigned int idx = 0; idx < num_elements_per_dim; idx++ ){
	  
	  // Index to the shared memory
	  const unsigned int shared_idx = dim*num_elements_per_dim+idx;
	  
	  // Batch idx (second slowest varying dimension)   
	  const unsigned int batch_idx = idx/num_elements_per_batch;
	  
	  // Local index to the image (or batch in our terminology)
	  const unsigned int idx_in_batch = idx-batch_idx*num_elements_per_batch;
    	  
	  if( stencil_image && stencil_image[idx_in_batch] > REAL(0) )
	    continue;

	  REAL phi = REAL(0);
	  REAL norm = REAL(0);
	  
	  typename reald<REAL,D>::Type derivatives;
	  
	  // Contributions from the spatial dimensions
	  //
	  
	  for( unsigned int d=0; d<D; d++ ){
	    derivatives.vec[d] = gradient_image[d*num_elements_per_dim+idx];
	    const unsigned int shared_idx_d = d*num_elements_per_dim+idx;
	    phi += (shared_mem[shared_idx_d]*derivatives.vec[d]);
	    norm += (derivatives.vec[d]*derivatives.vec[d]);
	  }
	  
	  // Contributions from the temporal dimension
	  //
	  
	  phi += gradient_image[D*num_elements_per_dim+idx];
	  
	  // Contribution from the intensity attentuation estimation
	  //
	  
	  phi -= shared_mem[D*num_elements_per_dim+idx];
	  
	  // Normalize
	  //
	  
	  phi /= ((alpha_/beta_)*(alpha_/beta_)+alpha_*alpha_+norm);
	  
	  // Form result displacement
	  //
	  
	  REAL result;
	  
	  if( dim<D )
	    result = shared_mem[shared_idx]-derivatives.vec[dim]*phi;
	  else
	    result = shared_mem[D*num_elements_per_dim+idx]+(alpha_/beta_)*(alpha_/beta_)*phi;
	  
	  // Clear the "termination" flag if the displacement field has changed above the threshold
	  //
	  
	  REAL delta = result-in_disp[dim*num_elements_per_dim+idx];
	  if( dim < D && delta*delta > disp_thresh_sqr )
	    continue_flag = 1;
	  
	  // Output result
	  //
	  
	  out_disp[dim*num_elements_per_dim+idx] = result;
	}
      }
      
      // Swap in/out buffers
      //
      
      hoNDArray<REAL> *tmp = ping;
      ping = pong;
      pong = tmp;

      // Check termination criteria
      //
      
      if( continue_flag == 0 ){
	if( this->output_mode_ >= hoOpticalFlowSolver<REAL,D>::OUTPUT_VERBOSE ) {
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

  template class EXPORTCPUREG hoCKOpticalFlowSolver<float,1>;
  template class EXPORTCPUREG hoCKOpticalFlowSolver<float,2>;
  template class EXPORTCPUREG hoCKOpticalFlowSolver<float,3>;
  template class EXPORTCPUREG hoCKOpticalFlowSolver<float,4>;

  template class EXPORTCPUREG hoCKOpticalFlowSolver<double,1>;
  template class EXPORTCPUREG hoCKOpticalFlowSolver<double,2>;
  template class EXPORTCPUREG hoCKOpticalFlowSolver<double,3>;
  template class EXPORTCPUREG hoCKOpticalFlowSolver<double,4>;  
}
